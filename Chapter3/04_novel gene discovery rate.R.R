# ============================================================
# Heaps pan-genome accumulation curve & novel gene discovery rate
# Based on Panaroo gene_presence_absence_roary.csv
# Four species: Borreliella burgdorferi / garinii / afzelii / bavariensis
# ============================================================

library(tidyverse)
library(readxl)
library(ggplot2)

# --- 1. Paths (replace with your local paths) ---
pav_path <- "path/to/panaroo_four_species_complete/gene_presence_absence_roary.csv"
s1_path  <- "path/to/Supplementary_Table_S1.xls"
out_dir  <- "path/to/output"

# --- 2. Color palette & species order ---
sp_colors <- c(
  "Borreliella burgdorferi"   = "#1F618D",
  "Borreliella garinii"       = "#E74C3C",
  "Borreliella afzelii"       = "#1ABC9C",
  "Borreliella bavariensis"   = "#5DADE2"
)
sp_levels <- names(sp_colors)

# --- 3. Load S1 metadata ---
s1 <- read_excel(s1_path) %>%
  mutate(
    Genome_ID = sub("_genomic$", "", Genome_ID),
    species = case_when(
      str_detect(species, "burgdorferi") ~ "Borreliella burgdorferi",
      str_detect(species, "garinii")     ~ "Borreliella garinii",
      str_detect(species, "afzelii")     ~ "Borreliella afzelii",
      str_detect(species, "bavariensis") ~ "Borreliella bavariensis",
      TRUE ~ species
    )
  ) %>%
  filter(species %in% sp_levels)

# --- 4. Load PAV matrix ---
pav <- read_csv(pav_path, show_col_types = FALSE)

meta_cols <- c("Gene", "Non-unique Gene name", "Annotation", "No. isolates",
               "No. sequences", "Avg group size", "Genome Fragment", 
               "Order within Fragment", "Accessory Fragment", 
               "Accessory Order with Fragment", "QC", "Min group size nuc",
               "Max group size nuc", "Avg group size nuc")
genome_cols <- setdiff(names(pav), meta_cols)

genome_map <- tibble(
  col_name = genome_cols,
  Genome_ID = sub("_genomic$", "", genome_cols)
)

valid_genomes <- genome_map %>%
  inner_join(s1 %>% select(Genome_ID, species), by = "Genome_ID")

n_per_species <- valid_genomes %>% count(species) %>% deframe()

# Binary PAV matrix (1 = present, 0 = absent)
pav_mat <- pav %>%
  select(Gene, all_of(valid_genomes$col_name)) %>%
  column_to_rownames("Gene") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0))) %>%
  as.matrix()
colnames(pav_mat) <- valid_genomes$Genome_ID

# Fix NA values if any
if (any(is.na(pav_mat))) {
  pav_mat[is.na(pav_mat)] <- 0
}

# --- 5. Heaps curve: random permutation (n = 100) ---
calc_heaps <- function(species_name, n_perm = 100) {
  gids <- valid_genomes %>% filter(species == species_name) %>% pull(Genome_ID)
  n <- length(gids)
  if (n < 5) return(NULL)
  
  cum_matrix <- matrix(0, nrow = n_perm, ncol = n)
  for (p in 1:n_perm) {
    perm <- sample(gids)
    pool <- character(0)
    for (i in 1:n) {
      present <- rownames(pav_mat)[pav_mat[, perm[i]] == 1]
      pool <- union(pool, present)
      cum_matrix[p, i] <- length(pool)
    }
  }
  tibble(
    species = species_name,
    n_genomes = 1:n,
    mean_cum = colMeans(cum_matrix),
    sd_cum = apply(cum_matrix, 2, sd)
  )
}

heaps_list <- map(sp_levels, calc_heaps, n_perm = 100)
heaps_df <- bind_rows(heaps_list) %>%
  mutate(species = factor(species, levels = sp_levels))

# --- 6. Power-law fitting: n = k * N^gamma ---
fit_heaps <- function(df) {
  fit <- nls(mean_cum ~ k * n_genomes^gamma, data = df, 
             start = list(k = max(df$mean_cum) / 2, gamma = 0.5),
             control = nls.control(maxiter = 500))
  ci <- confint2(fit)
  tibble(
    species = unique(df$species),
    gamma = coef(fit)["gamma"],
    gamma_low = ci["gamma", 1],
    gamma_high = ci["gamma", 2]
  )
}

heaps_fits <- heaps_df %>%
  group_by(species) %>%
  group_map(~ fit_heaps(.x), .keep = TRUE) %>%
  bind_rows() %>%
  mutate(species = factor(species, levels = sp_levels))

# --- 7. Plot: Heaps accumulation curve ---
legend_labels <- map_chr(sp_levels, function(sp) {
  n <- n_per_species[sp]
  g <- heaps_fits$gamma[heaps_fits$species == sp]
  sp_short <- gsub("Borreliella ", "B. ", sp)
  sprintf("italic('%s') ~ ' (n=%d, γ=%.3f)'", sp_short, n, g)
})
names(legend_labels) <- sp_levels
legend_expr <- map(legend_labels, ~ parse(text = .x))
names(legend_expr) <- sp_levels

p_heaps <- ggplot(heaps_df, aes(x = n_genomes, y = mean_cum, color = species, fill = species)) +
  geom_ribbon(aes(ymin = mean_cum - sd_cum, ymax = mean_cum + sd_cum),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = sp_colors, labels = legend_expr, name = "Species") +
  scale_fill_manual(values = sp_colors, labels = legend_expr, name = "Species") +
  labs(x = "Number of Genomes (N)", y = "Cumulative Gene Number (n)") +
  coord_cartesian(ylim = c(600, 1600)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = c(0.72, 0.32),
    legend.background = element_blank(),
    legend.text = element_text(size = 9.5),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.key = element_blank(),
    legend.key.height = unit(0.9, "lines")
  )

print(p_heaps)

# Optional save
# ggsave(file.path(out_dir, "heaps_accumulation_curve.png"), p_heaps,
#        width = 7, height = 6, dpi = 300, bg = "white")

# --- 8. Novel gene discovery rate (percentage) ---
calc_novel <- function(species_name, n_perm = 100) {
  gids <- valid_genomes %>% filter(species == species_name) %>% pull(Genome_ID)
  n <- length(gids)
  if (n < 5) return(NULL)
  
  mat <- matrix(0, nrow = n_perm, ncol = n)
  for (p in 1:n_perm) {
    perm <- sample(gids)
    pool <- character(0)
    for (i in 1:n) {
      present <- rownames(pav_mat)[pav_mat[, perm[i]] == 1]
      new_genes <- setdiff(present, pool)
      mat[p, i] <- length(new_genes)
      pool <- union(pool, present)
    }
  }
  tibble(
    species = species_name,
    n_genomes = 1:n,
    mean_novel = colMeans(mat),
    sd_novel = apply(mat, 2, sd)
  )
}

novel_list <- map(sp_levels, calc_novel, n_perm = 100)
novel_df <- bind_rows(novel_list)

# Mean genome size per species (for percentage conversion)
mean_genome_size <- tibble(
  Genome_ID = colnames(pav_mat),
  gene_count = colSums(pav_mat)
) %>%
  inner_join(valid_genomes, by = "Genome_ID") %>%
  group_by(species) %>%
  summarise(mean_size = mean(gene_count), .groups = "drop")

novel_df <- novel_df %>%
  mutate(species = as.character(species)) %>%
  left_join(
    mean_genome_size %>% mutate(species = as.character(species)),
    by = "species"
  ) %>%
  mutate(
    species = factor(species, levels = sp_levels),
    novel_pct = mean_novel / mean_size * 100,
    sd_pct = sd_novel / mean_size * 100
  )

# End-point data for annotation
novel_end <- novel_df %>%
  group_by(species) %>%
  slice_max(n_genomes, n = 1) %>%
  ungroup()

# Plot: Novel gene discovery rate (%)
legend_labels_b <- map_chr(sp_levels, function(sp) {
  n <- n_per_species[sp]
  sp_short <- gsub("Borreliella ", "B. ", sp)
  sprintf("italic('%s') ~ ' (n=%d)'", sp_short, n)
})
names(legend_labels_b) <- sp_levels
legend_expr_b <- map(legend_labels_b, ~ parse(text = .x))
names(legend_expr_b) <- sp_levels

p_novel <- ggplot(novel_df, aes(x = n_genomes, y = novel_pct, color = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) +
  annotate("text", x = 18, y = 1.12, label = "Saturation threshold (1%)", 
           color = "gray40", size = 3.5) +
  geom_segment(data = novel_end,
               aes(x = n_genomes, xend = n_genomes, y = 0, yend = 5, color = species),
               linetype = "dotted", linewidth = 0.8, alpha = 0.9) +
  geom_line(linewidth = 1.0) +
  geom_text(data = novel_end,
            aes(x = n_genomes + 1.5, y = 0.8,
                label = sprintf("+%.2f%%", novel_pct), color = species),
            hjust = 0, vjust = 0.5, size = 4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = sp_colors, labels = legend_expr_b, name = "Species") +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(x = "Number of Genomes (N)", y = "New Gene Discovery Rate (%)") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = c(0.78, 0.82),
    legend.background = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.key = element_blank(),
    legend.key.height = unit(0.9, "lines")
  )

# Uncomment to print/save
# print(p_novel)
# ggsave(file.path(out_dir, "novel_gene_discovery_rate.png"), p_novel,
#        width = 7, height = 6, dpi = 300, bg = "white")

# --- 9. Summary statistics ---
cat("\n=== Heaps power-law parameters ===\n")
heaps_fits %>%
  mutate(label = sprintf("%s: γ=%.3f (95%% CI: %.3f–%.3f)",
                         species, gamma, gamma_low, gamma_high)) %>%
  pull(label) %>% cat(sep = "\n")

cat("\n=== Final cumulative gene counts ===\n")
heaps_df %>%
  group_by(species) %>%
  slice_max(n_genomes, n = 1) %>%
  select(species, n_genomes, mean_cum) %>%
  mutate(mean_cum = round(mean_cum)) %>%
  print()

cat("\n=== End-point novel gene rate (%) ===\n")
novel_end %>%
  select(species, n_genomes, novel_pct) %>%
  mutate(label = sprintf("%s (N=%d): %.2f%%", species, n_genomes, novel_pct)) %>%
  pull(label) %>% cat(sep = "\n")

cat("\n=== Saturation status ===\n")
novel_end %>%
  mutate(status = ifelse(novel_pct < 1, "Saturated (<1%)", "Unsaturated (>1%)")) %>%
  select(species, novel_pct, status) %>%
  print()