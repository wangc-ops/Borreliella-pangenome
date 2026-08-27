# ============================================================
# fullsample_pan_genome_analysis.R
# Full-sample Heaps curve and novel gene discovery rate
# for four Borreliella species (93/75/39/34 genomes)
# Upload version — replace YOUR_PATH_TO before running
# ============================================================
library(tidyverse)
library(readxl)
library(ggplot2)

# --- user-defined paths ---
pav_path <- "YOUR_PATH_TO/gene_presence_absence_roary.csv"
s1_path  <- "YOUR_PATH_TO/Supplementary_Table_S1.xls"
out_dir  <- "YOUR_OUTPUT_DIRECTORY"

sp_colors <- c(
  "Borreliella burgdorferi"   = "#1F618D",
  "Borreliella garinii"       = "#E74C3C",
  "Borreliella afzelii"       = "#1ABC9C",
  "Borreliella bavariensis"   = "#5DADE2"
)
sp_levels <- names(sp_colors)
N_PERM    <- 100

# --- read S1 ---
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
  ) %>% filter(species %in% sp_levels)

# --- read PAV matrix ---
pav <- read_csv(pav_path, show_col_types = FALSE)
meta_cols <- c("Gene", "Non-unique Gene name", "Annotation", "No. isolates",
               "No. sequences", "Avg group size", "Genome Fragment", 
               "Order within Fragment", "Accessory Fragment", 
               "Accessory Order with Fragment", "QC", "Min group size nuc",
               "Max group size nuc", "Avg group size nuc")
genome_cols <- setdiff(names(pav), meta_cols)

genome_map <- tibble(col_name = genome_cols, Genome_ID = sub("_genomic$", "", genome_cols))
valid_genomes <- genome_map %>% inner_join(s1 %>% select(Genome_ID, species), by = "Genome_ID")
n_per_species <- valid_genomes %>% count(species) %>% deframe()

pav_mat <- pav %>%
  select(Gene, all_of(valid_genomes$col_name)) %>%
  column_to_rownames("Gene") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0))) %>%
  as.matrix()
colnames(pav_mat) <- valid_genomes$Genome_ID
if (any(is.na(pav_mat))) pav_mat[is.na(pav_mat)] <- 0

# --- mean genome size ---
mean_genome_size <- tibble(
  Genome_ID = colnames(pav_mat),
  gene_count = colSums(pav_mat)
) %>%
  inner_join(valid_genomes, by = "Genome_ID") %>%
  group_by(species) %>%
  summarise(mean_size = mean(gene_count), .groups = "drop")

# ============================================================
# Panel 1: Full-sample Heaps curve
# ============================================================
calc_heaps_full <- function(species_name, n_perm = N_PERM) {
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
  tibble(species = species_name, n_genomes = 1:n,
         mean_cum = colMeans(cum_matrix), sd_cum = apply(cum_matrix, 2, sd))
}

heaps_list <- map(sp_levels, calc_heaps_full)
heaps_df   <- bind_rows(heaps_list) %>% mutate(species = factor(species, levels = sp_levels))

fit_heaps <- function(df) {
  fit <- nls(mean_cum ~ k * n_genomes^gamma, data = df, 
             start = list(k = max(df$mean_cum)/2, gamma = 0.5),
             control = nls.control(maxiter = 500))
  ci <- tryCatch(confint(fit), error = function(e) matrix(NA, nrow=2, ncol=2))
  tibble(species = unique(df$species), gamma = coef(fit)["gamma"],
         gamma_low = ci["gamma", 1], gamma_high = ci["gamma", 2], k = coef(fit)["k"])
}

heaps_fits <- heaps_df %>% group_by(species) %>% 
  group_map(~ fit_heaps(.x), .keep = TRUE) %>% bind_rows() %>%
  mutate(species = factor(species, levels = sp_levels))

cat("=== Full-sample power-law fit ===\n")
print(heaps_fits)

legend_labels_a <- map_chr(sp_levels, function(sp) {
  n <- n_per_species[sp]
  g <- heaps_fits$gamma[heaps_fits$species == sp]
  sp_short <- gsub("Borreliella ", "B. ", sp)
  sprintf("italic(\'%s\') ~ \' (n=%d, \u03B3=%.3f)\'", sp_short, n, g)
})
names(legend_labels_a) <- sp_levels
legend_expr_a <- map(legend_labels_a, ~ parse(text = .x))
names(legend_expr_a) <- sp_levels

heaps_plot <- ggplot(heaps_df, aes(x = n_genomes, y = mean_cum, color = species, fill = species)) +
  geom_ribbon(aes(ymin = mean_cum - sd_cum, ymax = mean_cum + sd_cum), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = sp_colors, labels = legend_expr_a, name = "Species") +
  scale_fill_manual(values = sp_colors, labels = legend_expr_a, name = "Species") +
  labs(x = "Number of Genomes (N)", y = "Cumulative Gene Number (n)") +
  coord_cartesian(ylim = c(600, 1600)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        legend.position = c(0.72, 0.32), legend.background = element_blank(),
        legend.text = element_text(size = 9.5), legend.title = element_text(size = 10.5, face = "bold"),
        legend.key = element_blank(), legend.key.height = unit(0.9, "lines"))

print(heaps_plot)
ggsave(file.path(out_dir, "heaps_curve_fullsample.png"), heaps_plot, 
       width = 7, height = 6, dpi = 300, bg = "white")

# ============================================================
# Panel 2: Full-sample novel gene discovery rate (Algorithm A + full vlines)
# ============================================================
calc_novel_full <- function(species_name, n_perm = N_PERM) {
  gids <- valid_genomes %>% filter(species == species_name) %>% pull(Genome_ID)
  n <- length(gids)
  if (n < 5) return(NULL)
  
  sp_size <- mean_genome_size %>% filter(species == species_name) %>% pull(mean_size)
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
  tibble(species = species_name, n_genomes = 1:n,
         mean_novel = colMeans(mat), sd_novel = apply(mat, 2, sd), mean_size = sp_size)
}

novel_list <- map(sp_levels, calc_novel_full)
novel_df <- bind_rows(novel_list) %>%
  mutate(species = factor(species, levels = sp_levels),
         novel_pct = mean_novel / mean_size * 100,
         sd_pct = sd_novel / mean_size * 100)

novel_end <- novel_df %>% group_by(species) %>% slice_max(n_genomes, n = 1) %>% ungroup() %>%
  mutate(
    label = sprintf("+%.2f%%", novel_pct),
    y_label = case_when(
      species == "Borreliella garinii"       ~ 0.85,
      species == "Borreliella afzelii"     ~ 0.60,
      species == "Borreliella burgdorferi" ~ 0.40,
      species == "Borreliella bavariensis" ~ 0.20,
      TRUE ~ 0.8
    ),
    x_label = n_genomes + 1.5
  )

cat("\n=== Full-sample end-point novel gene rate ===\n")
print(novel_end %>% select(species, n_genomes, novel_pct, label))

legend_labels_b <- map_chr(sp_levels, function(sp) {
  n <- n_per_species[sp]
  sp_short <- gsub("Borreliella ", "B. ", sp)
  sprintf("italic(\'%s\') ~ \' (n=%d)\'", sp_short, n)
})
names(legend_labels_b) <- sp_levels
legend_expr_b <- map(legend_labels_b, ~ parse(text = .x))
names(legend_expr_b) <- sp_levels

novel_plot <- ggplot(novel_df, aes(x = n_genomes, y = novel_pct, color = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) +
  annotate("text", x = 95, y = 1.15, label = "Saturation threshold (1%)", 
           color = "gray40", size = 3.5, hjust = 1) +
  geom_segment(data = novel_end,
               aes(x = n_genomes, xend = n_genomes, y = 0, yend = 5, color = species),
               linetype = "dotted", linewidth = 0.8, alpha = 0.9) +
  geom_line(linewidth = 1.0) +
  geom_text(data = novel_end,
            aes(x = x_label, y = y_label, label = label, color = species),
            hjust = 0, vjust = 0.5, size = 3.8, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = sp_colors, labels = legend_expr_b, name = "Species") +
  labs(x = "Number of Genomes (N)", y = "New Gene Discovery Rate (%)") +
  coord_cartesian(xlim = c(0, 105), ylim = c(0, 5.2)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 5, 1), expand = c(0, 0)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        legend.position = c(0.78, 0.82), legend.background = element_blank(),
        legend.text = element_text(size = 10), legend.title = element_text(size = 11, face = "bold"),
        legend.key = element_blank(), legend.key.height = unit(0.9, "lines"))

print(novel_plot)
ggsave(file.path(out_dir, "novel_gene_rate_fullsample.png"), novel_plot, 
       width = 7, height = 5, dpi = 300, bg = "white")

cat("\n=== Saturation status (threshold = 1%) ===\n")
novel_end %>% mutate(status = ifelse(novel_pct < 1, "Saturated (<1%)", "Unsaturated (>1%)"),
                     label = sprintf("%s: %.2f%% (%s)", species, novel_pct, status)) %>%
  pull(label) %>% cat(sep = "\n")