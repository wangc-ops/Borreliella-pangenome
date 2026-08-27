# ============================================================
# CDS copy number stratification by pangenome tier (species-level)
# Stacked bar: actual CDS copies per tier
# ============================================================

library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)

# --- 1. Paths (replace with your local paths) ---
pav_file     <- "path/to/panaroo_four_species_complete/gene_presence_absence_roary.csv"
species_file <- "path/to/Supplementary_Table_S1.xls"
out_dir      <- "path/to/output"

# --- 2. Parameters ---
sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_short_levels <- c("B. burgdorferi", "B. garinii", "B. afzelii", "B. bavariensis")

tier_colors <- c(
  "Core"      = "#D9534F",
  "Soft-core" = "#F0AD4E",
  "Shell"     = "#5B9BD5",
  "Cloud"     = "#A0CBE8"
)

# --- 3. Read PAV ---
pav <- fread(pav_file, header = TRUE)

meta_cols <- c("Gene", "Non-unique Gene name", "Annotation", "No. isolates",
               "No. sequences", "Avg group size", "Genome Fragment", 
               "Order within Fragment", "Accessory Fragment", 
               "Accessory Order with Fragment", "QC", "Min group size nuc",
               "Max group size nuc", "Avg group size nuc")
strain_cols <- setdiff(names(pav), meta_cols)

pav[, (strain_cols) := lapply(.SD, function(x) as.integer(!is.na(x) & x > 0)), .SDcols = strain_cols]

# --- 4. Read S1 ---
s1 <- as.data.table(read_excel(species_file))
s1[, Genome_ID := sub("_genomic$", "", Genome_ID)]
s1[, species := fcase(
  grepl("burgdorferi", species), "Borreliella burgdorferi",
  grepl("garinii", species),     "Borreliella garinii",
  grepl("afzelii", species),     "Borreliella afzelii",
  grepl("bavariensis", species), "Borreliella bavariensis",
  default = species
)]
s1 <- s1[species %in% sp_levels]

# --- 5. Match ---
genome_map <- data.table(raw_col = strain_cols, Genome_ID = sub("_genomic$", "", strain_cols))
matched_dt <- genome_map[s1, on = "Genome_ID", nomatch = 0L]
pav_sub <- pav[, c("Gene", matched_dt$raw_col), with = FALSE]

# --- 6. Global tier ---
global_freq <- rowMeans(pav_sub[, -"Gene"])
gene_meta <- data.table(
  Gene = pav_sub$Gene,
  global_freq = global_freq,
  tier = fifelse(global_freq >= 0.99, "Core",
                 fifelse(global_freq >= 0.95, "Soft-core",
                         fifelse(global_freq >= 0.15, "Shell", "Cloud")))
)
gene_meta[, tier := factor(tier, levels = names(tier_colors))]

# --- 7. CDS copy number per tier per species ---
species_tier_cds <- map_dfr(sp_levels, function(sp) {
  gids <- matched_dt[species == sp, raw_col]
  sp_mat <- pav_sub[, c("Gene", gids), with = FALSE]
  cds_counts <- rowSums(sp_mat[, -"Gene"])
  
  data.table(Gene = sp_mat$Gene, cds = cds_counts) %>%
    left_join(gene_meta, by = "Gene") %>%
    group_by(tier) %>%
    summarise(total_cds = sum(cds), .groups = "drop") %>%
    mutate(
      species = sp,
      total_all = sum(total_cds),
      prop = total_cds / total_all * 100,
      species_short = factor(gsub("Borreliella ", "B. ", sp), levels = sp_short_levels)
    )
})

# --- 8. Plot ---
totals_cds <- species_tier_cds %>%
  distinct(species_short, total_all) %>%
  mutate(label_top = sprintf("Sigma == %d", total_all))

p_cds <- ggplot(species_tier_cds, aes(x = species_short, y = prop, fill = tier)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", prop)), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "black") +
  geom_text(data = totals_cds, aes(x = species_short, y = 102, label = label_top),
            inherit.aes = FALSE, size = 3.5, fontface = "bold", vjust = 0, parse = TRUE) +
  scale_fill_manual(values = tier_colors, name = "Gene category") +
  scale_x_discrete(expand = c(0.15, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 108), breaks = seq(0, 100, 25)) +
  labs(x = "Species", y = "Proportion of CDS copies (%)") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(face = "italic", size = 11, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9.5),
    legend.key = element_rect(color = NA),
    legend.key.size = unit(0.8, "lines"),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5)
  )

print(p_cds)

# Optional save
# ggsave(file.path(out_dir, "cds_copy_tier_composition.png"), p_cds,
#        width = 6.5, height = 5.5, dpi = 300, bg = "white")

# Optional supplementary table
# fwrite(species_tier_cds, file.path(out_dir, "cds_copy_tier_data.csv"))