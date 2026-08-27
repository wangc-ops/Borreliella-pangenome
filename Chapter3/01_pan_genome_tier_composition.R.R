# ============================================================
# pan_genome_tier_composition.R
n# Pan-genome tier pie chart and gene family frequency distribution
# for four Borreliella species (241 genomes)
# Upload version — replace YOUR_PATH_TO before running
# ============================================================
library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)

# --- user-defined paths ---
pav_file     <- "YOUR_PATH_TO/gene_presence_absence_roary.csv"
species_file <- "YOUR_PATH_TO/Supplementary_Table_S1.xls"
out_dir      <- "YOUR_OUTPUT_DIRECTORY"

# --- parameters ---
sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_short_levels <- c("B. burgdorferi", "B. garinii", "B. afzelii", "B. bavariensis")

tier_colors <- c(
  "Core"      = "#D9534F",
  "Soft-core" = "#F0AD4E",
  "Shell"     = "#5B9BD5",
  "Cloud"     = "#A0CBE8"
)

# --- read PAV ---
pav <- fread(pav_file, header = TRUE)

meta_cols <- c("Gene", "Non-unique Gene name", "Annotation", "No. isolates",
               "No. sequences", "Avg group size", "Genome Fragment", 
               "Order within Fragment", "Accessory Fragment", 
               "Accessory Order with Fragment", "QC", "Min group size nuc",
               "Max group size nuc", "Avg group size nuc")
strain_cols <- setdiff(names(pav), meta_cols)

pav[, (strain_cols) := lapply(.SD, function(x) as.integer(!is.na(x) & x > 0)), .SDcols = strain_cols]

# --- read S1 ---
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

# --- match ---
genome_map <- data.table(raw_col = strain_cols, Genome_ID = sub("_genomic$", "", strain_cols))
matched_dt <- genome_map[s1, on = "Genome_ID", nomatch = 0L]
pav_sub <- pav[, c("Gene", matched_dt$raw_col), with = FALSE]

n_total <- ncol(pav_sub) - 1

# --- global tier ---
global_freq <- rowMeans(pav_sub[, -"Gene"])
gene_meta <- data.table(
  Gene = pav_sub$Gene,
  global_freq = global_freq,
  n_isolates = rowSums(pav_sub[, -"Gene"]),
  tier = fifelse(global_freq >= 0.99, "Core",
                 fifelse(global_freq >= 0.95, "Soft-core",
                         fifelse(global_freq >= 0.15, "Shell", "Cloud")))
)
gene_meta[, tier := factor(tier, levels = names(tier_colors))]

# ============================================================
# Panel 1: Pan-genome tier composition (pie chart)
# ============================================================

tier_summary <- gene_meta %>%
  count(tier) %>%
  mutate(
    prop = n / sum(n) * 100,
    label = sprintf("%s\nn = %d (%.1f%%)", tier, n, prop)
  )

pie_plot <- ggplot(tier_summary, aes(x = "", y = n, fill = tier)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = tier_colors, name = NULL) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.2, color = "black") +
  labs(title = NULL) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key = element_rect(color = NA),
    legend.key.size = unit(0.8, "lines"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

print(pie_plot)
ggsave(file.path(out_dir, "pan_genome_tier_pie.png"), pie_plot,
       width = 5, height = 4, dpi = 300, bg = "white")

# ============================================================
# Panel 2: Gene family frequency distribution (histogram)
# ============================================================

gene_meta[, freq_bin := fcase(
  n_isolates == n_total, "100%",
  n_isolates >= 0.95 * n_total, "95\u201399%",
  n_isolates >= 0.75 * n_total, "75\u201394%",
  n_isolates >= 0.50 * n_total, "50\u201374%",
  n_isolates >= 0.25 * n_total, "25\u201349%",
  n_isolates >= 0.15 * n_total, "15\u201324%",
  n_isolates >= 0.05 * n_total, "5\u201314%",
  n_isolates >= 1, "1\u20134%",
  default = "0%"
)]

freq_order <- c("100%", "95\u201399%", "75\u201394%", "50\u201374%", "25\u201349%", 
                "15\u201324%", "5\u201314%", "1\u20134%", "0%")
gene_meta[, freq_bin := factor(freq_bin, levels = freq_order)]

freq_summary <- gene_meta %>%
  count(freq_bin) %>%
  mutate(prop = n / sum(n) * 100)

freq_plot <- ggplot(freq_summary, aes(x = freq_bin, y = n, fill = freq_bin)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, prop)), vjust = -0.3, size = 3, color = "black") +
  scale_fill_manual(
    values = c(
      "100%" = "#D9534F", "95\u201399%" = "#E07A3E", "75\u201394%" = "#F0AD4E",
      "50\u201374%" = "#5B9BD5", "25\u201349%" = "#7AB8E0", "15\u201324%" = "#A0CBE8",
      "5\u201314%" = "#C8DEE8", "1\u20134%" = "#E0EEF5", "0%" = "#F5F5F5"
    ),
    guide = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Frequency (proportion of 241 genomes)",
    y = "Number of gene families"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 10, r = 5, b = 5, l = 5)
  )

print(freq_plot)
ggsave(file.path(out_dir, "gene_family_frequency_distribution.png"), freq_plot,
       width = 6, height = 4.5, dpi = 300, bg = "white")

# ============================================================
# Summary
# ============================================================
cat("\n=== Pan-genome tier summary ===\n")
print(tier_summary)

cat("\n=== Frequency distribution summary ===\n")
print(freq_summary)