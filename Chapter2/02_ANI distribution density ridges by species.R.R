# ============================================================
# ANI distribution density ridges by species
# Input: long-format ANI table + summary statistics table
# Output: species-stratified density plot with mean lines and CLD labels
# ============================================================

library(ggplot2)
library(ggridges)
library(dplyr)

# --- 1. File paths (modify as needed) ---
ani_data_path <- "path/to/ani_intra_data.csv"      # long format: genome1, genome2, ani, species
stats_path <- "path/to/ani_summary_stats.csv"      # species, mean_ani, sd_ani, n, cld
out_pdf <- "ani_distribution_ridges.pdf"
out_png <- "ani_distribution_ridges.png"

# --- 2. Load data ---
ani_intra <- read.csv(ani_data_path, stringsAsFactors = FALSE)
ani_stats_sorted <- read.csv(stats_path, stringsAsFactors = FALSE)

# --- 3. Color scheme and species order ---
sp_colors <- c(
  "Borreliella burgdorferi" = "#1F618D",
  "Borreliella garinii"     = "#E74C3C",
  "Borreliella afzelii"     = "#1ABC9C",
  "Borreliella bavariensis" = "#5DADE2"
)

sp_levels <- c("Borreliella bavariensis", "Borreliella afzelii", 
               "Borreliella garinii", "Borreliella burgdorferi")

ani_intra$species <- factor(ani_intra$species, levels = sp_levels)

# --- 4. Plot ---
title_text <- "Welch's ANOVA: F(3.0, 1831.7) = 1730.24, p < 0.001, partial eta-sq = 0.16"

p_ani <- ggplot(ani_intra, aes(x = ani, y = species, fill = species)) +
  geom_density_ridges(alpha = 0.9, scale = 1.8, rel_min_height = 0.001,
                      linewidth = 0.3, color = "white") +
  # Mean lines
  geom_segment(data = ani_stats_sorted,
               aes(x = mean_ani, xend = mean_ani, 
                   y = as.numeric(factor(species, levels = sp_levels)), 
                   yend = as.numeric(factor(species, levels = sp_levels)) + 1.3),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.6, linetype = "dashed", alpha = 0.6) +
  # Sample size and mean±SD labels
  geom_text(data = ani_stats_sorted,
            aes(x = 95.2, 
                y = as.numeric(factor(species, levels = sp_levels)) + 0.5,
                label = sprintf("n = %s  %.2f%% ± %.2f%%", 
                                format(n, big.mark = ","), mean_ani, sd_ani)),
            inherit.aes = FALSE,
            color = "black", size = 3.2, hjust = 0, fontface = "plain") +
  # Games-Howell CLD labels
  geom_text(data = ani_stats_sorted,
            aes(x = 99.8, 
                y = as.numeric(factor(species, levels = sp_levels)) + 0.55,
                label = cld),
            inherit.aes = FALSE,
            color = "black", size = 5, fontface = "bold", hjust = 1) +
  scale_fill_manual(values = sp_colors) +
  scale_y_discrete(labels = function(x) gsub("Borreliella ", "B. ", x),
                   expand = expansion(add = c(0.3, 1.2))) +
  scale_x_continuous(limits = c(95, 100.2), breaks = seq(95, 100, 1),
                     expand = c(0, 0)) +
  labs(x = "Average Nucleotide Identity (%)", y = NULL, title = title_text) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic", size = 12, color = "black"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length.y = unit(2, "mm"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length.x = unit(2, "mm"),
    plot.title = element_text(size = 10.5, face = "plain", hjust = 0.5, 
                              color = "black", margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

ggsave(out_pdf, p_ani, width = 7.5, height = 5.5, dpi = 300)
ggsave(out_png, p_ani, width = 7.5, height = 5.5, dpi = 300)

print(p_ani)