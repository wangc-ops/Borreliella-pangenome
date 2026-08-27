# ============================================================================
# Script: Heaps_curve.R
# Description: Pan-genome accumulation curve and Heaps law fitting for 
#              Borreliella (Lyme disease spirochetes) based on 287 genomes.
# Input:       Panaroo gene_presence_absence_roary.csv
# Output:      Heaps curve plot
# Requirements: R >= 4.0; packages: data.table, ggplot2, dplyr, tidyr
# Parameters:  100 random permutations; power-law model n = A * N^gamma
# Note:        Set INPUT_FILE to your local path before running.
# ============================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- User-configurable path ---
INPUT_FILE <- "./data/gene_presence_absence_roary.csv"
# Example for Windows absolute path (replace with your own):
# INPUT_FILE <- "D:/your_folder/Panaroo/gene_presence_absence_roary.csv"

# --- Load data ---
dt <- fread(INPUT_FILE, header = TRUE)
strain_cols <- names(dt)[15:ncol(dt)]
n <- length(strain_cols)

mat <- as.matrix(dt[, ..strain_cols] != "")
rownames(mat) <- dt$Gene

# --- 100 permutations ---
set.seed(42)
n_perm <- 100
pan_mat <- matrix(0, nrow = n, ncol = n_perm)
core_mat <- matrix(0, nrow = n, ncol = n_perm)

for(p in 1:n_perm) {
  idx <- sample(n)
  mat_perm <- mat[, idx, drop = FALSE]
  for(k in 1:n) {
    sub <- mat_perm[, 1:k, drop = FALSE]
    pan_mat[k, p] <- sum(rowSums(sub) > 0)
    core_mat[k, p] <- sum(rowSums(sub) == k)
  }
}

# --- Summary statistics ---
plot_data <- data.frame(
  N = 1:n,
  Pan_mean = rowMeans(pan_mat),
  Pan_sd = apply(pan_mat, 1, sd),
  Core_mean = rowMeans(core_mat),
  Core_sd = apply(core_mat, 1, sd)
)

# --- Heaps law fitting ---
heaps_fit <- nls(Pan_mean ~ A * N^gamma, data = plot_data,
                 start = list(A = 1000, gamma = 0.2),
                 control = nls.control(maxiter = 200))

A_coef <- coef(heaps_fit)["A"]
gamma_coef <- coef(heaps_fit)["gamma"]

heaps_line <- data.frame(
  N = seq(1, n, length.out = 500),
  Genes = predict(heaps_fit, newdata = data.frame(N = seq(1, n, length.out = 500)))
)

# --- Plotting ---
plot_long <- plot_data %>%
  pivot_longer(cols = c(Pan_mean, Core_mean), names_to = "Type", values_to = "Genes") %>%
  mutate(
    Type = gsub("_mean", "", Type),
    SD = ifelse(Type == "Pan", Pan_sd, Core_sd)
  )

total_pan <- max(plot_data$Pan_mean)
final_core <- plot_data$Core_mean[n]

p <- ggplot() +
  geom_ribbon(data = filter(plot_long, Type == "Pan"),
              aes(x = N, ymin = Genes - SD, ymax = Genes + SD),
              fill = "#B8D4E3", alpha = 0.55) +
  geom_ribbon(data = filter(plot_long, Type == "Core"),
              aes(x = N, ymin = Genes - SD, ymax = Genes + SD),
              fill = "#E8A5A5", alpha = 0.55) +
  geom_line(data = heaps_line, aes(x = N, y = Genes),
            color = "#5B9BD5", linewidth = 1.0, linetype = "dashed", alpha = 0.9) +
  geom_line(data = filter(plot_long, Type == "Pan"),
            aes(x = N, y = Genes, color = "Pan"), linewidth = 1.3) +
  geom_line(data = filter(plot_long, Type == "Core"),
            aes(x = N, y = Genes, color = "Core"), linewidth = 1.3) +
  annotate("text", x = n * 1.03, y = total_pan,
           label = paste0(round(total_pan)), color = "#2E5C8A", size = 5.5, fontface = "bold", hjust = 0) +
  annotate("text", x = n * 1.03, y = final_core,
           label = paste0(round(final_core)), color = "#C0504D", size = 5.5, fontface = "bold", hjust = 0) +
  annotate("text", x = n * 0.55, y = total_pan * 1.08,
           label = paste0("Heaps law: n = ", round(A_coef), " × N^", sprintf("%.2f", gamma_coef)),
           color = "#5B9BD5", size = 4.5, fontface = "plain") +
  annotate("text", x = n * 0.72, y = final_core + (total_pan - final_core) * 0.08,
           label = paste0("Moderately open (γ = ", sprintf("%.2f", gamma_coef), ")"),
           color = "black", size = 5, fontface = "bold") +
  scale_color_manual(values = c("Core" = "#C0504D", "Pan" = "#2E5C8A"),
                     labels = c("Core", "Pan")) +
  scale_x_continuous(expand = c(0.02, 0), limits = c(0, n * 1.15),
                     breaks = c(0, 100, 200, 300)) +
  scale_y_continuous(expand = c(0.02, 0), limits = c(0, total_pan * 1.15),
                     breaks = seq(0, 4000, 1000)) +
  labs(x = "Number of Genomes", y = "Number of Gene Families", color = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8),
    legend.position = c(0.82, 0.72),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(20, 80, 20, 20)
  ) +
  annotate("text", x = -n * 0.05, y = total_pan * 1.12,
           label = "(a)", size = 6, fontface = "bold")

print(p)