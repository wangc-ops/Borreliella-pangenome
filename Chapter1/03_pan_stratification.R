# ============================================================================
# Script: pan_stratification.R
# Description: Pan-genome stratification analysis for Borreliella genus.
#              (1) Gene family frequency distribution with classification 
#                  thresholds (Core/Soft core/Shell/Cloud);
#              (2) Stratification pie chart showing proportional composition;
#              (3) CDS number boxplot across pan-genome categories with 
#                  Welch's ANOVA and Games-Howell post-hoc tests.
# Input:       Panaroo gene_presence_absence_roary.csv
# Output:      frequency_distribution.png, stratification_pie.png, 
#              CDS_boxplot.png
# Requirements: R >= 4.0; packages: data.table, ggplot2, dplyr
# ============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

# --- User-configurable path ---
# Place Panaroo output gene_presence_absence_roary.csv in ./data/ directory,
# or modify the path below to your local location.
INPUT_FILE <- "./data/gene_presence_absence_roary.csv"

# --- Color scheme (consistent with manuscript) ---
COL_CORE     <- "#D9534F"
COL_SOFT     <- "#F0AD4E"
COL_SHELL    <- "#5B9BD5"
COL_CLOUD    <- "#A0CBE8"

# --- 1. Load PAV matrix ---
dt <- fread(INPUT_FILE, header = TRUE)
strain_cols <- names(dt)[15:ncol(dt)]
n <- length(strain_cols)

cat("Loaded", n, "genomes and", nrow(dt), "gene families.\n")

# --- 2. Calculate presence frequency and stratification ---
dt[, presence := rowSums(.SD != "") / n, .SDcols = strain_cols]
dt[, n_genomes := rowSums(.SD != ""), .SDcols = strain_cols]

dt[, pan_cat := fcase(
  presence >= 0.99, "Core",
  presence >= 0.95 & presence < 0.99, "Soft core",
  presence >= 0.15 & presence < 0.95, "Shell",
  presence < 0.15, "Cloud"
)]

# Stratification thresholds
th_shell <- floor(0.15 * n)
th_soft  <- floor(0.95 * n)
th_core  <- floor(0.99 * n)

# --- 3. Frequency distribution histogram ---
freq_dist <- dt[, .N, by = n_genomes][order(n_genomes)]
freq_dist[, bar_color := fcase(
  n_genomes < th_shell, COL_CLOUD,
  n_genomes < th_soft,  COL_SHELL,
  n_genomes < th_core,  COL_SOFT,
  default =             COL_CORE
)]

p_hist <- ggplot(freq_dist, aes(x = n_genomes, y = N, fill = bar_color)) +
  geom_col(width = 1, show.legend = FALSE) +
  scale_fill_identity() +
  geom_vline(xintercept = c(th_shell, th_soft, th_core), 
             linetype = "dashed", color = "gray50", linewidth = 0.8) +
  labs(x = "Number of Genomes", y = "Number of gene families") +
  scale_x_continuous(expand = c(0, 0), limits = c(-5, n + 5),
                     breaks = c(0, 50, 100, 150, 200, 250, n)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 650),
                     breaks = seq(0, 600, 200)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.margin = margin(10, 10, 10, 10)
  )

# --- 4. Stratification pie chart ---
pie_data <- dt[, .N, by = pan_cat]
pie_data[, pan_cat := factor(pan_cat, levels = c("Core", "Soft core", "Shell", "Cloud"))]
pie_data[, pct := round(N / sum(N) * 100, 1)]
pie_data[, label := paste0(pan_cat, "\n", N, " (", pct, "%)")]

p_pie <- ggplot(pie_data, aes(x = "", y = N, fill = pan_cat)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  geom_text(aes(x = 1.25, label = label), 
            position = position_stack(vjust = 0.5), 
            size = 3.2, fontface = "bold", color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = c("Core" = COL_CORE, "Soft core" = COL_SOFT, 
                               "Shell" = COL_SHELL, "Cloud" = COL_CLOUD)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(40, 40, 40, 40)
  )

# --- 5. CDS number boxplot per category ---
counts <- dt[, lapply(.SD, function(x) sum(x != "")), by = pan_cat, .SDcols = strain_cols]
df_c <- melt(counts, id.vars = "pan_cat", variable.name = "strain", value.name = "count")
df_c[, pan_cat := factor(pan_cat, levels = c("Core", "Soft core", "Shell", "Cloud"))]

# Welch's ANOVA
welch_res <- oneway.test(count ~ pan_cat, data = df_c, var.equal = FALSE)
p_title <- ifelse(welch_res$p.value < 0.001, "p < 0.001", 
                  paste0("p = ", format(welch_res$p.value, digits = 2)))

cat("\nWelch's ANOVA:", p_title, "\n")

# Games-Howell post-hoc test
games_howell <- function(data, group_col, value_col) {
  groups <- levels(data[[group_col]])
  res <- data.frame(group1 = character(), group2 = character(), 
                    diff = numeric(), p.value = numeric(), stringsAsFactors = FALSE)
  for(i in 1:(length(groups)-1)) {
    for(j in (i+1):length(groups)) {
      g1 <- groups[i]; g2 <- groups[j]
      x1 <- data[data[[group_col]] == g1][[value_col]]
      x2 <- data[data[[group_col]] == g2][[value_col]]
      m1 <- mean(x1); m2 <- mean(x2)
      v1 <- var(x1); v2 <- var(x2)
      n1 <- length(x1); n2 <- length(x2)
      se <- sqrt(v1/n1 + v2/n2)
      t_stat <- (m1 - m2) / se
      df <- (v1/n1 + v2/n2)^2 / (v1^2/(n1^2*(n1-1)) + v2^2/(n2^2*(n2-1)))
      pval <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)
      res <- rbind(res, data.frame(group1 = g1, group2 = g2, diff = m1-m2, p.value = pval))
    }
  }
  res
}

gh_res <- games_howell(df_c, "pan_cat", "count")
cat("\nGames-Howell pairwise comparisons:\n")
print(gh_res)

# Letter assignment (a-d) based on descending mean
group_means <- df_c[, .(mean_count = mean(count)), by = pan_cat]
setorder(group_means, -mean_count)
group_means[, letter := letters[1:.N]]

# Letter position (above max value + 50 to avoid overlap)
letter_pos <- df_c[, .(ypos = max(count) + 50), by = pan_cat]
letter_pos <- merge(letter_pos, group_means[, .(pan_cat, letter)], by = "pan_cat")

# Plot
p_c <- ggplot(df_c, aes(x = pan_cat, y = count, fill = pan_cat)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "black", stroke = 1) +
  geom_text(data = letter_pos, aes(x = pan_cat, y = ypos, label = letter, fill = NULL),
            size = 5, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("Core" = COL_CORE, "Soft core" = COL_SOFT, 
                               "Shell" = COL_SHELL, "Cloud" = COL_CLOUD)) +
  labs(x = "Pan-genome category", y = "CDS number", 
       title = paste0("Welch's ANOVA ", p_title)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 11, hjust = 0.5, face = "plain"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_y_continuous(limits = c(-50, max(df_c$count) + 100),
                     breaks = seq(0, 800, 200),
                     expand = c(0, 0))

# --- 6. Preview figures ---
print(p_hist)
print(p_pie)
print(p_c)

# To save figures, uncomment below:
# ggsave("frequency_distribution.png", p_hist, width = 6, height = 4, dpi = 300)
# ggsave("stratification_pie.png", p_pie, width = 4, height = 4, dpi = 300)
# ggsave("CDS_boxplot.png", p_c, width = 5, height = 4.5, dpi = 300)