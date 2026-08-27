# ============================================================================
# Script: ANI_distribution.R
# Description: Average Nucleotide Identity (ANI) distribution histogram.
#              Compares intra-species vs inter-species ANI values with 
#              95% species boundary threshold.
# Input:       1. FastANI all-vs-all output (all pairs, one per line)
#              2. Species annotation table (strain ID -> species name)
# Output:      ANI_distribution_histogram.png
# Requirements: R >= 4.0; packages: data.table, ggplot2
# External:    FastANI v1.33 (for generating input ANI matrix if needed)
# ============================================================================

library(data.table)
library(ggplot2)

# --- User-configurable paths ---
# Place input files in ./data/ directory before running.
# FastANI output format: Query Ref ANI Mapped Total (no header)
ANI_FILE     <- "./data/fastani_all_vs_all.out"
SPECIES_FILE <- "./data/species_annotation.csv"

# --- Color scheme ---
COL_INTRA  <- "#4472C4"   # Deep blue: intra-species
COL_INTER  <- "#D9534F"   # Red: inter-species
COL_LINE   <- "#D9534F"   # Threshold line

# --- 1. Load FastANI output ---
# Expected format (no header): Query_path Ref_path ANI% Mapped_fragments Total_fragments
ani <- fread(ANI_FILE, header = FALSE, 
             col.names = c("Query", "Ref", "ANI", "Mapped", "Total"))

# Extract basename and remove .fna extension for strain ID matching
ani[, Query_ID := basename(tools::file_path_sans_ext(Query))]
ani[, Ref_ID   := basename(tools::file_path_sans_ext(Ref))]

# Keep only unique pairs (upper triangle, no self-comparisons)
ani <- ani[Query_ID < Ref_ID]

cat("Total unique comparisons loaded:", nrow(ani), "\n")

# --- 2. Load species annotation ---
# Expected CSV format (with header): Genome_ID,Species
sp <- fread(SPECIES_FILE, header = TRUE)
setnames(sp, c("Genome_ID", "Species"))

# --- 3. Merge species info and classify ---
ani <- merge(ani, sp[, .(Genome_ID, Sp1 = Species)], 
             by.x = "Query_ID", by.y = "Genome_ID", all.x = TRUE)
ani <- merge(ani, sp[, .(Genome_ID, Sp2 = Species)], 
             by.x = "Ref_ID", by.y = "Genome_ID", all.x = TRUE)

# Remove unannotated pairs (if any)
ani <- ani[!is.na(Sp1) & !is.na(Sp2)]

ani[, Comparison := fifelse(Sp1 == Sp2, "Intra-species", "Inter-species")]

# --- 4. Summary statistics ---
n_intra <- ani[Comparison == "Intra-species", .N]
n_inter <- ani[Comparison == "Inter-species", .N]
mean_intra <- ani[Comparison == "Intra-species", mean(ANI)]
mean_inter <- ani[Comparison == "Inter-species", mean(ANI)]

cat("Intra-species:", n_intra, "| Mean:", round(mean_intra, 2), "%\n")
cat("Inter-species:", n_inter, "| Mean:", round(mean_inter, 2), "%\n")

# --- 5. Plot ---
p <- ggplot(ani, aes(x = ANI, fill = Comparison)) +
  geom_histogram(binwidth = 0.3, color = "white", linewidth = 0.2, 
                 position = "identity", alpha = 0.85) +
  geom_vline(xintercept = 95, linetype = "dashed", color = COL_LINE, linewidth = 1) +
  
  annotate("text", x = 95.3, y = Inf, label = "Species boundary\n(95% ANI)", 
           color = COL_LINE, size = 3.5, fontface = "bold", hjust = 0, vjust = 1.5) +
  
  annotate("text", x = mean_inter + 1, y = 5500, 
           label = paste0("Inter-species\nmean: ", round(mean_inter, 2), "%"), 
           color = COL_INTER, size = 3.5, fontface = "bold", hjust = 0) +
  annotate("text", x = mean_intra - 1, y = 4500, 
           label = paste0("Intra-species\nmean: ", round(mean_intra, 2), "%"), 
           color = "#2E5C8A", size = 3.5, fontface = "bold", hjust = 1) +
  
  scale_fill_manual(values = c("Intra-species" = COL_INTRA, "Inter-species" = COL_INTER),
                    labels = c(paste0("Intra-species (n = ", n_intra, ")"),
                               paste0("Inter-species (n = ", n_inter, ")"))) +
  scale_x_continuous(limits = c(88, 100.5), breaks = seq(88, 100, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 6500), breaks = seq(0, 6000, 2000), expand = c(0, 0)) +
  labs(x = "Average Nucleotide Identity (%)", y = "Number of Comparisons", fill = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    legend.position = c(0.35, 0.92),
    legend.background = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# To save figure, uncomment:
# ggsave("./output/ANI_distribution_histogram.png", p, width = 6, height = 4.5, dpi = 300)