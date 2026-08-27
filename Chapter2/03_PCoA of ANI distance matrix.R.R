# ============================================================
# PCoA of ANI distance matrix
# Input: FastANI output + genome list + species metadata
# Output: circular scatter plot with 95% confidence ellipses
# ============================================================

library(ape)
library(ggplot2)
library(dplyr)
library(readr)

# --- 1. File paths (modify as needed) ---
fastani_path <- "path/to/fastani_output.out"
genome_list_path <- "path/to/genome_list.txt"
metadata_path <- "path/to/species_metadata.csv"   # Genome_ID, species
out_pdf <- "pcoa_ani_distance.pdf"
out_png <- "pcoa_ani_distance.png"

# --- 2. Color scheme ---
sp_colors <- c(
  "Borreliella burgdorferi" = "#1F618D",
  "Borreliella garinii"     = "#E74C3C",
  "Borreliella afzelii"     = "#1ABC9C",
  "Borreliella bavariensis" = "#5DADE2"
)
target_species <- names(sp_colors)

# --- 3. Load genome list ---
genome_list <- read_lines(genome_list_path)
genome_list <- trimws(genome_list)
genome_list <- genome_list[genome_list != ""]
genome_ids <- sub("_genomic.*$", "", basename(genome_list))

# --- 4. Load FastANI and build symmetric matrix ---
ani_raw <- read_delim(fastani_path, delim = "\t", col_names = FALSE, 
                      show_col_types = FALSE, progress = FALSE)
colnames(ani_raw) <- c("query", "ref", "ani", "frac")

ani_raw$query <- sub("_genomic.*$", "", basename(ani_raw$query))
ani_raw$ref <- sub("_genomic.*$", "", basename(ani_raw$ref))

n <- length(genome_ids)
ani_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(ani_matrix) <- genome_ids
colnames(ani_matrix) <- genome_ids
diag(ani_matrix) <- 100

for(i in 1:nrow(ani_raw)) {
  q <- ani_raw$query[i]
  r <- ani_raw$ref[i]
  if(q %in% genome_ids && r %in% genome_ids) {
    ani_matrix[q, r] <- ani_raw$ani[i]
    if(!is.na(ani_matrix[r, q])) {
      ani_matrix[q, r] <- ani_matrix[r, q] <- mean(c(ani_matrix[q, r], ani_matrix[r, q]), na.rm = TRUE)
    } else {
      ani_matrix[r, q] <- ani_raw$ani[i]
    }
  }
}

# Fill missing values with row means if any
if(sum(is.na(ani_matrix)) > 0) {
  for(i in 1:n) {
    row_mean <- mean(ani_matrix[i, ], na.rm = TRUE)
    ani_matrix[i, is.na(ani_matrix[i, ])] <- row_mean
    ani_matrix[is.na(ani_matrix[, i]), i] <- row_mean
  }
  diag(ani_matrix) <- 100
}

# --- 5. PCoA ---
ani_dist <- as.dist(1 - ani_matrix / 100)
pcoa_result <- pcoa(ani_dist)

pcoa_df <- as.data.frame(pcoa_result$vectors[, 1:2])
colnames(pcoa_df) <- c("PCo1", "PCo2")
pcoa_df$Genome <- rownames(pcoa_df)

eig <- pcoa_result$values$Eigenvalues
prop_var <- eig / sum(eig) * 100
pc1_lab <- sprintf("PCo1 (%.1f%%)", prop_var[1])
pc2_lab <- sprintf("PCo2 (%.1f%%)", prop_var[2])

# --- 6. Add species annotation ---
meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
meta <- meta %>% 
  filter(species %in% target_species) %>%
  mutate(Genome_ID = sub("_genomic.*$", "", Genome_ID))

pcoa_df$species <- meta$species[match(pcoa_df$Genome, meta$Genome_ID)]
pcoa_df <- pcoa_df[!is.na(pcoa_df$species), ]
pcoa_df$species <- factor(pcoa_df$species, levels = target_species)

# --- 7. Plot ---
p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = species, fill = species)) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.12, linewidth = 0.3) +
  geom_point(size = 2.5, alpha = 0.8, shape = 21, color = "black", stroke = 0.4) +
  scale_color_manual(values = sp_colors, name = "Species") +
  scale_fill_manual(values = sp_colors, name = "Species") +
  labs(x = pc1_lab, y = pc2_lab, title = "PCoA of ANI distance matrix") +
  scale_x_continuous(breaks = pretty(pcoa_df$PCo1, n = 5)) +
  scale_y_continuous(breaks = pretty(pcoa_df$PCo2, n = 5)) +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2.5, "mm"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold", 
                              color = "black", margin = margin(b = 10)),
    legend.position = c(0.88, 0.12),
    legend.justification = c("right", "bottom"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8, face = "italic"),
    legend.key.size = unit(3, "mm"),
    legend.background = element_rect(fill = "white", colour = "grey80", size = 0.2),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

ggsave(out_pdf, p_pcoa, width = 6.5, height = 5.5, dpi = 300)
ggsave(out_png, p_pcoa, width = 6.5, height = 5.5, dpi = 300)

print(p_pcoa)