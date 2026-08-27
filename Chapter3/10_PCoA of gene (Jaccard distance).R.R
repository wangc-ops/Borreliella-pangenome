# ============================================================
# PCoA of gene presence/absence (Jaccard distance)
# Convex hulls instead of ellipses for boundary indication
# Species: Borreliella burgdorferi / garinii / afzelii / bavariensis
# ============================================================

library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(vegan)

# --- 1. Paths (replace with your local paths) ---
pav_file      <- "path/to/panaroo_four_species_complete/gene_presence_absence_roary.csv"
species_file  <- "path/to/Supplementary_Table_S1.xls"
out_dir       <- "path/to/output"

# --- 2. Parameters ---
sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_short_levels <- c("B. burgdorferi", "B. garinii", "B. afzelii", "B. bavariensis")

sp_colors <- c(
  "Borreliella burgdorferi"  = "#1F618D",
  "Borreliella garinii"      = "#E74C3C",
  "Borreliella afzelii"      = "#1ABC9C",
  "Borreliella bavariensis"  = "#5DADE2"
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

# --- 5. Match genomes ---
genome_map <- data.table(raw_col = strain_cols, Genome_ID = sub("_genomic$", "", strain_cols))
matched_dt <- genome_map[s1, on = "Genome_ID", nomatch = 0L]
pav_sub <- pav[, c("Gene", matched_dt$raw_col), with = FALSE]

# --- 6. PCoA ---
# Transpose: rows = genomes, cols = genes
pav_t <- as.data.frame(t(pav_sub[, -"Gene"]))
colnames(pav_t) <- pav_sub$Gene
pav_t$Genome_ID <- rownames(pav_t)

pav_t <- merge(pav_t, matched_dt[, .(raw_col, species)], 
               by.x = "Genome_ID", by.y = "raw_col", all.x = TRUE)

mat <- as.matrix(pav_t[, -c("Genome_ID", "species")])
jac_dist <- vegdist(mat, method = "jaccard", binary = TRUE)

pcoa <- cmdscale(jac_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  Genome_ID = pav_t$Genome_ID,
  species = pav_t$species,
  PCo1 = pcoa$points[, 1],
  PCo2 = pcoa$points[, 2]
)

# Variance explained
eig <- pcoa$eig
var_explained <- round(eig[1:2] / sum(eig) * 100, 1)

# --- 7. Convex hulls (stable boundary indication) ---
hulls <- pcoa_df %>%
  group_by(species) %>%
  slice(chull(PCo1, PCo2)) %>%
  ungroup()

# --- 8. Plot ---
fig_s4d <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2)) +
  geom_polygon(data = hulls, aes(fill = species), alpha = 0.08, color = NA, show.legend = FALSE) +
  geom_point(aes(fill = species), shape = 21, size = 2.5, color = "black", alpha = 0.85) +
  scale_fill_manual(values = sp_colors, labels = sp_short_levels, name = "Species") +
  labs(
    x = sprintf("PCo1 (%.1f%%)", var_explained[1]),
    y = sprintf("PCo2 (%.1f%%)", var_explained[2]),
    title = "PCoA of gene presence/absence (Jaccard)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.85, 0.15),
    legend.background = element_blank(),
    legend.text = element_text(face = "italic"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(out_dir, "pcoa_jaccard_pangenome.png"), fig_s4d, 
       width = 8, height = 7, dpi = 300)