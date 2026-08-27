# ============================================================
# Strain-level PAV heatmap with pan-genome tier annotation
# Tiers: Core (>=99%), Soft-core (95-99%), Shell (15-95%), Cloud (<15%)
# ============================================================
library(data.table); library(readxl); library(ggplot2)

pav_file     <- "path/to/gene_presence_absence_roary.csv"
species_file <- "path/to/Supplementary_Table_S1.xls"
out_dir      <- "path/to/output"

species_order_y <- c("Borreliella bavariensis", "Borreliella afzelii",
                     "Borreliella garinii", "Borreliella burgdorferi")

# --- Read PAV matrix & binarize ---
pav <- fread(pav_file, header = TRUE)
meta_cols <- c("Gene", "Non-unique Gene name", "Annotation", "No. isolates",
               "No. sequences", "Avg group size", "Genome Fragment",
               "Order within Fragment", "Accessory Fragment",
               "Accessory Order with Fragment", "QC", "Min group size nuc",
               "Max group size nuc", "Avg group size nuc")
strain_cols <- setdiff(names(pav), meta_cols)
strain_ids  <- sub("_genomic$", "", strain_cols)
pav[, (strain_cols) := lapply(.SD, function(x) as.integer(!is.na(x) & x > 0)),
    .SDcols = strain_cols]

# --- Species annotation ---
s1 <- as.data.table(read_excel(species_file))
s1[, Genome_ID := sub("_genomic$", "", Genome_ID)]
s1[, species := fcase(
  grepl("burgdorferi", species), "Borreliella burgdorferi",
  grepl("garinii", species),     "Borreliella garinii",
  grepl("afzelii", species),     "Borreliella afzelii",
  grepl("bavariensis", species), "Borreliella bavariensis",
  default = species)]
s1 <- s1[species %in% species_order_y]

matched <- strain_ids %in% s1$Genome_ID
pav_sub <- pav[, c("Gene", strain_cols[matched]), with = FALSE]
s1_sub  <- s1[Genome_ID %in% strain_ids[matched]]

# --- Global frequency & tiers ---
global_freq <- rowMeans(pav_sub[, -"Gene"])
gene_meta <- data.table(
  Gene = pav_sub$Gene, global_freq = global_freq,
  category = fifelse(global_freq >= 0.99, "Core",
                     fifelse(global_freq >= 0.95, "Soft-core",
                             fifelse(global_freq >= 0.15, "Shell", "Cloud"))))
gene_meta <- gene_meta[order(-global_freq)]
gene_order <- gene_meta$Gene

# --- Strain ordering (species, then completeness) ---
genome_to_species <- setNames(s1_sub$species, s1_sub$Genome_ID)
strain_stats <- data.table(
  genome_id = strain_cols[matched],
  species = factor(genome_to_species[strain_ids[matched]], levels = species_order_y))
miss_counts <- colSums(pav_sub[, strain_cols[matched], with = FALSE] == 0)
strain_stats[, miss_count := miss_counts[genome_id]]
strain_stats <- strain_stats[order(species, miss_count)]
strain_order <- strain_stats$genome_id

# --- Tier & species boundaries ---
category_counts <- gene_meta[, .N, by = category]
category_counts[, category := factor(category, levels = c("Core", "Soft-core", "Shell", "Cloud"))]
category_counts <- category_counts[order(category)]
category_counts[, `:=`(cum_end = cumsum(N), cum_start = cumsum(N) - N + 1)]
category_counts[, mid := (cum_start + cum_end) / 2]

species_boundaries <- c(); species_labels_y <- c()
for (sp in species_order_y) {
  idx <- which(strain_stats$species == sp)
  species_labels_y[sp] <- mean(idx)
  if (sp != tail(species_order_y, 1))
    species_boundaries <- c(species_boundaries, max(idx) + 0.5)
}
n_genes <- length(gene_order); n_strains <- length(strain_order)

# --- Melt & plot ---
pav_long <- melt(pav_sub[, c("Gene", strain_order), with = FALSE],
                 id.vars = "Gene", variable.name = "genome_id", value.name = "presence")
pav_long[, genome_id := factor(genome_id, levels = strain_order)]
pav_long[, Gene := factor(Gene, levels = gene_order)]
pav_long[, presence := factor(presence, levels = c(1, 0), labels = c("Present", "Absent"))]

pav_colors <- c("Present" = "#D9534F", "Absent" = "#A0CBE8")
p_pav <- ggplot(pav_long, aes(x = Gene, y = genome_id, fill = presence)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = pav_colors, name = NULL) +
  scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 10) +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10, color = "gray20", margin = margin(t = 8)),
        axis.title.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, "cm"),
        legend.text = element_text(size = 9), legend.margin = margin(t = -5),
        plot.margin = margin(t = 35, r = 10, b = 5, l = 60)) +
  coord_cartesian(clip = "off") +
  annotate("segment", x = 1, xend = n_genes,
           y = species_boundaries, yend = species_boundaries,
           color = "white", linewidth = 0.5, lineend = "square") +
  annotate("text", x = -40, y = species_labels_y["Borreliella bavariensis"],
           label = "italic(B.~bavariensis)", parse = TRUE, size = 3.5, hjust = 1) +
  annotate("text", x = -40, y = species_labels_y["Borreliella afzelii"],
           label = "italic(B.~afzelii)", parse = TRUE, size = 3.5, hjust = 1) +
  annotate("text", x = -40, y = species_labels_y["Borreliella garinii"],
           label = "italic(B.~garinii)", parse = TRUE, size = 3.5, hjust = 1) +
  annotate("text", x = -40, y = species_labels_y["Borreliella burgdorferi"],
           label = "italic(B.~burgdorferi)", parse = TRUE, size = 3.5, hjust = 1) +
  annotate("segment",
           x = category_counts[category == "Core", cum_start],
           xend = category_counts[category == "Core", cum_end],
           y = n_strains + 3, yend = n_strains + 3, color = "black", linewidth = 0.4) +
  annotate("text", x = category_counts[category == "Core", mid], y = n_strains + 6,
           label = "Core", size = 3, fontface = "bold") +
  annotate("segment",
           x = category_counts[category == "Soft-core", cum_start],
           xend = category_counts[category == "Soft-core", cum_end],
           y = n_strains + 3, yend = n_strains + 3, color = "black", linewidth = 0.4) +
  annotate("text", x = category_counts[category == "Soft-core", mid], y = n_strains + 6,
           label = "Soft-core", size = 3, fontface = "bold") +
  annotate("segment",
           x = category_counts[category == "Shell", cum_start],
           xend = category_counts[category == "Shell", cum_end],
           y = n_strains + 3, yend = n_strains + 3, color = "black", linewidth = 0.4) +
  annotate("text", x = category_counts[category == "Shell", mid], y = n_strains + 6,
           label = "Shell", size = 3, fontface = "bold") +
  annotate("segment",
           x = category_counts[category == "Cloud", cum_start],
           xend = category_counts[category == "Cloud", cum_end],
           y = n_strains + 3, yend = n_strains + 3, color = "black", linewidth = 0.4) +
  annotate("text", x = category_counts[category == "Cloud", mid], y = n_strains + 6,
           label = "Cloud", size = 3, fontface = "bold") +
  annotate("segment",
           x = c(category_counts[category == "Soft-core", cum_start],
                 category_counts[category == "Shell", cum_start],
                 category_counts[category == "Cloud", cum_start]),
           xend = c(category_counts[category == "Soft-core", cum_start],
                    category_counts[category == "Shell", cum_start],
                    category_counts[category == "Cloud", cum_start]),
           y = n_strains + 0.5, yend = n_strains + 3,
           color = "black", linewidth = 0.3, linetype = "dashed") +
  labs(x = "Gene families (ranked by global frequency)", y = NULL)

ggsave(file.path(out_dir, "pav_heatmap_strain_level.pdf"), p_pav,
       width = 14, height = 8, dpi = 300, bg = "white")