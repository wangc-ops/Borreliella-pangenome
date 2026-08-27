# ============================================================
# Host-interaction gene PAV heatmap
# ComplexHeatmap with species top-bar + functional sidebar
# Species: Borreliella burgdorferi / garinii / afzelii / bavariensis
# ============================================================

library(data.table)
library(readxl)
library(ComplexHeatmap)
library(circlize)

# --- 1. Paths (replace with your local paths) ---
vf_pav_file  <- "path/to/vf_pav_matrix.csv"
species_file <- "path/to/Supplementary_Table_S1.xls"
out_dir      <- "path/to/output"

# --- 2. Parameters ---
sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_short <- c("B. burgdorferi", "B. garinii", "B. afzelii", "B. bavariensis")

# --- 3. Read S1 metadata ---
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

# --- 4. Read VF PAV matrix ---
vf_pav <- fread(vf_pav_file)

meta_cols <- c("Gene", "Non-unique Gene name", "Annotation", "vfg_id", "VF_Name", 
               "VFcategory", "Function", "pident", "evalue")
genome_cols <- setdiff(names(vf_pav), meta_cols)

pav_mat <- as.matrix(vf_pav[, ..genome_cols])
rownames(pav_mat) <- vf_pav$Gene

# --- 5. Genome ordering by species ---
genome_map <- data.table(col = genome_cols, Genome_ID = sub("_genomic$", "", genome_cols))
genome_map <- merge(genome_map, s1[, .(Genome_ID, species)], by = "Genome_ID", all.x = TRUE)
genome_map[, species := factor(species, levels = sp_levels)]

col_order <- genome_map[order(species), col]
pav_mat <- pav_mat[, col_order]
genome_species <- genome_map[order(species), species]

# --- 6. Functional classification (manual, based on Annotation) ---
gene_annot <- vf_pav[, .(Gene, Annotation)]
gene_annot[, func_class := fcase(
  grepl("flagell|chemotaxis|motility", tolower(Annotation)), 
  "Flagella & Motility",
  grepl("dbp|revA|fibronectin|adhesin|integrin|p66|decorin|bbk32", tolower(Annotation)), 
  "Adhesin & ECM-binding",
  grepl("complement|crasp|cspZ|outer surface|osp|erp|vls|variable|antigen", tolower(Annotation)), 
  "Immune Evasion & Surface",
  default = "Metabolic & Other"
)]

# Refine specific genes by name
gene_annot[Gene %in% c("dbpA","dbpB","revA","group_1042","group_602","group_977"), 
           func_class := "Adhesin & ECM-binding"]
gene_annot[Gene %in% c("cspZ","ospA","ospB","ospC","ospE","vlsE1",
                       "group_1155","group_1156","group_1121"), 
           func_class := "Immune Evasion & Surface"]
gene_annot[Gene %in% c("erpB1~~~erpC","erpC","erpC~~~erpA~~~erpP","erpC~~~erpM",
                       "erpG","erpH","erpQ","erpX"), 
           func_class := "Immune Evasion & Surface"]

# --- 7. Gene category: Core/shared vs species-specific ---
vf_long <- melt(vf_pav[, c("Gene", genome_cols), with = FALSE], 
                id.vars = "Gene", variable.name = "col", value.name = "presence")
vf_long <- merge(vf_long, genome_map[, .(col, species)], by = "col")

gene_sp <- vf_long[, .(
  in_burg = any(species == "Borreliella burgdorferi" & presence == 1),
  in_gar  = any(species == "Borreliella garinii"     & presence == 1),
  in_afz  = any(species == "Borreliella afzelii"     & presence == 1),
  in_bav  = any(species == "Borreliella bavariensis" & presence == 1)
), by = Gene]

gene_sp[, category := fcase(
  in_burg & !in_gar & !in_afz & !in_bav, "B. burgdorferi-specific",
  !in_burg & in_gar & !in_afz & !in_bav, "B. garinii-specific",
  !in_burg & !in_gar & in_afz & !in_bav, "B. afzelii-specific",
  !in_burg & !in_gar & !in_afz & in_bav, "B. bavariensis-specific",
  default = "Core (shared)"
)]

# --- 8. Merge & filter ---
gene_meta <- merge(gene_sp, gene_annot[, .(Gene, func_class)], by = "Gene")

filter_mode <- "host_focus"  # "host_focus" or "all"

if (filter_mode == "host_focus") {
  gene_prev <- vf_long[, .(min_prev = min(presence), max_prev = max(presence)), by = Gene]
  gene_meta <- merge(gene_meta, gene_prev, by = "Gene")
  keep <- gene_meta[, func_class %in% c("Adhesin & ECM-binding", "Immune Evasion & Surface") |
                      (func_class == "Flagella & Motility" & max_prev < 1) |
                      (func_class == "Metabolic & Other" & max_prev < 1)]
  gene_meta <- gene_meta[keep]
  cat("Host-focus mode:", nrow(gene_meta), "genes\n")
} else {
  cat("Full mode:", nrow(gene_meta), "genes\n")
}

cat_levels <- c("Core (shared)", "B. burgdorferi-specific", "B. garinii-specific", 
                "B. afzelii-specific", "B. bavariensis-specific")
gene_meta[, category := factor(category, levels = cat_levels)]
setorder(gene_meta, category, func_class, Gene)

pav_mat <- pav_mat[gene_meta$Gene, ]

# --- 9. Colors ---
sp_colors <- c(
  "Borreliella burgdorferi" = "#1F618D",
  "Borreliella garinii"     = "#E74C3C",
  "Borreliella afzelii"     = "#1ABC9C",
  "Borreliella bavariensis" = "#5DADE2"
)
cat_colors <- c(
  "Core (shared)"           = "grey40",
  "B. burgdorferi-specific" = "#1F618D",
  "B. garinii-specific"     = "#E74C3C",
  "B. afzelii-specific"     = "#1ABC9C",
  "B. bavariensis-specific" = "#5DADE2"
)
func_colors <- c(
  "Flagella & Motility"       = "#5B9BD5",
  "Adhesin & ECM-binding"     = "#D9534F",
  "Immune Evasion & Surface"  = "#F0AD4E",
  "Metabolic & Other"         = "grey80"
)

# --- 10. Top annotation ---
top_anno <- HeatmapAnnotation(
  Species = genome_species,
  col = list(Species = sp_colors),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    title = "Species", 
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels = sp_short,
    labels_gp = gpar(fontsize = 8, fontface = "italic")
  )
)

# --- 11. Left annotation (clean, no text labels) ---
left_anno <- rowAnnotation(
  Category = gene_meta$category,
  Function = gene_meta$func_class,
  col = list(Category = cat_colors, Function = func_colors),
  show_annotation_name = c(FALSE, FALSE),
  width = unit(8, "mm"),
  show_legend = c(TRUE, TRUE),
  annotation_legend_param = list(
    Category = list(title = "Gene category", title_gp = gpar(fontsize = 9)),
    Function = list(title = "Function", title_gp = gpar(fontsize = 9))
  )
)

# --- 12. Heatmap ---
ht <- Heatmap(
  pav_mat,
  name = "PAV",
  col = c("0" = "white", "1" = "#D9534F"),
  top_annotation = top_anno,
  left_annotation = left_anno,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(10, "cm"),
  column_split = genome_species,
  column_gap = unit(1, "mm"),
  column_title = "Virulence gene",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_title = "Genomes",
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_title_rot = 90,
  column_title_rot = 0,
  row_split = gene_meta$category,
  row_gap = unit(2, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "Status", 
    at = c(0, 1), 
    labels = c("Absent", "Present"),
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  )
)

# --- 13. Draw ---
draw(ht)

# --- 14. Save (uncomment after checking) ---
# pdf(file.path(out_dir, "host_interaction_gene_pav_heatmap.pdf"), width = 10, height = 8)
# draw(ht)
# dev.off()