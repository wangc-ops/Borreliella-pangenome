########################################################################
## 02_heatmap_pav_carriers.R
## Antigen PAV heatmap for plasmid-sequence carriers (13 families).
## Depends on: 01_build_pav_matrices.R (results/antigen_pav_matrix_carriers.csv)
## Rows:    13 families (plasmid-borne on top; chromosome-encoded controls
##          BmpA_B/P66 at the bottom)
## Columns: carrier genomes, blocked by species (block titles + gaps),
##          ward.D2 clustering on binary distance within each block
## Annotation: top = Species band + cp32 count bars; right = per-family
##             number of positive genomes
## Output: results/heatmap_pav_carriers.pdf
########################################################################
suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(grid)
})

pav <- fread("results/antigen_pav_matrix_carriers.csv")
spcol <- intersect(names(pav), c("sp", "species"))[1]
setnames(pav, spcol, "sp")

to_sp_code <- function(x) {
  x <- trimws(tolower(as.character(x)))
  if (all(x %in% c("bb", "gar", "afz", "bav"))) return(x)
  out <- rep(NA_character_, length(x))
  out[grepl("burgdorferi", x)] <- "bb"
  out[grepl("garinii", x)]     <- "gar"
  out[grepl("afzelii", x)]     <- "afz"
  out[grepl("bavariensis", x)] <- "bav"
  if (any(is.na(out))) stop("Unmappable species name(s): ",
                            paste(unique(x[is.na(out)]), collapse = ", "))
  out
}
pav[, sp := to_sp_code(sp)]

## cp32 counts (top annotation bars)
cp <- fread("data/cp32_like_burden.csv")
stopifnot(all(c("genome_id", "cp32_count") %in% names(cp)))
cp[, genome := sub("_genomic$", "", genome_id)]
cpv <- as.numeric(cp$cp32_count[match(pav$genome, cp$genome)])
stopifnot(!any(is.na(cpv)), length(cpv) == nrow(pav))
cat("cp32 counts matched for", length(cpv), "genomes\n")

fams <- c("OspC","OspA","DbpA","DbpB","OspB","OspD",
          "Erp","Mlp","CspA_Z","BBK32","VlsE",
          "BmpA_B","P66")
stopifnot(all(fams %in% names(pav)))
X <- as.matrix(pav[, ..fams]); rownames(X) <- pav$genome

pal    <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
sp_lab <- c(bb = "B. burgdorferi", gar = "B. garinii",
            afz = "B. afzelii †", bav = "B. bavariensis")
sp_fac    <- factor(pav$sp, levels = c("bb", "gar", "afz", "bav"))
sp_pretty <- factor(unname(sp_lab[as.character(sp_fac)]), levels = unname(sp_lab))
pal2      <- setNames(pal, unname(sp_lab))

## Column order: species blocks + within-block clustering
ord <- unlist(lapply(split(seq_len(nrow(X)), sp_fac), function(idx) {
  if (length(idx) > 2) {
    hc <- hclust(dist(X[idx, , drop = FALSE], method = "binary"), method = "ward.D2")
    idx[hc$order]
  } else idx
}))

row_split <- factor(ifelse(fams %in% c("BmpA_B", "P66"),
                           "Chromosome-encoded", "Plasmid-borne"),
                    levels = c("Plasmid-borne", "Chromosome-encoded"))

fam_pos <- colSums(X)[fams]   # right-side bars: positives per family

ht <- Heatmap(t(X)[fams, ord, drop = FALSE],
              name = "Detected",
              col  = c("0" = "#EBEBEB", "1" = "#1F618D"),
              rect_gp = gpar(col = "white", lwd = 0.4),
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = row_split, row_title = NULL,
              column_split = sp_pretty[ord], column_title = NULL,
              column_gap = unit(1.5, "mm"),
              show_column_names = FALSE,
              row_names_gp = gpar(fontface = "italic", fontsize = 9),
              top_annotation = HeatmapAnnotation(
                Species = sp_pretty[ord],
                `cp32 count` = anno_barplot(
                  cpv[ord], bar_width = 0.9,
                  gp = gpar(fill = "#636363", col = NA),
                  height = unit(15, "mm")),
                col = list(Species = pal2),
                show_annotation_name = TRUE,
                annotation_name_gp = gpar(fontsize = 8),
                show_legend = FALSE),
              right_annotation = rowAnnotation(
                `Positives (n)` = anno_barplot(
                  fam_pos, bar_width = 0.9,
                  gp = gpar(fill = "#1F618D", col = NA),
                  width = unit(18, "mm")),
                show_annotation_name = TRUE,
                annotation_name_gp = gpar(fontsize = 8)),
              heatmap_legend_param = list(title = "Detected"))
ht <- ht %v% HeatmapAnnotation(
  block = anno_block(gp = gpar(fill = "white", col = NA),
                     labels = levels(sp_pretty),
                     labels_gp = gpar(fontface = "bold.italic", fontsize = 9),
                     height = unit(6, "mm")),
  which = "column", show_annotation_name = FALSE)

## On-screen preview
draw(ht, padding = unit(c(2, 14, 2, 2), "mm"))

## Report: positives per family (carrier subset)
cat("\n== Positives per family (carrier subset) ==\n")
print(data.table(family = fams, positives = fam_pos))

## Save final PDF
pdf("results/heatmap_pav_carriers.pdf", width = 12, height = 5.8)
draw(ht, padding = unit(c(2, 14, 2, 2), "mm"))
dev.off()
cat("Saved: results/heatmap_pav_carriers.pdf\n")
