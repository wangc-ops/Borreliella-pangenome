########################################################################
## 08_heatmap_pav_all_genomes.R
## Antigen PAV heatmap for the full genome set (deposition-bias view).
## Depends on: 01_build_pav_matrices.R (results/antigen_pav_matrix_all.csv)
## Key visual: genomes without deposited plasmid sequences are blank for
## all plasmid-borne families but retain the chromosome-encoded controls.
## Columns: genomes blocked by species; within each block, Deposited
## first then Not deposited; ward.D2 clustering within each subgroup.
## Top annotation: Species + deposition Status.
## Sanity check printed: plasmid-family rate ~ 0 and chromosome-family
## rate ~ 1 in non-carriers.
## Output: results/heatmap_pav_all_genomes.pdf
########################################################################
suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(grid)
})

w_all <- fread("results/antigen_pav_matrix_all.csv")
spcol <- intersect(names(w_all), c("sp", "species"))[1]
setnames(w_all, spcol, "sp")

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
w_all[, sp := to_sp_code(sp)]

fams <- c("OspC","OspA","DbpA","DbpB","OspB","OspD",
          "Erp","Mlp","CspA_Z","BBK32","VlsE",
          "BmpA_B","P66")
stopifnot(all(fams %in% names(w_all)))
X <- as.matrix(w_all[, ..fams]); rownames(X) <- w_all$genome

pal    <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
sp_lab <- c(bb = "B. burgdorferi", gar = "B. garinii",
            afz = "B. afzelii †", bav = "B. bavariensis")
sp_fac    <- factor(w_all$sp, levels = c("bb", "gar", "afz", "bav"))
sp_pretty <- factor(unname(sp_lab[as.character(sp_fac)]), levels = unname(sp_lab))
pal2      <- setNames(pal, unname(sp_lab))
status    <- factor(ifelse(w_all$carrier == 1, "Deposited", "Not deposited"),
                    levels = c("Deposited", "Not deposited"))

## Column order: species -> deposition status (Deposited first) -> within-group clustering
grp_key <- data.table(idx = seq_len(nrow(X)), sp = sp_fac, st = status)
setorder(grp_key, sp, st)
ord <- unlist(lapply(split(grp_key$idx, list(grp_key$sp, grp_key$st)), function(idx) {
  if (length(idx) > 2) {
    hc <- hclust(dist(X[idx, , drop = FALSE], method = "binary"), method = "ward.D2")
    idx[hc$order]
  } else idx
}))

row_split <- factor(ifelse(fams %in% c("BmpA_B", "P66"),
                           "Chromosome-encoded", "Plasmid-borne"),
                    levels = c("Plasmid-borne", "Chromosome-encoded"))

ht <- Heatmap(t(X)[fams, ord, drop = FALSE],
              name = "Detected",
              col  = c("0" = "#EBEBEB", "1" = "#1F618D"),
              rect_gp = gpar(col = "white", lwd = 0.3),
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = row_split, row_title = NULL,
              column_split = sp_pretty[ord], column_title = NULL,
              column_gap = unit(1.5, "mm"),
              show_column_names = FALSE,
              row_names_gp = gpar(fontface = "italic", fontsize = 9),
              top_annotation = HeatmapAnnotation(
                Species = sp_pretty[ord],
                Status  = status[ord],
                col = list(Species = pal2,
                           Status = c("Deposited" = "#7FB3D5",
                                      "Not deposited" = "#E59866")),
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

## Sanity check: plasmid-family rate ~ 0 / chromosome-family rate ~ 1 in non-carriers
pl  <- setdiff(fams, c("BmpA_B", "P66"))
chr <- c("BmpA_B", "P66")
chk <- w_all[, .(n = .N,
                 plasmid_fam_rate = round(mean(rowMeans(.SD[, ..pl])), 3),
                 chrom_fam_rate   = round(mean(rowMeans(.SD[, ..chr])), 3)),
             by = carrier, .SDcols = c(pl, chr)]
cat("\n== Check: carrier = FALSE should show plasmid_fam_rate = 0, chrom_fam_rate ~ 1 ==\n")
print(chk)

## Save final PDF
pdf("results/heatmap_pav_all_genomes.pdf", width = 14, height = 5.5)
draw(ht, padding = unit(c(2, 14, 2, 2), "mm"))
dev.off()
cat("Saved: results/heatmap_pav_all_genomes.pdf\n")
