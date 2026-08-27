########################################################################
## 04_pcoa_antigen_profiles.R
## PCoA of antigen profiles (carrier subset; variable plasmid-borne
## families; binary Jaccard) + PERMANOVA / betadisper on bb/gar/bav.
## Depends on: 01_build_pav_matrices.R (results/antigen_pav_matrix_carriers.csv)
## Final design: point shape + colour per species; convex hulls outline
## each group (covariance-degenerate groups cannot support ellipses;
## hulls carry no statistical assumptions); test statistics are printed
## in the plot caption.
## afz is display-only and daggered in the legend.
## Outputs: results/pcoa_coordinates.csv, results/pcoa_antigen_profiles.pdf
########################################################################
suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ggplot2)
})

dir.create("results", showWarnings = FALSE)

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

fam_all <- c("OspA","OspB","OspC","OspD","DbpA","DbpB","Erp","Mlp",
             "CspA_Z","BBK32","VlsE","BmpA_B","P66")
plasmid_fams <- setdiff(fam_all, c("BmpA_B", "P66"))   # drop chromosome controls
X <- as.matrix(pav[, ..plasmid_fams]); rownames(X) <- pav$genome

keep <- colSums(X) > 0 & colSums(X) < nrow(X)
cat("Invariant families removed:", paste(plasmid_fams[!keep], collapse = ", "), "\n")
cat("Variable families analysed:", paste(colnames(X)[keep], collapse = ", "), "\n")
Xv <- X[, keep, drop = FALSE]

dm  <- vegdist(Xv, method = "jaccard", binary = TRUE)
pc  <- cmdscale(dm, k = 2, eig = TRUE)
pct <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
cat(sprintf("PCo1 = %.1f%%, PCo2 = %.1f%%\n", pct[1], pct[2]))

coords <- data.table(genome = rownames(pc$points), sp = pav$sp,
                     PCo1 = pc$points[, 1], PCo2 = pc$points[, 2])
fwrite(coords, "results/pcoa_coordinates.csv")

## ---- Statistics (bb/gar/bav; afz excluded) ----
idx    <- pav$sp %in% c("bb", "gar", "bav")
dm_mat <- as.matrix(dm)
grp    <- droplevels(factor(pav$sp[idx]))
ad <- adonis2(as.dist(dm_mat[idx, idx]) ~ grp)
cat("\n== PERMANOVA (bb/gar/bav) ==\n"); print(ad)
bd <- betadisper(as.dist(dm_mat[idx, idx]), grp)
bd_an <- anova(bd)
cat("== betadisper ==\n"); print(bd_an)
cat("Mean within-group distances:\n"); print(tapply(bd$distances, bd$group, mean))

p_lab <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
anno <- sprintf("PERMANOVA (bb/gar/bav, n = %d): R² = %.3f, F = %.1f, %s\nbetadisper: F = %.2f, %s",
                sum(idx), ad$R2[1], ad$F[1], p_lab(ad$`Pr(>F)`[1]),
                bd_an$`F value`[1], p_lab(bd_an$`Pr(>F)`[1]))
cat("\nPlot caption text:\n", anno, "\n")

pal    <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
shp    <- c(bb = 16, gar = 17, afz = 15, bav = 18)   # circle/triangle/square/diamond
sp_lab <- c(bb = "B. burgdorferi", gar = "B. garinii",
            afz = "B. afzelii †", bav = "B. bavariensis")
coords[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"))]

## ---- Convex hulls per group (duplicated points are fine; chull picks vertices) ----
hulls <- coords[, .SD[grDevices::chull(PCo1, PCo2), ], by = sp]

p <- ggplot(coords, aes(PCo1, PCo2)) +
  geom_polygon(data = hulls,
               aes(fill = sp, color = sp, group = sp),
               alpha = 0.12, linewidth = 0.6, linetype = "dashed",
               show.legend = FALSE) +
  geom_point(aes(color = sp, shape = sp), size = 2.4, alpha = 0.55) +
  scale_color_manual(values = pal, labels = sp_lab, name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_shape_manual(values = shp, labels = sp_lab, name = NULL) +
  labs(x = sprintf("PCo1 (%.1f%%)", pct[1]),
       y = sprintf("PCo2 (%.1f%%)", pct[2]),
       caption = anno) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.4))) +
  theme_classic(base_size = 10) +
  theme(legend.text  = element_text(face = "italic"),
        plot.caption = element_text(hjust = 0, size = 8, color = "grey20"),
        axis.line    = element_line(color = "black"),
        panel.grid   = element_blank())
print(p)
ggsave("results/pcoa_antigen_profiles.pdf", p,
       width = 6.5, height = 5.2, device = cairo_pdf)
cat("Saved: results/pcoa_antigen_profiles.pdf, results/pcoa_coordinates.csv\n")
