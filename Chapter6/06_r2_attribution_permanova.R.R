# ============================================================
# 06. Variance attribution of between-species structure (Fig 6d):
#     full PAV vs plasmid-associated accessory families (R2 comparison)
# Input : data/strain_plasmid_status_master_corrected.csv
#         data/Table_S18_plasmid_antigen_function_matrix.csv
#         external: gene_presence_absence_roary.csv -> EXT_PAV
# Output: results/fig6d_permanova_r2_final.csv
#         results/fig6d_r2_nulls.csv
#         results/fig6d_r2_attribution_bar.pdf
# Note  : question -- how much of the accessory-genome species structure is
#         explained by the plasmid-borne fraction? Controls: the 66
#         chromosomal accessory families plus a size-matched random null
#         (showing that small-set R2 inflation is a size artifact).
#         permutations: this test (and the antigen-spectrum PCoA) uses 999
#         permutations per M&M; the 9,999-permutation PERMANOVA of the
#         population-structure analysis is a different test. This script = 999.
#         flag convention: original S18 flags (no re-classification) ->
#         1,427 plasmid-associated (95.6%) and 66 chromosomal accessory
#         families; all R2 anchors below are locked from the actual run
#         outputs reported in the main text; the first run validates them,
#         no backfill needed.
#         n = 129 (bb/gar/bav plasmid carriers); afz excluded from PERMANOVA
#         runtime: 3 x 999-permutation PERMANOVA + betadisper ~ 1-2 min;
#         the null model (999 draws x 99 permutations; the draw count of 999
#         follows M&M, the per-draw permutation count of 99 is unspecified
#         there) is an internal calibration and does not enter the main text
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(vegan) })

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

# External input (too large to ship in data/; set the path before running)
EXT_PAV <- "external/gene_presence_absence_roary.csv"
stopifnot(file.exists(EXT_PAV))

NPERM <- 999   # 999 permutations for this test per M&M (the 9,999 of the
               # population-structure PERMANOVA is a different test)

# Flag switch (fixed at FALSE; do not enable):
#   FALSE = original S18 flags -> 1,427 plasmid / 66 chromosomal (main text & M&M)
#   TRUE  = audit-based re-classification of hsdR~~~parA and group_39 ->
#           1,425 plasmid / 68 chromosomal; archived for reference only
APPLY_FLIP <- FALSE

# ---------- master / 129 analysis genomes ----------
to_sp_code <- function(x) {
  out <- rep(NA_character_, length(x))
  out[grepl("garinii",     x, ignore.case = TRUE)] <- "gar"
  out[grepl("afzelii",     x, ignore.case = TRUE)] <- "afz"
  out[grepl("bavariensis", x, ignore.case = TRUE)] <- "bav"
  out[grepl("burgdorferi", x, ignore.case = TRUE) & is.na(out)] <- "bb"
  out
}
master <- fread(file.path(DIR_DATA, "strain_plasmid_status_master_corrected.csv"))
master[, sp := to_sp_code(species)]
stopifnot(!any(is.na(master$sp)))
an <- master[has_plasmid_corrected == TRUE & sp %in% c("bb", "gar", "bav")]
stopifnot(nrow(an) == 129)

# ---------- binary PAV matrix ----------
pav_raw <- fread(EXT_PAV)
meta_cols <- c("Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
               "Avg sequences per locus", "Avg sequences per isolate",
               "Genome Fragment", "Order within Fragment",
               "Accessory Fragment", "Accessory Order with Fragment", "QC",
               "Min group size nuc", "Max group size nuc", "Avg group size nuc")
genome_cols <- setdiff(names(pav_raw), meta_cols)
gid <- sub("_genomic$", "", genome_cols)
stopifnot(length(genome_cols) == 241, all(gid %in% master$Genome_ID))
raw <- as.matrix(pav_raw[, ..genome_cols])
mat <- t(ifelse(is.na(raw) | raw == "" | raw == "0", 0, 1))
rownames(mat) <- gid
colnames(mat) <- pav_raw$Gene
stopifnot(all(mat %in% c(0, 1)))

# ---------- tiers (241-genome convention; direct PAV counting) ----------
freq_all <- colSums(mat)
tier241 <- ifelse(freq_all >= 0.99 * 241, "Core",
           ifelse(freq_all >= 0.95 * 241, "Softcore",
           ifelse(freq_all >= 0.15 * 241, "Shell", "Cloud")))
stopifnot(sum(tier241 == "Core") == 763, sum(tier241 == "Softcore") == 9,
          sum(tier241 == "Shell") == 453, sum(tier241 == "Cloud") == 1040)
acc_fams <- names(tier241)[tier241 %in% c("Shell", "Cloud")]
stopifnot(length(acc_fams) == 1493)

# ---------- S18 flags (re-classification governed by APPLY_FLIP, see header) ----------
s18 <- fread(file.path(DIR_DATA, "Table_S18_plasmid_antigen_function_matrix.csv"))
if (APPLY_FLIP) {
  s18[Gene_Family %in% c("hsdR~~~parA", "group_39"), Plasmid_Associated := 0L]
}
plas_true <- s18$Plasmid_Associated %in% c(TRUE, "TRUE", 1, "1")
fam_plasmid_acc <- intersect(acc_fams, s18$Gene_Family[plas_true])
fam_chrom_acc   <- setdiff(acc_fams, s18$Gene_Family[plas_true])
if (APPLY_FLIP) {
  stopifnot(length(fam_plasmid_acc) == 1425, length(fam_chrom_acc) == 68)
} else {
  stopifnot(length(fam_plasmid_acc) == 1427, length(fam_chrom_acc) == 66)
}
n_chrom <- length(fam_chrom_acc)   # 66 (original flags) or 68 (re-classified);
                                   # everything downstream is parameterized on this

keep_g <- an$Genome_ID
sp_vec <- factor(an$sp[match(keep_g, an$Genome_ID)], levels = c("bb", "gar", "bav"))
sets <- list(full_pav    = colnames(mat),
             plasmid_acc = fam_plasmid_acc,
             chrom_acc   = fam_chrom_acc)

# ---------- PERMANOVA + betadisper ----------
set.seed(1)
run_one <- function(fams, label) {
  sub <- mat[keep_g, fams, drop = FALSE]
  d <- vegdist(sub, method = "jaccard", binary = TRUE)
  ad <- adonis2(d ~ sp_vec, permutations = NPERM)
  bd <- betadisper(d, sp_vec)          # harmless warning: negative square roots set to 0
  pt <- permutest(bd, permutations = NPERM)
  data.table(set = label, n_families = length(fams),
             R2 = round(ad$R2[1], 4), F = round(ad$F[1], 2),
             P_permanova = ad$`Pr(>F)`[1], P_betadisper = pt$tab$`Pr(>F)`[1])
}
res <- rbindlist(lapply(names(sets), function(nm) run_one(sets[[nm]], nm)))
res_sens <- run_one(setdiff(colnames(mat), fam_chrom_acc), paste0("full_excl_", n_chrom))
res_final <- rbind(res[set %in% c("full_pav", "plasmid_acc")], res_sens)
print(res)         # chrom_acc kept for reference (audited; not plotted)
print(res_final)
fwrite(res_final, file.path(DIR_RES, "fig6d_permanova_r2_final.csv"))
# R2 anchors (actual run outputs; main-text convention 1,427/66, 999 permutations):
stopifnot(abs(res[set == "full_pav", R2] - 0.671) < 0.001,
          abs(res[set == "plasmid_acc", R2] - 0.566) < 0.001,
          res[set == "full_pav", P_permanova] <= 0.001,
          res[set == "plasmid_acc", P_permanova] <= 0.001)

# ---------- size-matched null for the chromosomal accessory set
#   (internal calibration, not in main text; 999 draws x 99 permutations) ----------
null_r2 <- function(pool, k, times = 999) {
  replicate(times, {
    repeat {
      fs <- sample(pool, k)
      sub <- mat[keep_g, fs, drop = FALSE]
      if (all(rowSums(sub) > 0)) break
    }
    d <- vegdist(sub, method = "jaccard", binary = TRUE)
    adonis2(d ~ sp_vec, permutations = 99)$R2[1]
  })
}
set.seed(1)
r2_null <- null_r2(sets$plasmid_acc, n_chrom)
cat(sprintf("chrom_acc (%d families) R2 = %s | size-matched plasmid-pool null median = %.3f, max = %.3f\n",
            n_chrom, res[set == "chrom_acc", R2], median(r2_null), max(r2_null)))
fwrite(data.table(null_r2 = r2_null), file.path(DIR_RES, "fig6d_r2_nulls.csv"))
# Null anchor (actual run output, 66-draw convention): median R2 ~ 0.561
# (999 draws; tolerance widened to 0.01 because each draw's adonis2 uses 99
#  permutations, so the median is robust to that noise but not bit-exact)
if (!APPLY_FLIP) stopifnot(abs(median(r2_null) - 0.561) < 0.01)

# ---------- Fig 6d: two-bar comparison ----------
set_labs <- c(full_pav    = paste0("Full PAV\n(", format(length(sets$full_pav), big.mark = ","), " families)"),
              plasmid_acc = paste0("Plasmid-associated accessory\n(",
                                   format(length(sets$plasmid_acc), big.mark = ","), ")"))
res_plot <- res[set %in% c("full_pav", "plasmid_acc")]
res_plot[, set := factor(set, levels = names(set_labs))]
res_plot[, p_lab := ifelse(P_permanova <= 0.001, "italic(P) == 0.001",
                           sprintf("italic(P) == %.3f", P_permanova))]

p <- ggplot(res_plot, aes(x = set, y = R2)) +
  geom_col(fill = c("grey35", "grey70"), color = "black", width = 0.55) +
  geom_text(aes(label = p_lab), parse = TRUE, vjust = -0.5, size = 3.5) +
  scale_x_discrete(labels = set_labs) +
  scale_y_continuous(limits = c(0, max(res_plot$R2) * 1.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = expression("Variance explained by species (" * R^2 * ")")) +
  theme_classic(base_size = 11) +
  theme(axis.line  = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
print(p)
ggsave(file.path(DIR_RES, "fig6d_r2_attribution_bar.pdf"), p, width = 5, height = 4.5)
cat("=== 06 done: R2 attribution table + null model + Fig 6d exported (anchors locked) ===\n")
