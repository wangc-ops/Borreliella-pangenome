# ============================================================
# 08. Table S19: provenance audit of the 66 chromosomal accessory families
# Input : data/strain_plasmid_status_master_corrected.csv
#         data/Table_S18_plasmid_antigen_function_matrix.csv
#         external: gene_presence_absence_roary.csv -> EXT_PAV
# Output: results/tableS19_chromosomal_accessory_audit.csv
# Note  : flag convention = original S18 flags (no re-classification) ->
#         66 families; the six-class count anchors (16/14/11/9/8/8) below are
#         locked from the actual run outputs, no backfill needed.
#         The re-classification path (68 families) remains available via the
#         APPLY_FLIP switch for reference only; do not enable.
# ============================================================
suppressPackageStartupMessages(library(data.table))

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

APPLY_FLIP <- FALSE   # original S18 flags (consistent with script 06);
                      # TRUE = re-classified convention, archived for reference only

# External input (too large to ship in data/; set the path before running)
EXT_PAV <- "external/gene_presence_absence_roary.csv"
stopifnot(file.exists(EXT_PAV))

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

# ---------- binary PAV matrix (same rule as script 06) ----------
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

# ---------- tiers + chromosomal accessory set (governed by APPLY_FLIP) ----------
freq_all <- colSums(mat)
tier241 <- ifelse(freq_all >= 0.99 * 241, "Core",
           ifelse(freq_all >= 0.95 * 241, "Softcore",
           ifelse(freq_all >= 0.15 * 241, "Shell", "Cloud")))
acc_fams <- names(tier241)[tier241 %in% c("Shell", "Cloud")]
s18 <- fread(file.path(DIR_DATA, "Table_S18_plasmid_antigen_function_matrix.csv"))
if (APPLY_FLIP) {
  s18[Gene_Family %in% c("hsdR~~~parA", "group_39"), Plasmid_Associated := 0L]
}
plas_true <- s18$Plasmid_Associated %in% c(TRUE, "TRUE", 1, "1")
fam_chrom_acc <- setdiff(acc_fams, s18$Gene_Family[plas_true])
stopifnot(length(fam_chrom_acc) == if (APPLY_FLIP) 68 else 66)

# ---------- audit table ----------
audit <- data.table(family = fam_chrom_acc)
audit <- merge(audit, pav_raw[, .(family = Gene, Annotation)], by = "family")
audit[, tier := tier241[family]]
audit[, n_total := colSums(mat[, audit$family, drop = FALSE])]
for (s in c("bb", "gar", "bav", "afz")) {
  gs <- master[sp == s, Genome_ID]
  audit[, (paste0("freq_", s)) := round(colMeans(mat[gs, audit$family, drop = FALSE]), 2)]
}
batch_ids <- master[has_plasmid_corrected == FALSE &
                    grepl("^GCF_02657|^GCF_02664", Genome_ID), Genome_ID]
stopifnot(length(batch_ids) == 98)
audit[, n_batch := colSums(mat[batch_ids, audit$family, drop = FALSE])]

audit[, flag := fcase(
  n_total <= 2 & n_batch >= 1, "batch singleton",
  n_total <= 2,                "rare singleton",
  freq_bb == 0 & freq_gar == 0 & freq_bav == 0, "afz-specific variant",
  !grepl("^group_", family),   "conserved housekeeping variant",
  n_total >= 30 & (pmax(freq_bb, freq_gar, freq_bav) -
                   pmin(freq_bb, freq_gar, freq_bav)) >= 0.5,
                               "lineage-block variant",
  default = "rare cloud (unresolved)"
)]
print(audit[, .N, by = flag][order(-N)])
if (!APPLY_FLIP) {
  # 66-family convention: anchors locked from the actual run outputs
  fc <- audit[, .N, by = flag]
  stopifnot(fc[flag == "batch singleton", N] == 16,
            fc[flag == "lineage-block variant", N] == 14,
            fc[flag == "rare cloud (unresolved)", N] == 11,
            fc[flag == "afz-specific variant", N] == 9,
            fc[flag == "rare singleton", N] == 8,
            fc[flag == "conserved housekeeping variant", N] == 8)
} else {
  # 68-family convention (reference only): classification of the two
  # re-classified families --
  #   hsdR~~~parA: no "group_" prefix and n_total = 60, must fall into
  #                conserved housekeeping (8 -> 9)
  #   group_39   : n_total = 59, falls into lineage-block or unresolved
  #                depending on its species-frequency profile
  # TODO(backfill after rerun): lock the six-class counts as assertions
}
fwrite(audit[order(flag, -n_total)],
       file.path(DIR_RES, "tableS19_chromosomal_accessory_audit.csv"))

# ---------- complementary-locus check: thrS split (expected exactly 1 copy
#   per genome across all 241 genomes) ----------
cat("thrS locus copies per genome (group_468 + thrS):\n")
print(table(rowSums(mat[, c("group_468", "thrS")])))
stopifnot(all(rowSums(mat[, c("group_468", "thrS")]) == 1))
cat(sprintf("=== 08 done: Table S19 audit table exported (%d-family convention, APPLY_FLIP=%s) ===\n",
            length(fam_chrom_acc), APPLY_FLIP))
