# ============================================================
# 07. Table S18 v2: regenerate the Antigen / Antigen_Family columns with the
#     13-family mapping
# Input : data/Table_S18_plasmid_antigen_function_matrix.csv  (previous version;
#                                                             flag columns verified)
#         data/antigen13_genes.csv                            (from script 01)
#         data/coupling_13family_135.csv                      (from script 04; cross-check)
#         external: gene_presence_absence_roary.csv -> EXT_PAV
# Output: results/Table_S18_plasmid_antigen_function_matrix_v2.csv
# Note  : the antigen-presence counts of the previous S18 use an outdated
#         convention (ospA 90 / ospC 38 / dbpA_B 90 / ospB 88 / ospD 45,
#         inconsistent with the final 134-135 carriers); the Function_Category
#         and Plasmid_Associated columns were audited and are kept unchanged
# ============================================================
suppressPackageStartupMessages(library(data.table))

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

# External input (too large to ship in data/; set the path before running)
EXT_PAV <- "external/gene_presence_absence_roary.csv"
stopifnot(file.exists(EXT_PAV))

a13     <- fread(file.path(DIR_DATA, "antigen13_genes.csv"))
s18_old <- fread(file.path(DIR_DATA, "Table_S18_plasmid_antigen_function_matrix.csv"))
cpl     <- fread(file.path(DIR_DATA, "coupling_13family_135.csv"))
stopifnot(all(c("genome", "locus", "fam13") %in% names(a13)))
stopifnot(uniqueN(a13$fam13) == 13, nrow(s18_old) == 2265)

# ---------- 1. PAV to long format: family x genome x locus (";" split) ----------
raw <- fread(EXT_PAV)
meta_cols <- c("Gene","Non-unique Gene name","Annotation","No. isolates",
               "No. sequences","Avg sequences per isolate","Genome Fragment",
               "Order within Fragment","Accessory Fragment","Accessory Order with Fragment",
               "QC","Min group size nuc","Max group size nuc","Avg group size nuc")
genome_cols <- setdiff(names(raw), meta_cols)
stopifnot(length(genome_cols) == 241, all(endsWith(genome_cols, "_genomic")))
stripped <- sub("_genomic$", "", genome_cols)
stopifnot(uniqueN(stripped) == 241, all(a13$genome %in% stripped))

long <- melt(raw, id.vars = "Gene", measure.vars = genome_cols,
             variable.name = "genome", value.name = "cell")
long[, genome := sub("_genomic$", "", as.character(genome))]
long <- long[!is.na(cell) & cell != ""]
long <- long[, .(locus = trimws(unlist(strsplit(cell, "[;\t]")))),
             by = .(Gene_Family = Gene, genome)]

# ---------- 2. antigen locus -> Panaroo family ----------
m <- merge(a13, long, by = c("genome", "locus"))
match_rate <- nrow(m) / nrow(a13)
cat(sprintf("antigen gene match rate: %d/%d = %.1f%%\n", nrow(m), nrow(a13), match_rate * 100))
unmatched <- a13[!paste(genome, locus) %in% paste(m$genome, m$locus)]
if (nrow(unmatched) > 0) {
  cat("unmatched antigen genes (first 10 rows; inspect manually):\n"); print(head(unmatched, 10))
}
stopifnot(match_rate > 0.9)

# ---------- 3. family -> fam13 majority vote + conflict report ----------
votes <- m[, .(n = .N), by = .(Gene_Family, fam13)]
conflict_tbl <- votes[, .(n_fam13_types = .N), by = Gene_Family]
best <- votes[order(Gene_Family, -n), .SD[1], by = Gene_Family]
best <- merge(best, conflict_tbl, by = "Gene_Family")
best[, conflict := n_fam13_types > 1]
cat("Panaroo families mapped to antigens:", nrow(best),
    " | with multi-fam13 conflicts:", sum(best$conflict), "\n")
print(best[, .(n_panaroo_families = .N), by = fam13][order(-n_panaroo_families)])
if (any(best$conflict)) {
  cat("conflict family details (majority vote takes the top count):\n")
  print(votes[Gene_Family %in% best[conflict == TRUE, Gene_Family]][order(Gene_Family, -n)])
}

# ---------- 4. cross-check against the coupling table: per-genome antigen
#   gene count == ag_total ----------
genes_pg <- a13[, .(n_genes = .N), by = genome]
xc <- merge(cpl[, .(genome, ag_total)], genes_pg, by = "genome")
mism <- xc[n_genes != ag_total]
if (nrow(mism) > 0) { cat("mismatched genomes:\n"); print(mism) }
stopifnot(nrow(xc) == 135, nrow(mism) == 0)
cat("check passed: per-genome antigen gene counts identical to ag_total\n")

# ---------- 5. localization spot-check (printed for review, non-blocking:
#   P66/BmpA_B expected chromosomal, the rest plasmid) ----------
loc_check <- merge(best, s18_old[, .(Gene_Family, Plasmid_Associated)], by = "Gene_Family")
print(loc_check[, .(n_families = .N), by = .(fam13, Plasmid_Associated)][order(fam13, Plasmid_Associated)])

# ---------- 6. update S18 and export ----------
s18v2 <- copy(s18_old)
s18v2[, c("Antigen", "Antigen_Family") := NULL]
s18v2 <- merge(s18v2, best[, .(Gene_Family, Antigen_Family = fam13)],
               by = "Gene_Family", all.x = TRUE)
s18v2[is.na(Antigen_Family), Antigen_Family := "Non-antigen"]
s18v2[, Antigen := as.integer(Antigen_Family != "Non-antigen")]
stopifnot(nrow(s18v2) == 2265)

chg <- table(old = s18_old$Antigen, new = s18v2$Antigen)
cat("old -> new Antigen column changes:\n"); print(chg)
cat("antigen families under the new mapping:", s18v2[, sum(Antigen)], "\n")

fwrite(s18v2, file.path(DIR_RES, "Table_S18_plasmid_antigen_function_matrix_v2.csv"))
cat("=== 07 done: Table S18 v2 written ===\n")
