########################################################################
## 01_build_pav_matrices.R
## Build the 13-family presence/absence (PAV) matrices for all genomes
## and for plasmid-sequence carriers, plus per-species detection rates.
## RUN FIRST: scripts 02/04/06/07/08 depend on the matrices written here.
##
## Expected layout: run from the project root with inputs in data/.
## Inputs  (data/):    strain_plasmid_status_master_corrected.csv
##                     antigen_pav_reconciled.csv
##                     antigen_copy_numbers.csv
##                     antigen_pav_long.tsv
## Outputs (results/): antigen_pav_matrix_carriers.csv
##                     antigen_pav_matrix_all.csv
##                     antigen_detection_rates_by_species.csv
########################################################################
suppressPackageStartupMessages(library(data.table))

dir.create("results", showWarnings = FALSE)

master_path <- "data/strain_plasmid_status_master_corrected.csv"
recon_path  <- "data/antigen_pav_reconciled.csv"
copy_path   <- "data/antigen_copy_numbers.csv"
long_path   <- "data/antigen_pav_long.tsv"

to_pres <- function(x) {
  if (is.logical(x))   return(x)
  if (is.character(x)) return(tolower(trimws(x)) %in% c("1", "true"))
  as.numeric(x) > 0
}

## Map full species names to codes (the master table stores full names);
## passes through unchanged if already coded.
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

## ---- 1) Reconciled long table: 6 osp/dbp families ----
rec <- fread(recon_path)
stopifnot(all(c("Genome_ID", "antigen", "detected") %in% names(rec)))
rec_long <- rec[, .(family  = antigen,
                    genome  = sub("_genomic$", "", Genome_ID),
                    present = to_pres(detected))]
norm_fam <- c("ospa"="OspA", "ospb"="OspB", "ospc"="OspC",
              "ospd"="OspD", "dbpa"="DbpA", "dbpb"="DbpB")
rec_long[, family := fifelse(tolower(family) %in% names(norm_fam),
                             unname(norm_fam[tolower(family)]), family)]

## ---- 2) Copy-number table: Erp / Mlp (present = copies > 0) ----
cn <- fread(copy_path)
stopifnot(all(c("Genome_ID", "erp_copies", "mlp_copies") %in% names(cn)))
cn_long <- rbind(
  cn[, .(family = "Erp", genome = sub("_genomic$", "", Genome_ID),
         present = as.numeric(erp_copies) > 0)],
  cn[, .(family = "Mlp", genome = sub("_genomic$", "", Genome_ID),
         present = as.numeric(mlp_copies) > 0)])

## ---- 3) Legacy long table: 5 retained families ----
old <- fread(long_path)
stopifnot(all(c("family", "genome_id", "presence") %in% names(old)))
keep5 <- c("CspA_Z", "BBK32", "VlsE", "BmpA_B", "P66")
old_long <- old[family %in% keep5,
                .(family, genome = sub("_genomic$", "", genome_id),
                  present = to_pres(presence))]

## ---- 4) Merge + hard assertion (all 13 families must be present) ----
all_long <- unique(rbind(rec_long, cn_long, old_long), by = c("family", "genome"))
expected <- c("OspA","OspB","OspC","OspD","DbpA","DbpB","Erp","Mlp",
              "CspA_Z","BBK32","VlsE","BmpA_B","P66")
miss <- setdiff(expected, unique(all_long$family))
if (length(miss) > 0) stop("[assertion failed] Missing families: ",
                           paste(miss, collapse = ", "))
cat("All 13 families present:", paste(sort(unique(all_long$family)), collapse = ", "), "\n")

## ---- 5) Wide table + merge with master (species mapped, asserted) ----
mast <- fread(master_path)
stopifnot(all(c("Genome_ID", "species", "has_plasmid_corrected") %in% names(mast)))
mast[, genome  := sub("_genomic$", "", Genome_ID)]
mast[, sp      := to_sp_code(species)]
mast[, carrier := to_pres(has_plasmid_corrected)]
stopifnot(all(mast$sp %in% c("bb", "gar", "afz", "bav")))
cat("Full dataset per species (expect bb 93 / gar 75 / afz 39 / bav 34):\n")
print(mast[, .N, by = sp][order(sp)])

wide <- dcast(all_long, genome ~ family, value.var = "present", fill = FALSE)
wide <- merge(mast[, .(genome, sp, carrier)], wide, by = "genome")
stopifnot(nrow(wide) == 241)
fam_cols <- intersect(expected, names(wide))
wide[, (fam_cols) := lapply(.SD, as.numeric), .SDcols = fam_cols]

w_carriers <- wide[carrier == TRUE]
cat("Carrier subset n =", nrow(w_carriers), " | full dataset n =", nrow(wide), "\n")
cat("Carrier subset per species (expect bb 91 / gar 16 / afz 6 / bav 22):\n")
print(w_carriers[, .N, by = sp][order(sp)])

fwrite(w_carriers, "results/antigen_pav_matrix_carriers.csv")
fwrite(wide,       "results/antigen_pav_matrix_all.csv")

## ---- 6) Per-species detection rates (carrier subset) ----
fmt_kn <- function(k, n) sprintf("%d/%d (%.1f%%)", k, n, 100 * k / n)
rates <- rbindlist(lapply(expected, function(f)
  w_carriers[, .(family = f, k = sum(get(f)), n = .N), by = sp]))
R <- dcast(rates, family ~ sp, value.var = c("k", "n"))
R[, `:=`(bb  = fmt_kn(k_bb,  n_bb),
         gar = fmt_kn(k_gar, n_gar),
         afz = fmt_kn(k_afz, n_afz),
         bav = fmt_kn(k_bav, n_bav))]
tot <- rates[, .(k = sum(k), n = sum(n)), by = family]
R <- merge(R, tot, by = "family")
R[, total := fmt_kn(k, n)]
R[, family := factor(family, levels = expected)]
setorder(R, family)
R_out <- R[, .(family, bb, gar, `afz†` = afz, bav, total)]
cat("\n== Detection rates, carrier subset, k/n (%) ==\n")
print(R_out)
fwrite(R_out, "results/antigen_detection_rates_by_species.csv")
cat("Saved: results/antigen_detection_rates_by_species.csv\n")
