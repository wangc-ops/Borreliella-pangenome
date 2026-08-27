########################################################################
## 07_vf_crossvalidation.R
## Cross-validation of antigen detection against VFDB virulence factor
## gene groups (VFGs) in the carrier subset, plus genome-level checks
## that unions of lineage-split VFGs reproduce the curated detection.
## Depends on: 01_build_pav_matrices.R (results/antigen_pav_matrix_carriers.csv)
## Statistics: 2 x 3 Fisher exact tests on bb/gar/bav with BH correction;
##             afz (n = 6) is display-only (daggered column).
## Inputs  (data/):    vf_detection_by_species_reconciled.csv
##                     vf_pav_matrix.csv
## Outputs (results/): vf_detection_by_species_carriers.csv
########################################################################
suppressPackageStartupMessages(library(data.table))

dir.create("results", showWarnings = FALSE)

vf_path  <- "data/vf_detection_by_species_reconciled.csv"
vfm_path <- "data/vf_pav_matrix.csv"
pav_path <- "results/antigen_pav_matrix_carriers.csv"
out_csv  <- "results/vf_detection_by_species_carriers.csv"

## ---- Part 1: per-species VFG detection table with tests ----
vf <- fread(vf_path)
stopifnot(all(c("vfg_id","full","full_n","subset","subset_n","Gene","Annotation") %in% names(vf)))
map <- c("91" = "bb", "16" = "gar", "6" = "afz", "22" = "bav")
key <- as.character(vf$subset_n)
stopifnot(all(key %in% names(map)))
vf[, sp := unname(map[key])]
chk <- vf[, .(n = .N, nsp = uniqueN(sp)), by = vfg_id]
stopifnot(all(chk$n == 4), all(chk$nsp == 4))
exp_full <- c(bb = 93, gar = 75, afz = 39, bav = 34)
vf[, bad := full_n != unname(exp_full[sp])]
stopifnot(sum(vf$bad) == 0, all(vf$subset <= vf$full))
vf[, bad := NULL]
cat("Checks passed:", uniqueN(vf$vfg_id), "VFGs x 4 species (carrier-subset scale)\n")

test_one <- function(d) {
  d3 <- d[sp %in% c("bb", "gar", "bav")]
  if (all(d3$subset == d3$subset_n)) return(list(p = NA_real_, cat0 = "invariant_universal"))
  if (all(d3$subset == 0))           return(list(p = NA_real_, cat0 = "invariant_absent"))
  mat <- matrix(c(d3$subset, d3$subset_n - d3$subset), nrow = 2, byrow = TRUE)
  list(p = tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_),
       cat0 = "tested")
}
res <- vf[, test_one(.SD), by = .(vfg_id, Gene, Annotation)]
res[, p_adj := p.adjust(p, method = "BH")]
res[, category := fifelse(cat0 != "tested", cat0,
                          fifelse(!is.na(p_adj) & p_adj < 0.05, "species_biased", "intermediate"))]

fmt <- function(k, n) sprintf("%d/%d (%.1f%%)", k, n, 100 * k / n)
wide <- dcast(vf, vfg_id ~ sp, value.var = c("subset", "subset_n"))
wide[, `:=`(bb_str  = fmt(subset_bb,  subset_n_bb),
            gar_str = fmt(subset_gar, subset_n_gar),
            afz_str = fmt(subset_afz, subset_n_afz),
            bav_str = fmt(subset_bav, subset_n_bav))]
VFT <- merge(res, wide[, .(vfg_id, bb_str, gar_str, afz_str, bav_str)], by = "vfg_id")
setnames(VFT, c("bb_str","gar_str","afz_str","bav_str"),
         c("bb (k/n, %)","gar (k/n, %)","afz† (k/n, %)","bav (k/n, %)"))
setorderv(VFT, c("category", "p_adj"), na.last = TRUE)
fwrite(VFT[, .(vfg_id, Gene, Annotation,
               `bb (k/n, %)`, `gar (k/n, %)`, `afz† (k/n, %)`, `bav (k/n, %)`,
               fisher_p = signif(p, 3), p_adj_BH = signif(p_adj, 3), category)],
       out_csv)
cat("Saved:", out_csv, "\n")
cat("\n== Category summary ==\n"); print(VFT[, .N, by = category][order(-N)])
cat("\n== Genes of interest ==\n")
for (g in c("cspZ", "cspA|CRASP", "dbpB", "ospB", "ospC", "bbk32", "guaA|GMP", "guaB|IMP")) {
  hit <- VFT[grepl(g, Gene, ignore.case = TRUE) | grepl(g, Annotation, ignore.case = TRUE)]
  if (nrow(hit) == 0) { cat(g, ": not found\n"); next }
  cat("-- ", g, " --\n", sep = "")
  print(hit[, .(vfg_id, Gene, `bb (k/n, %)`, `gar (k/n, %)`,
                `afz† (k/n, %)`, `bav (k/n, %)`, p_adj, category)])
}

## ---- Part 2: genome-level agreement, VFG unions vs curated detection ----
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

vfm <- fread(vfm_path)
pav <- fread(pav_path)
spcol <- intersect(names(pav), c("sp", "species"))[1]
setnames(pav, spcol, "sp")
pav[, sp := to_sp_code(sp)]

targets <- list(OspB = c("ospB", "group_1304"),
                DbpB = c("dbpB", "group_1337"),
                OspC = c("ospC", "group_1568"))
meta_cols <- c("Gene","Non-unique Gene name","Annotation",
               "vfg_id","VF_Name","VFcategory","Function","pident","evalue")
genome_cols <- setdiff(names(vfm), meta_cols)
cat("\nUnion check: number of genome columns =", length(genome_cols), "(expect 241)\n")
long <- melt(vfm, id.vars = "Gene", measure.vars = genome_cols,
             variable.name = "genome", value.name = "val")
long[, genome := sub("_genomic$", "", genome)]
long[, val := as.numeric(suppressWarnings(as.numeric(val)) > 0)]

for (lc in names(targets)) {
  hit <- long[Gene %in% targets[[lc]]]
  cat("\n== ", lc, " ==  union of: ", paste(unique(hit$Gene), collapse = " + "), "\n", sep = "")
  uni <- hit[, .(u = as.numeric(any(val == 1))), by = genome]
  cur <- pav[, .(genome, sp, curated = as.numeric(get(lc)))]
  m   <- merge(cur, uni, by = "genome", all.x = TRUE)
  m[is.na(u), u := 0]
  cat("Genome-level agreement:", sum(m$u == m$curated), "/", nrow(m), "\n")
  print(m[, .(union_pos = sum(u), curated_pos = sum(curated), n = .N), by = sp])
  disc <- m[u != curated]
  if (nrow(disc) > 0) { cat("Discordant genomes:\n"); print(disc) }
}
