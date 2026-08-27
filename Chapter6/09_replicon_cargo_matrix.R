# ============================================================
# 09. Table S20: cargo-replicon coupling matrix + co-carriage tests
# Input : data/gff_cds.rds / gff_contigs.rds        (from script 01)
#         data/antigen13_genes.csv                  (from script 01)
#         data/Table_S18_plasmid_antigen_function_matrix_v2.csv  (from script 07)
#         data/strain_plasmid_status_master_corrected.csv
#         external: gene_presence_absence_roary.csv -> EXT_PAV
# Output: results/tableS20_replicon_cargo_matrix.csv  (matrix)
#         results/replicon_cargo_tests.csv            (paired Wilcoxon + Fisher,
#                                                      with >=10 kb sensitivity)
#         results/replicon_cargo_tests_family.csv     (family-deduplicated convention)
# Note  : replicon types are anchored by antigen marker genes; hitchhiking
#         paralogs on mixed contigs are resolved by a length-priority rule;
#         cargo composition goes locus -> PAV family -> S18 v2 category with
#         Antigen priority override (same convention as Fig 6a)
# ============================================================
suppressPackageStartupMessages(library(data.table))

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

# External input (too large to ship in data/; set the path before running)
EXT_PAV <- "external/gene_presence_absence_roary.csv"
stopifnot(file.exists(EXT_PAV))

cds   <- readRDS(file.path(DIR_DATA, "gff_cds.rds"))       # genome, seqid, locus, product
ct    <- readRDS(file.path(DIR_DATA, "gff_contigs.rds"))   # genome, seqid, len, replicon
a13   <- fread(file.path(DIR_DATA, "antigen13_genes.csv")) # genome, locus, contig, fam13
s18v2 <- fread(file.path(DIR_DATA, "Table_S18_plasmid_antigen_function_matrix_v2.csv"))

# ---------- 1. locus -> PAV family ----------
raw <- fread(EXT_PAV)
meta_cols <- c("Gene","Non-unique Gene name","Annotation","No. isolates",
               "No. sequences","Avg sequences per isolate","Genome Fragment",
               "Order within Fragment","Accessory Fragment","Accessory Order with Fragment",
               "QC","Min group size nuc","Max group size nuc","Avg group size nuc")
genome_cols <- setdiff(names(raw), meta_cols)
stopifnot(length(genome_cols) == 241, all(endsWith(genome_cols, "_genomic")))
long <- melt(raw, id.vars = "Gene", measure.vars = genome_cols,
             variable.name = "genome", value.name = "cell")
long[, genome := sub("_genomic$", "", as.character(genome))]
long <- long[!is.na(cell) & cell != ""]
fam_map <- long[, .(locus = trimws(unlist(strsplit(cell, "[;\t]")))),
                by = .(Gene_Family = Gene, genome)]

# ---------- 2. attach category to every gene (Antigen priority override,
#   same convention as Fig 6a) ----------
ag_fams <- s18v2[Antigen == 1, Gene_Family]
cds2 <- merge(cds, fam_map, by = c("genome", "locus"), all.x = TRUE)
cat("CDS family match rate:", sprintf("%.1f%%", 100 * mean(!is.na(cds2$Gene_Family))), "\n")
cds2 <- merge(cds2, s18v2[, .(Gene_Family, category = Function_Category)],
              by = "Gene_Family", all.x = TRUE)
cds2[Gene_Family %in% ag_fams, category := "Antigen"]
cds2 <- cds2[!is.na(category)]
stopifnot(all(cds2$category %in% c("Metabolic/Other", "Surface-associated",
                                   "Hypothetical", "MGE/Plasmid", "Antigen")))

# ---------- 3. replicon anchoring + length-priority resolution of mixed types ----------
type_vec <- c(Erp = "cp32-like", Mlp = "cp32-like",
              OspA = "lp54-like", OspB = "lp54-like",
              DbpA = "lp54-like", DbpB = "lp54-like",
              OspC = "cp26", VlsE = "lp28-1-like",
              BBK32 = "lp36-like", OspD = "lp38-like",
              CspA_Z = "other-Ag", P66 = "chromosomal", BmpA_B = "chromosomal")
ag_ctg <- a13[, .(types = list(unique(type_vec[fam13]))), by = .(genome, contig)]
ag_ctg <- merge(ag_ctg, ct, by.x = c("genome", "contig"), by.y = c("genome", "seqid"))
stopifnot(nrow(ag_ctg) > 0)

resolve_type <- function(t, len) {
  t2 <- setdiff(t, "chromosomal")
  if (length(t2) == 0) return("chromosomal")
  if (length(t2) == 1) return(t2)
  if ("lp54-like"   %in% t2 && len > 40000) return("lp54-like")   # lp54 ~54 kb
  if ("lp28-1-like" %in% t2 && len < 40000) return("lp28-1-like") # vlsE passenger spillover
  if ("cp32-like"   %in% t2) return("cp32-like")                  # erp/mlp bulk
  paste(sort(t2), collapse = "+")
}
ag_ctg[, replicon_type_raw := sapply(types, function(t) {
  t2 <- setdiff(t, "chromosomal")
  if (length(t2) <= 1) return(if (length(t2) == 0) "chromosomal" else t2)
  "mixed"
})]
ag_ctg[, replicon_type := mapply(resolve_type, types, len)]
cat("mixed-type resolution, before -> after:\n")
print(ag_ctg[, .N, by = .(replicon_type_raw, replicon_type)][order(replicon_type_raw, -N)])

chr_ok <- ag_ctg[replicon_type == "chromosomal", mean(replicon == "chromosome")]
pla_ok <- ag_ctg[replicon_type != "chromosomal", mean(replicon == "plasmid")]
cat(sprintf("localization concordance: chromosomal antigens -> chromosome %.1f%% | others -> plasmid %.1f%%\n",
            100 * chr_ok, 100 * pla_ok))
stopifnot(chr_ok > 0.95, pla_ok > 0.90)

# B31 benchmark (GCF_000008685.2 = B. burgdorferi B31 complete genome;
# non-empty assertions)
b31 <- ag_ctg[genome == "GCF_000008685.2_ASM868v2", .N, by = replicon_type]
get_n <- function(tp) b31[replicon_type == tp, sum(N)]
cat("B31 replicon types:\n"); print(b31)
stopifnot(get_n("cp26") == 1, get_n("lp54-like") == 1,
          get_n("lp28-1-like") == 1,
          get_n("cp32-like") >= 4 & get_n("cp32-like") <= 14,
          get_n("lp36-like") == 1, get_n("lp38-like") == 1)
cat("sanity: B31 benchmark passed\n")

# ---------- 4. cargo matrix: replicon type x functional category ----------
ag_key <- ag_ctg[replicon == "plasmid", .(genome, contig, replicon_type)]
cargo <- merge(cds2, ag_key, by.x = c("genome", "seqid"), by.y = c("genome", "contig"))
mat_c <- dcast(cargo[, .N, by = .(replicon_type, category)],
               replicon_type ~ category, fill = 0)
cat("\ncargo matrix (gene counts):\n"); print(mat_c)
# Final anchors (main-text run values): antigen genes cp32-like 2,497 /
#   lp54-like 647 / cp26 136 / lp28-1-like + lp38-like + lp36-like total 116
stopifnot(mat_c[replicon_type == "cp32-like", Antigen] == 2497,
          mat_c[replicon_type == "lp54-like", Antigen] == 647,
          mat_c[replicon_type == "cp26", Antigen] == 136,
          mat_c[replicon_type %in% c("lp28-1-like", "lp38-like", "lp36-like"),
                sum(Antigen)] == 116)
fwrite(mat_c, file.path(DIR_RES, "tableS20_replicon_cargo_matrix.csv"))

# ---------- 5. co-carriage tests (bb/gar/bav, n=129) ----------
master <- fread(file.path(DIR_DATA, "strain_plasmid_status_master_corrected.csv"))
spl <- tolower(master$species)
master[, sp := fcase(grepl("garinii", spl), "gar", grepl("afzelii", spl), "afz",
                     grepl("bavariensis", spl), "bav", grepl("burgdorferi", spl), "bb")]
keep_g <- master[sp %in% c("bb", "gar", "bav"), Genome_ID]

pla_ctg <- ct[replicon == "plasmid" & genome %in% keep_g, .(genome, seqid, len)]
pla_ctg[, is_ag := paste(genome, seqid) %in% paste(ag_key$genome, ag_key$contig)]
cmp <- merge(cds2, pla_ctg, by.x = c("genome", "seqid"), by.y = c("genome", "seqid"))

# 5a. gene-copy (CDS) convention
run_tests <- function(dt, tag) {
  pg <- dt[, .(hyp_ag = sum(category == "Hypothetical" & is_ag),
               tot_ag = sum(is_ag),
               hyp_na = sum(category == "Hypothetical" & !is_ag),
               tot_na = sum(!is_ag)), by = genome]
  pg <- pg[tot_ag > 0 & tot_na > 0]   # genomes need both contig classes for pairing
  pg[, `:=`(share_ag = hyp_ag / tot_ag, share_na = hyp_na / tot_na)]
  wt <- suppressWarnings(wilcox.test(pg$share_ag, pg$share_na, paired = TRUE))
  ft <- fisher.test(matrix(c(pg[, sum(hyp_ag)], pg[, sum(tot_ag - hyp_ag)],
                             pg[, sum(hyp_na)], pg[, sum(tot_na - hyp_na)]),
                           nrow = 2, byrow = TRUE))
  data.table(test = tag, n_genomes = nrow(pg),
             median_share_ag = median(pg$share_ag),
             median_share_na = median(pg$share_na),
             wilcox_P = wt$p.value,
             fisher_OR = unname(ft$estimate), fisher_P = ft$p.value)
}
res <- rbind(run_tests(cmp, "CDS-level, all plasmid contigs"),
             run_tests(cmp[len >= 10000], "CDS-level, contigs >= 10 kb"))
print(res)
fwrite(res, file.path(DIR_RES, "replicon_cargo_tests.csv"))

# 5b. gene-family (deduplicated) convention -- excludes dilution by
# multicopy antigen arrays
run_tests_fam <- function(dt, tag) {
  pg <- dt[, .(hyp_ag = uniqueN(Gene_Family[category == "Hypothetical" & is_ag]),
               fam_ag = uniqueN(Gene_Family[is_ag]),
               hyp_na = uniqueN(Gene_Family[category == "Hypothetical" & !is_ag]),
               fam_na = uniqueN(Gene_Family[!is_ag])), by = genome]
  pg <- pg[fam_ag > 0 & fam_na > 0]
  pg[, `:=`(share_ag = hyp_ag / fam_ag, share_na = hyp_na / fam_na)]
  wt <- suppressWarnings(wilcox.test(pg$share_ag, pg$share_na, paired = TRUE))
  data.table(test = tag, n_genomes = nrow(pg),
             median_share_ag = median(pg$share_ag),
             median_share_na = median(pg$share_na),
             wilcox_P = wt$p.value)
}
res_fam <- rbind(run_tests_fam(cmp, "family-level, all plasmid contigs"),
                 run_tests_fam(cmp[len >= 10000], "family-level, contigs >= 10 kb"))
print(res_fam)
fwrite(res_fam, file.path(DIR_RES, "replicon_cargo_tests_family.csv"))

# Hard assertions: direction/significance + final value locks (main-text run values)
#   family convention: all contigs 10.2% vs 22.0% (n=124); >=10 kb 10.7% vs 28.6%
#   CDS convention:    all contigs 8.6% vs 19.0%
stopifnot(res$median_share_ag < res$median_share_na,
          all(res$wilcox_P < 1e-10), all(res$fisher_OR < 0.5),
          all(res_fam$median_share_ag < res_fam$median_share_na),
          all(res_fam$wilcox_P < 1e-10))
stopifnot(res_fam[1, n_genomes] == 124,
          abs(res_fam[1, median_share_ag] - 0.102) < 0.0005,
          abs(res_fam[1, median_share_na] - 0.220) < 0.0005,
          abs(res_fam[2, median_share_ag] - 0.107) < 0.0005,
          abs(res_fam[2, median_share_na] - 0.286) < 0.0005,
          abs(res[1, median_share_ag] - 0.086) < 0.0005,
          abs(res[1, median_share_na] - 0.190) < 0.0005)
cat("sanity: co-carriage tests (CDS + family conventions) direction, significance and values all locked\n")
cat("=== 09 done: Table S20 + co-carriage tests exported ===\n")
