# ============================================================
# 01. Antigen 13-family gene-level table (foundation for all downstream analyses)
# Input : data/02_antigen_gene_data.tsv          (legacy antigen gene table, 11 families)
#         data/antigen_pav_final_241.csv         (final 13-family PAV, reference "fin")
#         data/antigen_copy_numbers.csv          (final copy numbers, reference "cop")
#         external: GFF3 directory (241 *_genomic.gff3)  -> EXT_GFF_DIR
#         external: gene_presence_absence_roary.csv (Panaroo) -> EXT_PAV
#         external: gene_data.csv (Panaroo) -> EXT_GENE_DATA
# Output: results/gff_contigs.rds / gff_cds.rds  (contig/CDS caches, shared by later scripts)
#         results/ag11_clean.rds                 (cleaned legacy 11-family table cache)
#         results/antigen13_genes.csv            (4,449 rows: genome, locus, contig, fam13, pfam, source)
#         results/antigen13_family_map.csv       (fam13 x Panaroo family listing)
# Note  : frozen conventions (all validated; do not modify):
#         - replicon rule: longest contig per genome = chromosome; no length
#           threshold for gene-level localization
#         - DbpA = union of 6 allelic mutually exclusive families
#           (dbpA + group_1329/1198/1175/1276/1298)
#         - OspD = the single ospD family (45/45 identical to fin;
#           group_107/group_755 are ospD-like, excluded; see M&M QC note)
#         - Erp = 65 families grepped "erp|crasp" (241/241 copy-number match; includes cspZ)
#         - Mlp = 23 families grepped "Mlp" minus group_436 (contaminant family)
#         - roary cells with multiple loci are split by ";" (not whitespace)
#         - 935T (GCF_000714705.1) ospB is an assembly/annotation artifact, removed
# ============================================================
suppressPackageStartupMessages(library(data.table))

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

# External inputs (too large to ship in data/; set these paths before running)
EXT_GFF_DIR   <- "external/gff3_four_species"
EXT_PAV       <- "external/gene_presence_absence_roary.csv"
EXT_GENE_DATA <- "external/gene_data.csv"
stopifnot(dir.exists(EXT_GFF_DIR), file.exists(EXT_PAV), file.exists(EXT_GENE_DATA))

## ============================================================
## Part 1: parse 241 GFF3 files -> contig / CDS caches
## ============================================================
fs <- list.files(EXT_GFF_DIR, pattern = "\\.gff3$", full.names = TRUE)
stopifnot(length(fs) == 241)

parse_one <- function(f) {
  gid <- sub("_genomic\\.gff3$", "", basename(f))
  ln  <- readLines(f)
  sr  <- grep("^##sequence-region", ln, value = TRUE)
  pr  <- strsplit(sub("^##sequence-region[ ]*", "", sr), "\\s+")
  ctg <- data.table(genome = gid,
                    seqid  = sapply(pr, `[`, 2),
                    len    = suppressWarnings(as.numeric(sapply(pr, `[`, 4))))
  body <- ln[!grepl("^#", ln) & nzchar(ln)]
  fld  <- strsplit(body, "\t", fixed = TRUE)
  cds  <- fld[sapply(fld, function(x) length(x) >= 9 && x[3] == "CDS")]
  list(ctg = ctg,
       cds = data.table(genome  = gid,
                        seqid   = sapply(cds, `[`, 1),
                        locus   = sub(".*locus_tag=([^;]+).*", "\\1", sapply(cds, `[`, 9)),
                        product = sub(".*product=([^;]+).*",  "\\1", sapply(cds, `[`, 9))))
}
res     <- lapply(fs, parse_one)
contigs <- rbindlist(lapply(res, `[[`, "ctg"))
cds     <- rbindlist(lapply(res, `[[`, "cds"))

# Replicon rule: longest contig = chromosome (all 241 genomes have a single
# chromosome contig, 852-923 kb, zero exceptions)
stopifnot(!any(is.na(contigs$len)))
contigs[, replicon := { r <- rep("plasmid", .N); r[which.max(len)] <- "chromosome"; r },
        by = genome]
cat("contigs:", nrow(contigs), " chromosomes:", sum(contigs$replicon == "chromosome"), "\n")
cat("chromosome length range:", paste(range(contigs[replicon == "chromosome", len]), collapse = " - "), "\n")
stopifnot(nrow(contigs) == 3292,
          sum(contigs$replicon == "chromosome") == 241,
          min(contigs[replicon == "chromosome", len]) == 852383,
          max(contigs[replicon == "chromosome", len]) == 922801)

saveRDS(contigs, file.path(DIR_RES, "gff_contigs.rds"))
saveRDS(cds,     file.path(DIR_RES, "gff_cds.rds"))

## ============================================================
## Part 2: legacy gene table -> cleaned 11-family table (ag11_clean)
##   OspA_B split into OspA/OspB; DbpA_B into DbpA/DbpB; Erp_OspE -> Erp
##   (which() indices keep alignment, avoiding logical-vector recycling;
##    zero residuals after splitting)
## ============================================================
ag <- fread(file.path(DIR_DATA, "02_antigen_gene_data.tsv"))
ag[, genome := sub("_genomic$", "", genome_id)]
ag[, fam13 := family]
gn <- tolower(ag$gene_name); de <- tolower(ag$description)
set_fam <- function(idx, pat, label) {
  hit <- idx[grepl(pat, gn[idx]) | grepl(pat, de[idx])]
  ag$fam13[hit] <<- label
  setdiff(idx, hit)
}
ia  <- which(ag$family == "OspA_B")
ia  <- set_fam(ia,  "ospa|surface protein a", "OspA")
ia  <- set_fam(ia,  "ospb|surface protein b", "OspB")
idb <- which(ag$family == "DbpA_B")
idb <- set_fam(idb, "dbpa|decorin.binding protein a", "DbpA")
idb <- set_fam(idb, "dbpb|decorin.binding protein b", "DbpB")
ag$fam13[ag$family == "Erp_OspE"] <- "Erp"
cat("split residuals: OspA_B", length(ia), " DbpA_B", length(idb), "\n")
stopifnot(length(ia) == 0, length(idb) == 0)
saveRDS(ag, file.path(DIR_RES, "ag11_clean.rds"))

## ============================================================
## Part 3: reference tables + PAV matrix + family lists
## ============================================================
fin <- fread(file.path(DIR_DATA, "antigen_pav_final_241.csv"))
fin[, genome := sub("_genomic$", "", genome)]
cop <- fread(file.path(DIR_DATA, "antigen_copy_numbers.csv"))

pav  <- fread(EXT_PAV)
gcol <- grep("_genomic$", names(pav), value = TRUE)
gn_  <- sub("_genomic$", "", gcol)
stopifnot(length(gcol) == 241, all(fin$genome %in% gn_))
mat  <- as.matrix(pav[, ..gcol])          # raw locus strings; no binarization
rownames(mat) <- pav$Gene

fam_list <- list(
  DbpA = c("dbpA", "group_1329", "group_1198", "group_1175", "group_1276", "group_1298"),
  OspD = "ospD",
  Erp  = pav[grepl("erp|crasp", Annotation, ignore.case = TRUE), Gene],
  Mlp  = setdiff(pav[grepl("Mlp", Annotation, ignore.case = TRUE), Gene], "group_436")
)
cat("Erp families:", length(fam_list$Erp), " Mlp families:", length(fam_list$Mlp), "\n")
stopifnot(length(fam_list$Erp) == 65, length(fam_list$Mlp) == 23)

## ============================================================
## Part 4: dual-channel gene-level table
##   channel A (oldtable): 9 families covered by the legacy table
##     (DbpA covered only 100/134 genomes and Erp 99/125 in the legacy
##      table, both incomplete -> channel B)
##   channel B (roary):    DbpA / OspD / Erp / Mlp, loci split by ";"
## ============================================================
keep9 <- c("OspA", "OspB", "OspC", "DbpB", "CspA_Z", "BBK32", "VlsE", "BmpA_B", "P66")
tabA <- ag[fam13 %in% keep9,
           .(genome, locus = annotation_id, contig = scaffold_name,
             fam13, pfam = NA_character_, source = "oldtable")]

extract_loci <- function(fams, label) {
  out <- list()
  for (f in fams) for (j in seq_along(gcol)) {
    x <- mat[f, j]
    if (is.na(x) || x == "" || x == "0") next
    ls <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
    out[[length(out) + 1]] <- data.table(genome = gn_[j], locus = ls,
                                         fam13 = label, pfam = f, source = "roary")
  }
  rbindlist(out)
}
tabB <- rbindlist(lapply(names(fam_list), function(lb) extract_loci(fam_list[[lb]], lb)))
tabB[, contig := NA_character_]

## ---- contig backfill, two channels: CDS channel first, gene_data fallback ----
cds_key <- cds[, .(locus, contig = seqid, genome)]
setkey(cds_key, locus, genome)
tabB[, contig := cds_key[.(tabB$locus, tabB$genome), contig]]
cat("unfilled after CDS channel:", sum(is.na(tabB$contig)), "/", nrow(tabB), "\n")
stopifnot(sum(is.na(tabB$contig)) == 68)

gd <- fread(EXT_GENE_DATA)
gd[, genome := sub("_genomic(\\.gff3)?$", "", basename(gff_file))]
gd_key <- unique(gd[, .(locus = annotation_id, genome, contig = scaffold_name)],
                 by = c("locus", "genome"))
setkey(gd_key, locus, genome)
ui <- which(is.na(tabB$contig))
tabB[ui, contig := gd_key[.(tabB$locus[ui], tabB$genome[ui]), contig]]
cat("still unfilled after both channels:", sum(is.na(tabB$contig)), "/", nrow(tabB), "\n")
stopifnot(!any(is.na(tabB$contig)))

## ============================================================
## Part 5: merge + 935T OspB removal
##   ospB of GCF_000714705.1_935T (GDDIIJ_01087, contig_7) is a
##   confirmed assembly/annotation artifact: the strain is not in the
##   OspB-positive set of fin
## ============================================================
tab <- rbind(tabA, tabB[, .(genome, locus, contig, fam13, pfam, source)])
extra <- setdiff(unique(tab[fam13 == "OspB", genome]), fin$genome[as.numeric(fin$OspB) == 1])
cat("extra OspB genomes:", paste(extra, collapse = ", "), "(expected: 935T only)\n")
stopifnot(identical(extra, "GCF_000714705.1_935T"))
n_before <- nrow(tab)
tab <- tab[!(fam13 == "OspB" & genome %in% extra)]
stopifnot(n_before - nrow(tab) == 1L)

## ============================================================
## Part 6: hard assertions -- 13-family carrier sets == fin
## ============================================================
fams13 <- c("OspA","OspB","OspC","OspD","DbpA","DbpB","Erp","Mlp",
            "CspA_Z","BBK32","VlsE","BmpA_B","P66")
expect_n <- c(OspA=134, OspB=117, OspC=135, OspD=45, DbpA=134, DbpB=134,
              Erp=125, Mlp=125, CspA_Z=17, BBK32=16, VlsE=16, BmpA_B=240, P66=241)
ok_all <- TRUE
for (f in fams13) {
  carriers <- unique(tab[fam13 == f, genome])
  pos  <- fin$genome[as.numeric(fin[[f]]) == 1]
  same <- identical(sort(carriers), sort(pos))
  ok_all <- ok_all && same && length(carriers) == expect_n[[f]]
  cat(sprintf("%-8s %3d vs fin %3d  %s\n", f, length(carriers), length(pos),
              ifelse(same, "✓", "✗✗✗")))
}
stopifnot(ok_all)

## ---- Erp/Mlp copy numbers == cop (241/241) ----
for (f in c("Erp", "Mlp")) {
  cnt <- tab[fam13 == f, .N, by = genome]
  ref <- cop[match(fin$genome, Genome_ID), if (f == "Erp") erp_copies else mlp_copies]
  v <- integer(241); names(v) <- fin$genome; v[cnt$genome] <- cnt$N
  cat(f, "copy-number match:", sum(v == ref), "/241\n")
  stopifnot(all(v == ref))
}

## ---- cross-family shared loci: 15, all CspA_Z x Erp(group_73); kept, not removed ----
dup <- tab[, .(nfam = uniqueN(fam13)), by = .(genome, locus)][nfam > 1]
cat("cross-family shared loci:", nrow(dup), "(expected 15, CspA_Z x Erp group_73)\n")
if (nrow(dup)) {
  shared <- tab[dup, on = .(genome, locus)][, .N, by = .(fam13, pfam)]
  print(shared)
  stopifnot(nrow(dup) == 15,
            nrow(shared[fam13 == "CspA_Z"]) == 1, nrow(shared[fam13 == "Erp"]) == 1,
            shared[fam13 == "Erp", pfam] == "group_73")
}

## ============================================================
## Part 7: export
## ============================================================
stopifnot(nrow(tab) == 4449, uniqueN(tab$fam13) == 13)
fwrite(tab, file.path(DIR_RES, "antigen13_genes.csv"))
fmap <- rbindlist(lapply(names(fam_list), function(lb)
  data.table(fam13 = lb, panaroo_family = fam_list[[lb]])))
fwrite(fmap, file.path(DIR_RES, "antigen13_family_map.csv"))

loc_tab <- tab[, .(loci = .N, genomes = uniqueN(genome)), by = fam13]
print(loc_tab[order(-loci)])
expect_loci <- c(Erp=1847L, BmpA_B=817L, Mlp=769L, P66=241L, OspC=136L, OspA=134L,
                 DbpB=134L, DbpA=134L, OspB=117L, OspD=45L, VlsE=42L, CspA_Z=17L, BBK32=16L)
got_loci <- setNames(loc_tab$loci, loc_tab$fam13)[names(expect_loci)]
stopifnot(all(got_loci == expect_loci))
cat("=== 01 done: antigen13_genes.csv (4,449 rows) + family map exported; all assertions passed ===\n")
