# ============================================================
# 03. Family-level plasmid association + Table S18 flag audit + tier statistics
#     (source of the main-text association sentence)
# Input : data/Table_S18_plasmid_antigen_function_matrix.csv (flag column audited)
#         data/gff_contigs.rds              (from script 01)
#         external: gene_data.csv                   -> EXT_GENE_DATA
#         external: gene_presence_absence_roary.csv -> EXT_PAV
# Output: results/s18_flag_audit.csv           (audit comparison + final flags)
#         results/tier_plasmid_association.csv (tier x association, two tier conventions)
# Note  : background -- plasmid-borne families are structurally capped at
#         135/241 = 56.0% presence (106 bare-chromosome genomes reflect a
#         deposition bias), so no core/soft-core family can be plasmid-associated;
#         the tier x association result is reported as a main-text sentence
#         rather than a panel.
#         Declared rule (M&M): a family is plasmid-associated when >50% of its
#         genes are plasmid-borne; thresholds 1/3, 1/2, 2/3 all give 94.8%.
#         Final flags = original S18 flags (no re-classification); the two
#         audit-discordant families (hsdR~~~parA 0.294, group_39 0.067) keep
#         their original flags, documented in s18_flag_audit.csv.
# ============================================================
suppressPackageStartupMessages(library(data.table))

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

# External inputs (too large to ship in data/; set these paths before running)
EXT_GENE_DATA <- "external/gene_data.csv"
EXT_PAV       <- "external/gene_presence_absence_roary.csv"
stopifnot(file.exists(EXT_GENE_DATA), file.exists(EXT_PAV))

s18 <- fread(file.path(DIR_DATA, "Table_S18_plasmid_antigen_function_matrix.csv"))
stopifnot(nrow(s18) == 2265)

## ============================================================
## Part 1: gene -> family -> location recomputation across all 241 genomes
## ============================================================
ct <- readRDS(file.path(DIR_DATA, "gff_contigs.rds"))
ct[, genome := sub("_genomic$", "", genome)]
chrom <- ct[, .SD[which.max(len)], by = genome][, .(genome, chrom_contig = seqid)]

gd <- fread(EXT_GENE_DATA)
gd[, genome := sub("_genomic(\\.gff3)?$", "", basename(gff_file))]
genes <- merge(gd[, .(genome, locus = annotation_id, contig = scaffold_name)],
               chrom, by = "genome")
genes[, location := fifelse(contig == chrom_contig, "Chromosome", "Plasmid")]

pav  <- fread(EXT_PAV)
gcol <- grep("_genomic$", names(pav), value = TRUE)
stopifnot(length(gcol) == 241)
rl <- melt(pav[, c("Gene", gcol), with = FALSE], id.vars = "Gene",
           variable.name = "genome_g", value.name = "cell")
rl <- rl[cell != "" & cell != "0" & !is.na(cell)]
l2f <- rl[, .(locus = trimws(unlist(strsplit(cell, ";", fixed = TRUE)))),
          by = .(family = Gene, genome = sub("_genomic$", "", genome_g))]
setkey(l2f, locus, genome)
genes[, family := l2f[.(genes$locus, genes$genome), family]]

gf <- genes[!is.na(family),
            .(n_pl = sum(location == "Plasmid"),
              n_ch = sum(location == "Chromosome")), by = family]
gf[, pl_frac := round(n_pl / (n_pl + n_ch), 3)]
gf[, assoc_ours := fifelse(pl_frac > 0.5, 1L, 0L)]

## ============================================================
## Part 2: S18 flag audit (concordance 1425/1427)
## ============================================================
cmp <- merge(s18[, .(Gene_Family, N_Genomes_Present, Plasmid_Associated)],
             gf, by.x = "Gene_Family", by.y = "family")
audit_tab <- table(S18 = cmp$Plasmid_Associated, ours = cmp$assoc_ours)
cat("\n== flag comparison (S18 rows x recomputed columns) ==\n"); print(audit_tab)
stopifnot(audit_tab["0", "0"] == 838, audit_tab["0", "1"] == 0,
          audit_tab["1", "0"] == 2,   audit_tab["1", "1"] == 1425)

disc <- cmp[Plasmid_Associated != assoc_ours]
cat("discordant families (S18 = plasmid, gene-level recomputation = chromosome majority):\n")
print(disc[, .(Gene_Family, N_Genomes_Present, n_pl, n_ch, pl_frac)])
stopifnot(identical(sort(disc$Gene_Family), sort(c("hsdR~~~parA", "group_39"))),
          abs(disc[Gene_Family == "hsdR~~~parA", pl_frac] - 0.294) < 0.001,
          abs(disc[Gene_Family == "group_39", pl_frac] - 0.067) < 0.001)

# Final classification = original S18 flags (no re-classification);
# assoc_ours / pl_frac remain in the output as the audit record.
cmp[, Plasmid_Associated_final := as.integer(Plasmid_Associated)]
fwrite(cmp[, .(Gene_Family, N_Genomes_Present,
               Plasmid_Associated_S18 = as.integer(Plasmid_Associated),
               n_pl, n_ch, pl_frac, assoc_ours, Plasmid_Associated_final)],
       file.path(DIR_RES, "s18_flag_audit.csv"))
cat("final plasmid-associated families (original S18 flags):", sum(cmp$Plasmid_Associated_final),
    "; under the >50% gene-level rule:", sum(cmp$assoc_ours),
    "(difference of 2 = hsdR~~~parA, group_39; audit record, flags kept)\n")

## ============================================================
## Part 3: tier x association + main-text sentence (two tier conventions
##   reported side by side; convention B goes to the main text)
##   Convention A: S18 freq column (<99% = accessory)
##   Convention B: direct PAV counting (Core>=99% / Soft>=95% / Shell>=15% / Cloud)
##   With the original S18 flags, convention B = 1,427 / 66 / 95.6% (main text).
##   The threshold-sensitivity check uses the gene-level rule (pl_frac) and
##   differs from the flag-based count by exactly the 2 audited families;
##   this documented difference does not affect the "unchanged" conclusion.
## ============================================================
# Convention A
cmp[, freq := N_Genomes_Present / 241]
cmp[, tier_A := fcase(freq >= 0.99, "Core", freq >= 0.95, "Soft",
                      freq >= 0.15, "Shell", default = "Cloud")]
tabA <- cmp[, .(n = .N, plasmid = sum(Plasmid_Associated_final),
                pct = round(100 * sum(Plasmid_Associated_final) / .N, 1)),
            by = tier_A][order(-n)]
cat("\n== convention A (S18 freq): tier x plasmid association ==\n"); print(tabA)
accA <- cmp[tier_A != "Core"]
sentA <- sprintf("Convention A: %d accessory families, %d plasmid-associated = %.1f%%; %d chromosomal",
                 nrow(accA), sum(accA$Plasmid_Associated_final),
                 100 * sum(accA$Plasmid_Associated_final) / nrow(accA),
                 nrow(accA) - sum(accA$Plasmid_Associated_final))
cat(sentA, "\n")
stopifnot(nrow(accA) == 1503, sum(accA$Plasmid_Associated_final) == 1427)  # 94.9%

# Threshold sensitivity (convention-A accessory set, gene-level rule)
for (th in c(1/3, 0.5, 2/3)) {
  cat(sprintf("threshold >%.3f: %.1f%% ", th,
              100 * sum(accA$pl_frac > th) / nrow(accA)))
}
cat("(all three thresholds should give 94.8%)\n")
stopifnot(all(sapply(c(1/3, 0.5, 2/3), function(th)
  abs(100 * sum(accA$pl_frac > th) / nrow(accA) - 94.8) < 0.05)))

# Convention B: direct PAV counting (same matrix rule as the fig6d script)
raw <- as.matrix(pav[, ..gcol])
pav_bin <- t(ifelse(is.na(raw) | raw == "" | raw == "0", 0, 1))
freq_all <- colSums(pav_bin)
tier_B <- ifelse(freq_all >= 0.99 * 241, "Core",
          ifelse(freq_all >= 0.95 * 241, "Softcore",
          ifelse(freq_all >= 0.15 * 241, "Shell", "Cloud")))
stopifnot(sum(tier_B == "Core") == 763, sum(tier_B == "Softcore") == 9,
          sum(tier_B == "Shell") == 453, sum(tier_B == "Cloud") == 1040)
flag_final <- setNames(cmp$Plasmid_Associated_final, cmp$Gene_Family)
accB <- names(tier_B)[tier_B %in% c("Shell", "Cloud")]
plB  <- sum(flag_final[accB])
sentB <- sprintf("Convention B: %d Shell+Cloud families, %d plasmid-associated = %.1f%%; %d chromosomal",
                 length(accB), plB, 100 * plB / length(accB), length(accB) - plB)
cat("\n", sentB, "\n")
stopifnot(length(accB) == 1493, plB == 1427)  # 95.6%; 66 chromosomal (main-text convention)

# Structural-constraint sentence
cat("\nCore+Softcore plasmid-associated families (should be 0 under both): A =",
    cmp[tier_A %in% c("Core", "Soft"), sum(Plasmid_Associated_final)],
    " B =", sum(flag_final[names(tier_B)[tier_B %in% c("Core", "Softcore")]]), "\n")
cat("highest presence among plasmid-associated families:",
    round(100 * max(cmp[Plasmid_Associated_final == 1, N_Genomes_Present]) / 241, 1),
    "% (structural ceiling = 135/241 = 56.0%)\n")

fwrite(rbind(tabA[, .(tier = tier_A, n, plasmid, pct, scale = "A_S18freq")],
             data.table(tier = c("Shell+Cloud"), n = length(accB), plasmid = plB,
                        pct = round(100 * plB / length(accB), 1), scale = "B_PAV")),
       file.path(DIR_RES, "tier_plasmid_association.csv"))
cat("=== 03 done: audit table + tier x association exported ===\n")
