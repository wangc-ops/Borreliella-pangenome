# ============================================================
# 02. Plasmid vs chromosome functional composition + Fisher enrichment (Fig 6a)
#   Part 1: rebuild composition base table + gene-level Fisher (authoritative)
#   Part 2: mean-composition stacked bars (Fig 6a) + table-level Fisher cross-check
# Input : data/strain_plasmid_status_master_corrected.csv  (master)
#         data/antigen13_genes.csv                         (from script 01)
#         data/gff_contigs.rds                             (from script 01)
#         data/Table_S18_plasmid_antigen_function_matrix.csv (Function_Category column audited)
#         external: gene_data.csv                   -> EXT_GENE_DATA
#         external: gene_presence_absence_roary.csv -> EXT_PAV
# Output: results/fig6a_function_composition_135.csv  (genome, location, category, N, sp, total, prop)
#             NB: total repeats across the 5 category rows of a genome-location;
#                 unique by (genome, location) before summing across categories
#         results/fig6a_fisher_results.csv            (both 135 and 129 Fisher sets)
#         results/fig6a_function_composition_stacked_bar.pdf
# Note  : 135 plasmid-carrying genomes only; Panaroo-clustered genes only
#         (2,823 Bakta-annotated genes not clustered by Panaroo are excluded,
#          declared here: plasmid 2,748 / chromosome 75)
#         main-text Fisher set = 129 (bb/gar/bav); the 135 set is a robustness note
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

# External inputs (too large to ship in data/; set these paths before running)
EXT_GENE_DATA <- "external/gene_data.csv"
EXT_PAV       <- "external/gene_presence_absence_roary.csv"
stopifnot(file.exists(EXT_GENE_DATA), file.exists(EXT_PAV))

## ============================================================
## Part 1: base-table rebuild (established builder logic)
## ============================================================
mast <- fread(file.path(DIR_DATA, "strain_plasmid_status_master_corrected.csv"))
mast[, genome := sub("_genomic$", "", Genome_ID)]
mast[, carrier := as.logical(has_plasmid_corrected)]
sp_map <- c("Borreliella burgdorferi" = "bb", "Borreliella garinii" = "gar",
            "Borreliella afzelii" = "afz", "Borreliella bavariensis" = "bav")
mast[, sp := sp_map[species]]
stopifnot(!any(is.na(mast$sp)))
carriers <- mast[carrier == TRUE, genome]
stopifnot(length(carriers) == 135)

# Replicon rule: longest contig = chromosome (same as script 01; no length
# threshold at gene level)
ct <- readRDS(file.path(DIR_DATA, "gff_contigs.rds"))
ct[, genome := sub("_genomic$", "", genome)]
chrom <- ct[, .SD[which.max(len)], by = genome][, .(genome, chrom_contig = seqid)]

gd <- fread(EXT_GENE_DATA)
gd[, genome := sub("_genomic(\\.gff3)?$", "", basename(gff_file))]
genes <- merge(gd[, .(genome, locus = annotation_id, contig = scaffold_name)],
               chrom, by = "genome")
genes[, location := fifelse(contig == chrom_contig, "Chromosome", "Plasmid")]
g135 <- genes[genome %in% carriers]

# Antigen flag (matched within genome, guarding against cross-genome locus coincidence)
ag <- fread(file.path(DIR_DATA, "antigen13_genes.csv"))
g135[, is_antigen := paste(genome, locus) %in% paste(ag$genome, ag$locus)]

# locus -> Panaroo family (";" split)
s18 <- fread(file.path(DIR_DATA, "Table_S18_plasmid_antigen_function_matrix.csv"))
pav <- fread(EXT_PAV)
gcol <- grep("_genomic$", names(pav), value = TRUE)
rl <- melt(pav[, c("Gene", gcol), with = FALSE], id.vars = "Gene",
           variable.name = "genome_g", value.name = "cell")
rl <- rl[cell != "" & cell != "0" & !is.na(cell)]
l2f <- rl[, .(locus = trimws(unlist(strsplit(cell, ";", fixed = TRUE)))),
          by = .(family = Gene, genome = sub("_genomic$", "", genome_g))]
setkey(l2f, locus, genome)
g135[, family := l2f[.(g135$locus, g135$genome), family]]

# Drop Panaroo-unclustered genes (declared in M&M: 2,823; plasmid 2,748 / chromosome 75)
g135c <- g135[!is.na(family)]
cat("clustered genes:", nrow(g135c),
    " unclustered excluded:", nrow(g135) - nrow(g135c),
    "(plasmid", nrow(g135[is.na(family) & location == "Plasmid"]),
    "/chromosome", nrow(g135[is.na(family) & location == "Chromosome"]), ")\n")
stopifnot(nrow(g135c) == 162355,
          nrow(g135) - nrow(g135c) == 2823,
          nrow(g135[is.na(family) & location == "Plasmid"]) == 2748,
          nrow(g135[is.na(family) & location == "Chromosome"]) == 75)

g135c[, category := s18$Function_Category[match(family, s18$Gene_Family)]]
stopifnot(!any(is.na(g135c$category)))
g135c[is_antigen == TRUE, category := "Antigen"]   # Antigen priority override

## ---- Fisher enrichment (gene-level, authoritative) ----
run_fisher <- function(dt, tag) {
  rbindlist(lapply(sort(unique(dt$category)), function(cc) {
    a  <- dt[category == cc & location == "Plasmid", .N]
    b  <- dt[category == cc & location == "Chromosome", .N]
    c1 <- dt[category != cc & location == "Plasmid", .N]
    d  <- dt[category != cc & location == "Chromosome", .N]
    ft <- fisher.test(matrix(c(a, c1, b, d), 2, byrow = TRUE))
    data.table(set = tag, category = cc, plasmid_in = a, plasmid_out = c1,
               chrom_in = b, chrom_out = d,
               plasmid_pct = round(100 * a / (a + c1), 2),
               chrom_pct   = round(100 * b / (b + d), 2),
               OR = round(unname(ft$estimate), 2),
               CI = paste0(round(ft$conf.int, 2), collapse = "-"),
               P  = formatC(ft$p.value, format = "e", digits = 2))
  }))
}
g129 <- g135c[genome %in% mast[sp != "afz" & carrier == TRUE, genome]]
stopifnot(uniqueN(g129$genome) == 129)
f135 <- run_fisher(g135c, "135")
f129 <- run_fisher(g129,  "129_bb/gar/bav")
cat("\n== Fisher (135, robustness) ==\n");  print(f135)
cat("\n== Fisher (129, main text) ==\n");   print(f129)
# Locked anchors (129): Antigen 14.97 [13.59-16.50] | MGE 60.10 | Surface 2.34 | Hypo 2.11 | Metab 0.23
stopifnot(abs(f129[category == "Antigen", OR] - 14.97) < 0.01,
          abs(f129[category == "MGE/Plasmid", OR] - 60.10) < 0.01,
          abs(f129[category == "Surface-associated", OR] - 2.34) < 0.01,
          abs(f129[category == "Hypothetical", OR] - 2.11) < 0.01,
          abs(f129[category == "Metabolic/Other", OR] - 0.23) < 0.01,
          f129[category == "Antigen", plasmid_in] == 3263,
          f129[category == "Antigen", chrom_in] == 483)
# Anchors (135 set): Antigen OR 14.79, plasmid 6.43% vs chrom 0.46%
stopifnot(abs(f135[category == "Antigen", OR] - 14.79) < 0.01,
          abs(f135[category == "Antigen", plasmid_pct] - 6.43) < 0.01,
          abs(f135[category == "Antigen", chrom_pct] - 0.46) < 0.01)
fwrite(rbind(f135, f129), file.path(DIR_RES, "fig6a_fisher_results.csv"))

## ---- composition base table ----
comp <- g135c[, .N, by = .(genome, location, category)]
comp <- merge(comp, mast[, .(genome, sp)], by = "genome")
comp[, total := sum(N), by = .(genome, location)]
comp[, prop  := N / total]
fwrite(comp, file.path(DIR_RES, "fig6a_function_composition_135.csv"))

# Anchors: plasmid genes per genome bb 421 / bav 399 / afz 340 / gar 211
pg  <- unique(comp[location == "Plasmid", .(genome, sp, total)])
pgm <- pg[, round(mean(total)), by = sp]
stopifnot(pgm[sp == "bb", V1] == 421, pgm[sp == "bav", V1] == 399,
          pgm[sp == "afz", V1] == 340, pgm[sp == "gar", V1] == 211)
# Anchors: pooled antigen fraction (135) plasmid 6.43% vs chromosome 0.46%
# (exactly one row per genome within a location x category group, so a plain
#  sum(total) is correct here; the 5x-duplication issue only arises when
#  summing across categories)
pool <- comp[, .(N = sum(N), total = sum(total)), by = .(location, category)]
pool[, pct := round(100 * N / total, 2)]
stopifnot(abs(pool[location == "Plasmid"    & category == "Antigen", pct] - 6.43) < 0.05,
          abs(pool[location == "Chromosome" & category == "Antigen", pct] - 0.46) < 0.05)
# Anchors: gar plasmid Hypothetical 13.3%, Antigen 3.83 +/- 1.15 (CV ~30%; main
# text must carry a qualifier)
gar_hyp <- comp[sp == "gar" & location == "Plasmid" & category == "Hypothetical", mean(prop) * 100]
gar_ag  <- comp[sp == "gar" & location == "Plasmid" & category == "Antigen", prop * 100]
stopifnot(abs(gar_hyp - 13.3) < 0.2, abs(mean(gar_ag) - 3.83) < 0.05, abs(sd(gar_ag) - 1.15) < 0.05)
cat("sanity: composition anchors all passed\n")

## ============================================================
## Part 2: Fig 6a stacked bars + table-level Fisher cross-check
## ============================================================
# Cross-check: recompute the 129 Fisher from the base table (unique total first);
# must agree exactly with the gene-level authoritative version.
# (reproduction pitfall: total repeats across the 5 rows of a genome-location,
#  so a naive sum(total) inflates the denominator 5-fold)
d129  <- comp[sp %in% c("bb", "gar", "bav")]
tot129 <- unique(d129[, .(genome, location, total)])[, .(T = sum(total)), by = location]
chk_or <- function(cc) {
  pin  <- d129[location == "Plasmid"    & category == cc, sum(N)]
  cin  <- d129[location == "Chromosome" & category == cc, sum(N)]
  pout <- tot129[location == "Plasmid", T]    - pin
  cout <- tot129[location == "Chromosome", T] - cin
  unname(fisher.test(matrix(c(pin, pout, cin, cout), 2, byrow = TRUE))$estimate)
}
cats <- c("Antigen", "MGE/Plasmid", "Surface-associated", "Hypothetical", "Metabolic/Other")
for (cc in cats) {
  stopifnot(abs(chk_or(cc) - f129[category == cc, OR]) < 0.01)
}
cat("sanity: table-level Fisher identical to gene-level authoritative version\n")

mean_comp <- comp[, .(mean_prop = mean(prop)), by = .(sp, location, category)]
cat_levels <- c("Metabolic/Other", "Surface-associated", "Hypothetical",
                "MGE/Plasmid", "Antigen")
# Final palette (Antigen purple #7D3C98 on top of stack)
cat_cols <- c("Metabolic/Other"   = "#5DADE2",
              "Surface-associated" = "#1ABC9C",
              "Hypothetical"       = "#E74C3C",
              "MGE/Plasmid"        = "#F39C12",
              "Antigen"            = "#7D3C98")
mean_comp[, category := factor(category, levels = cat_levels)]
mean_comp[, location := factor(location, levels = c("Plasmid", "Chromosome"))]
mean_comp[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"),
                         labels = c("B. burgdorferi", "B. garinii",
                                    "B. afzelii†", "B. bavariensis"))]
p <- ggplot(mean_comp, aes(x = sp, y = mean_prop, fill = category)) +
  geom_bar(stat = "identity", width = 0.65, color = "white", linewidth = 0.4) +
  facet_wrap(~location, scales = "free_y") +
  scale_fill_manual(values = cat_cols, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean fraction of clustered genes") +
  theme_bw(base_size = 11) +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_text(face = "italic", angle = 35, hjust = 1),
        strip.background = element_rect(fill = "grey92", color = "black"),
        legend.position = "bottom")
print(p)
ggsave(file.path(DIR_RES, "fig6a_function_composition_stacked_bar.pdf"),
       p, width = 8, height = 4.5)
cat("=== 02 done: base table + Fisher + Fig 6a exported ===\n")
