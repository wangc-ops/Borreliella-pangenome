# ============================================================
# 04. Antigen coupling index, between-species comparison (Fig 6b)
#   Part 1: rebuild coupling base table (established builder logic)
#   Part 2: Welch ANOVA + Games-Howell + automatic CLD + boxplot
# Input : data/antigen13_genes.csv                  (from script 01)
#         data/gff_contigs.rds                      (from script 01)
#         data/fig6a_function_composition_135.csv   (from script 02; Hypothetical rows only)
#         data/strain_plasmid_status_master_corrected.csv
# Output: results/coupling_13family_135.csv  (genome, ag_total, ag_plasmid, antigen_coupling,
#                                             hyp_chrom, hyp_plasmid, hyp_coupling, sp; 135 rows)
#         results/coupling_descriptives.csv
#         results/fig6b_antigen_coupling_boxplot.pdf
# Note  : metric definitions (M&M) --
#         antigen_coupling = plasmid-borne antigen genes / total antigen genes
#                            (gene level, 13 families)
#         hyp_coupling     = plasmid Hypothetical genes / total Hypothetical genes
#                            (control metric)
#         all 135 plasmid carriers plotted; tests restricted to bb/gar/bav
#         (n=129); afz (n=6) shown with dagger only
#         CLD: multcompLetters, remapped by descending mean (a = highest mean);
#         the old hardcoded a/b/b/c is abolished
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(rstatix); library(multcompView)
})

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

## ============================================================
## Part 1: coupling base-table rebuild
## ============================================================
mast <- fread(file.path(DIR_DATA, "strain_plasmid_status_master_corrected.csv"))
mast[, genome := sub("_genomic$", "", Genome_ID)]
sp_map <- c("Borreliella burgdorferi" = "bb", "Borreliella garinii" = "gar",
            "Borreliella afzelii" = "afz", "Borreliella bavariensis" = "bav")
mast[, sp := sp_map[species]]
stopifnot(!any(is.na(mast$sp)))
carriers <- mast[as.logical(has_plasmid_corrected) == TRUE]
stopifnot(nrow(carriers) == 135)

# Antigen side: antigen13 gene-level table x replicon rule (longest contig = chromosome)
a13 <- fread(file.path(DIR_DATA, "antigen13_genes.csv"))
ct  <- readRDS(file.path(DIR_DATA, "gff_contigs.rds"))
ct[, genome := sub("_genomic$", "", genome)]
chrom <- ct[, .SD[which.max(len)], by = genome][, .(genome, chrom_contig = seqid)]
a13 <- merge(a13, chrom, by = "genome")
a13[, on_plasmid := contig != chrom_contig]
agg <- a13[, .(ag_total = .N, ag_plasmid = sum(on_plasmid)), by = genome]
c135 <- merge(carriers[, .(genome, sp)], agg, by = "genome")
stopifnot(nrow(c135) == 135, all(c135$ag_total > 0))
c135[, antigen_coupling := ag_plasmid / ag_total]

# Hypothetical side: from the script-02 composition table (Panaroo-clustered genes)
comp <- fread(file.path(DIR_DATA, "fig6a_function_composition_135.csv"))
h <- dcast(comp[category == "Hypothetical"], genome ~ location,
           value.var = "N", fill = 0)
# Fix: genomes with zero plasmid Hypothetical genes may lose the column/NA
# after dcast -> fill 0
if (!"Plasmid" %in% names(h)) h[, Plasmid := 0L]
h[is.na(Plasmid), Plasmid := 0L]
h[is.na(Chromosome), Chromosome := 0L]
# Verified: GCF_038801875.1 is a genuine minimal plasmidome (26 clustered
# plasmid genes, 1 antigen), not a data error
c135 <- merge(c135,
              h[, .(genome, hyp_chrom = Chromosome, hyp_plasmid = Plasmid)],
              by = "genome")
c135[, hyp_coupling := hyp_plasmid / (hyp_chrom + hyp_plasmid)]
stopifnot(nrow(c135) == 135, !any(is.na(c135$hyp_coupling)))
setcolorder(c135, c("genome", "ag_total", "ag_plasmid", "antigen_coupling",
                    "hyp_chrom", "hyp_plasmid", "hyp_coupling", "sp"))
setkey(c135, genome)
fwrite(c135, file.path(DIR_RES, "coupling_13family_135.csv"))

## ============================================================
## Part 2: statistics + Fig 6b
## ============================================================
cpl <- copy(c135)
cpl[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"))]

# Tests use bb/gar/bav only (n=129); afz is plotted with a dagger only
d <- cpl[sp %in% c("bb", "gar", "bav")]
d[, sp := droplevels(sp)]
stopifnot(nrow(d) == 129)

# ---------- descriptive statistics ----------
desc_ant <- d[, .(mean = mean(antigen_coupling), sd = sd(antigen_coupling)), by = sp]
desc_hyp <- d[, .(mean = mean(hyp_coupling),     sd = sd(hyp_coupling)),     by = sp]
desc <- merge(desc_ant, desc_hyp, by = "sp", suffixes = c("_antigen", "_hyp"))
print(desc)
fwrite(desc, file.path(DIR_RES, "coupling_descriptives.csv"))
# afz display values (not in tests): antigen 0.777 +/- 0.125; hyp 0.309 +/- 0.114
print(cpl[sp == "afz", .(mean_antigen = mean(antigen_coupling),
                         sd_antigen   = sd(antigen_coupling),
                         mean_hyp     = mean(hyp_coupling),
                         sd_hyp       = sd(hyp_coupling))])

# ---------- Welch ANOVA + Games-Howell ----------
wa_ant <- welch_anova_test(d, antigen_coupling ~ sp)
gh_ant <- as.data.table(games_howell_test(d, antigen_coupling ~ sp))
wa_hyp <- welch_anova_test(d, hyp_coupling ~ sp)
gh_hyp <- as.data.table(games_howell_test(d, hyp_coupling ~ sp))
print(wa_ant); print(gh_ant)
print(wa_hyp); print(gh_hyp)

# ---------- hard assertions (final values locked) ----------
stopifnot(abs(wa_ant$F - 47.51) < 0.05, wa_ant$p < 1e-9, all(gh_ant$p.adj < 0.001))
stopifnot(abs(desc[sp == "bb",  mean_antigen] - 0.887) < 0.005,
          abs(desc[sp == "gar", mean_antigen] - 0.635) < 0.005,
          abs(desc[sp == "bav", mean_antigen] - 0.795) < 0.005)
stopifnot(abs(wa_hyp$F - 11.54) < 0.05, wa_hyp$p < 1e-3)
stopifnot(abs(desc[sp == "bb",  mean_hyp] - 0.482) < 0.005,
          abs(desc[sp == "gar", mean_hyp] - 0.378) < 0.005,
          abs(desc[sp == "bav", mean_hyp] - 0.508) < 0.005)
stopifnot(gh_hyp[group1 == "bb" & group2 == "bav", p.adj] > 0.05)  # bb-bav ns (0.469)
cat("sanity: Welch + Games-Howell hard assertions all passed\n")

# ---------- CLD (automatic; highest-mean group remapped to "a") ----------
make_cld <- function(gh_dt, means_vec) {
  pv <- gh_dt$p.adj
  names(pv) <- paste(gh_dt$group1, gh_dt$group2, sep = "-")
  raw <- multcompLetters(pv)$Letters
  o <- names(means_vec)[order(-means_vec)]
  uniq <- unique(raw[o])
  mp <- setNames(letters[seq_along(uniq)], uniq)
  data.table(sp = names(raw), cld = unname(mp[raw]))
}
mv <- setNames(desc_ant$mean, as.character(desc_ant$sp))
cld_ant <- make_cld(gh_ant, mv)   # expected: bb=a, bav=b, gar=c
cld_ant[, y := sapply(sp, function(s)
  max(d[sp == s, antigen_coupling], na.rm = TRUE)) + 0.07]
print(cld_ant)

# ---------- Fig 6b ----------
sp_col <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
sp_labels <- c(bb  = "italic(B.~burgdorferi)",
               gar = "italic(B.~garinii)",
               afz = "italic(B.~afzelii)*'†'",
               bav = "italic(B.~bavariensis)")
p <- ggplot(cpl, aes(x = sp, y = antigen_coupling, fill = sp)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.6) +
  geom_jitter(shape = 21, width = 0.12, size = 1.8,
              color = "black", stroke = 0.3) +
  geom_text(data = cld_ant, aes(x = sp, y = y, label = cld),
            inherit.aes = FALSE, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = sp_col, guide = "none") +
  scale_x_discrete(labels = function(x) parse(text = sp_labels[as.character(x)])) +
  scale_y_continuous(limits = c(0, 1.18),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Species", y = "Antigen coupling index") +
  theme_classic(base_size = 12) +
  theme(axis.line   = element_line(color = "black"),
        axis.text.x = element_text(size = 8.5, angle = 30, hjust = 1, vjust = 1),
        plot.margin = margin(t = 5, r = 10, b = 15, l = 5))
print(p)
ggsave(file.path(DIR_RES, "fig6b_antigen_coupling_boxplot.pdf"),
       p, width = 7.0, height = 4.5)
cat("=== 04 done: coupling base table + Fig 6b exported ===\n")
