# ============================================================
# 05. Plasmid cluster richness x antigen family richness dose response (Fig 6c)
# Input : data/strain_plasmid_status_master_corrected.csv
#         data/antigen13_genes.csv              (from script 01)
#         data/plasmid_clusters_cluster.tsv     (plasmid cluster table)
#         data/gff_contigs.rds                  (from script 01; span/contigs as candidate x)
# Output: results/fig6c_spearman_results.csv
#         results/fig6c_plasmid_antigen_dose_scatter.pdf
# Note  : question -- does a larger plasmidome (more clusters) carry a broader
#         antigen repertoire?
#         pre-registered criterion: panel kept only if the within-bb Spearman
#         correlation is significant (passed; values locked from actual output)
#         sensitivity: recomputed after dropping the minimal plasmidome
#         GCF_038801875.1 (n=90; declared in M&M)
#         all 135 carriers plotted; Spearman restricted to bb/gar/bav;
#         afz not tested
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

# ---------- species mapping ----------
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
carriers <- master[has_plasmid_corrected == TRUE]
stopifnot(nrow(carriers) == 135)

# ---------- y: antigen family richness (how many of the 13 families detected) ----------
ag <- fread(file.path(DIR_DATA, "antigen13_genes.csv"))
stopifnot(uniqueN(ag$fam13) == 13)
rich_ag <- unique(ag[, .(genome, fam13)])[, .(antigen_richness = uniqueN(fam13)), by = genome]

# ---------- candidate x: plasmid cluster count / plasmid span / plasmid contigs ----------
clu <- fread(file.path(DIR_DATA, "plasmid_clusters_cluster.tsv"), header = FALSE, sep = "\t")
clu[, genome_id := sub("\\|.*$", "", V2)]
stopifnot(mean(clu$genome_id %in% master$Genome_ID) > 0.95)
rich_clu <- unique(clu[, .(cluster = V1, genome_id)])[
  , .(plasmid_clusters = uniqueN(cluster)), by = genome_id]

gc <- readRDS(file.path(DIR_DATA, "gff_contigs.rds"))
rich_span <- gc[replicon == "plasmid",
                .(plasmid_span_bp = sum(len), plasmid_contigs = .N), by = genome]

# ---------- analysis table (135 carriers) ----------
dt <- carriers[, .(genome_id = Genome_ID, sp)]
dt <- merge(dt, rich_ag,   by.x = "genome_id", by.y = "genome", all.x = TRUE)
dt <- merge(dt, rich_clu,  by = "genome_id", all.x = TRUE)
dt <- merge(dt, rich_span, by.x = "genome_id", by.y = "genome", all.x = TRUE)
dt[is.na(antigen_richness), antigen_richness := 0L]
stopifnot(nrow(dt) == 135, !any(is.na(dt$plasmid_clusters)))

# ---------- Spearman: bb/gar/bav x three candidate x ----------
run_spear <- function(d, xvar, label) {
  ct_ <- suppressWarnings(cor.test(d[[xvar]], d$antigen_richness,
                                   method = "spearman", exact = FALSE))
  data.table(group = label, x_var = xvar, n = nrow(d),
             rho = round(unname(ct_$estimate), 3), P = signif(ct_$p.value, 3))
}
res <- rbindlist(lapply(c("bb", "gar", "bav"), function(s) {
  d <- dt[sp == s]
  rbindlist(lapply(c("plasmid_clusters", "plasmid_span_bp", "plasmid_contigs"),
                   function(v) run_spear(d, v, s)))
}))
print(res)
fwrite(res, file.path(DIR_RES, "fig6c_spearman_results.csv"))

# Pre-registered criterion + final value locks (main-text run values):
#   bb x plasmid_clusters: rho = 0.505, P = 3.4e-7 (n=91); gar: rho ~ 0, P > 0.05
verdict <- res[group == "bb" & x_var == "plasmid_clusters"]
stopifnot(abs(verdict$rho - 0.505) < 0.001, verdict$P < 1e-6, verdict$n == 91)
stopifnot(res[group == "gar" & x_var == "plasmid_clusters", P] > 0.05)

# ---------- sensitivity: drop minimal plasmidome GCF_038801875.1
#   (main text: rho = 0.487, P < 0.001) ----------
bb_no_out <- dt[sp == "bb" & !grepl("GCF_038801875.1", genome_id)]
ct_sens <- suppressWarnings(cor.test(bb_no_out$plasmid_clusters, bb_no_out$antigen_richness,
                                     method = "spearman", exact = FALSE))
cat("Sensitivity (n =", nrow(bb_no_out), "): rho =", round(unname(ct_sens$estimate), 3),
    ", P =", signif(ct_sens$p.value, 3), "\n")
stopifnot(nrow(bb_no_out) == 90,
          abs(round(unname(ct_sens$estimate), 3) - 0.487) < 0.001,
          ct_sens$p.value < 0.001)

# ---------- Fig 6c ----------
sp_cols <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
dt[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"))]
ann <- paste0("italic(B.~burgdorferi)*': Spearman '*rho == '",
              sprintf("%.2f", verdict$rho), "'*', '*italic(P) < 0.001'")

p <- ggplot(dt, aes(x = plasmid_clusters, y = antigen_richness)) +
  geom_jitter(aes(fill = sp), shape = 21, size = 1.8,
              width = 0.12, height = 0.12, alpha = 0.85, color = "black", stroke = 0.3) +
  geom_smooth(data = dt[sp == "bb"], method = "lm", se = TRUE,
              color = "black", linewidth = 0.6, fill = "grey70", alpha = 0.3) +
  annotate("text", x = -Inf, y = Inf, label = ann, parse = TRUE,
           hjust = -0.05, vjust = 1.5, size = 3.5) +
  scale_fill_manual(values = sp_cols,
                    labels = c("B. burgdorferi", "B. garinii",
                               "B. afzelii †", "B. bavariensis"), name = NULL) +
  scale_y_continuous(breaks = 0:13) +
  labs(x = "Plasmid cluster richness", y = "Antigen family richness (of 13)") +
  theme_classic(base_size = 11) +
  theme(legend.text = element_text(face = "italic"),
        axis.line  = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
print(p)
ggsave(file.path(DIR_RES, "fig6c_plasmid_antigen_dose_scatter.pdf"),
       p, width = 6.5, height = 4.5)
cat("=== 05 done: Spearman table + Fig 6c exported ===\n")
