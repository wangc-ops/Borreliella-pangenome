# plasmid_cluster_proxy_validation.R
# Validates de novo cluster counts as a conservative proxy of plasmidome
# structural complexity: (a) true replicon counts in Complete Genomes;
# (b) PFam32 partition loci across all 135 plasmid-carrying genomes.
# Inputs under ./data: Supplementary_Table_S1.xls,
#   strain_plasmid_status_master_corrected.csv, plasmid_clusters_cluster.tsv,
#   gff3/ (Bakta *_genomic.gff3)
library(data.table); library(readxl); library(stringr)
library(ggplot2); library(patchwork)

to_sp_code <- function(x) {
  out <- dplyr::case_when(
    x == "Borreliella burgdorferi" ~ "bb", x == "Borreliella garinii" ~ "gar",
    x == "Borreliella afzelii" ~ "afz", x == "Borreliella bavariensis" ~ "bav",
    TRUE ~ NA_character_)
  if (any(is.na(out))) stop("unmapped species")
  out
}

master <- fread("data/strain_plasmid_status_master_corrected.csv")
s1 <- as.data.table(read_excel("data/Supplementary_Table_S1.xls", sheet = 1))
s1[, gcf := str_extract(assembly, "GCF_\\d+\\.\\d+")]
sub <- merge(master[, .(gcf = str_extract(Genome_ID, "GCF_\\d+\\.\\d+"),
                        has_plasmid_corrected)],
             s1[, .(gcf, species, assembly_level)], by = "gcf")[has_plasmid_corrected == TRUE]
stopifnot(nrow(sub) == 135)
sub[, sp := to_sp_code(species)]

# Replicon counts (>5 kb non-chromosomal contigs) and PFam32 partition loci
files <- list.files("data/gff3", pattern = "_genomic\\.gff3$", full.names = TRUE)
names(files) <- str_extract(basename(files), "GCF_\\d+\\.\\d+")
files <- files[sub$gcf]; stopifnot(!any(is.na(files)))
parse_one <- function(f) {
  L <- readLines(f, warn = FALSE)
  sr <- L[startsWith(L, "##sequence-region")]
  lens <- as.numeric(str_extract(sr, "\\d+$"))
  attr <- str_split_fixed(L[!startsWith(L, "#")], "\t", 9)[, 9]
  list(n_plasmid_contig = sum(lens > 5000) - 1L,
       n_partition = sum(str_detect(attr, fixed(
         "product=chromosome replication/partitioning protein"))))
}
sub <- cbind(sub, rbindlist(lapply(files, parse_one)))
# Anchor: B31 must yield 21 replicons
stopifnot(sub[gcf == "GCF_000008685.2"]$n_plasmid_contig == 21)

# Cluster counts (cluster representatives <= 500 kb)
cl <- fread("data/plasmid_clusters_cluster.tsv", sep = "\t", header = FALSE,
            col.names = c("rep", "member"))
cl[, `:=`(gcf = str_extract(member, "GCF_\\d+\\.\\d+"),
          rep_len = as.numeric(str_extract(rep, "(?<=length=)\\d+")))]
nc <- cl[gcf %in% sub$gcf & rep_len <= 5e5, .(n_cluster = uniqueN(rep)), by = gcf]
m <- merge(sub, nc, by = "gcf", all.x = TRUE)
m[is.na(n_cluster), n_cluster := 0]
stopifnot(nrow(m) == 135, all(m$n_cluster > 0))
fwrite(m[, .(gcf, sp, assembly_level, n_plasmid_contig, n_partition, n_cluster)],
       "plasmid_cluster_proxy_validation_source.tsv")

# Statistics
cg <- m[assembly_level == "Complete Genome"]
print(cor.test(cg$n_plasmid_contig, cg$n_cluster, method = "spearman"))  # rho ~ 0.987
print(cor.test(m$n_partition, m$n_cluster, method = "spearman"))        # rho ~ 0.895
print(m[, .(n = .N, rho = round(cor(n_partition, n_cluster, method = "spearman"), 3)),
        by = assembly_level])

# Figure
sp_cols <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
sp_labs <- c(bb = "B. burgdorferi", gar = "B. garinii",
             afz = "B. afzelii", bav = "B. bavariensis")
m[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"))]
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))
mk <- function(dat, xvar, xlab, rho_lab) {
  ggplot(dat, aes(x = .data[[xvar]], y = n_cluster, fill = sp)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey60", linewidth = 0.5) +
    geom_point(shape = 21, size = 1.8, color = "black", stroke = 0.3,
               alpha = 0.85,
               position = position_jitter(width = 0.15, height = 0.15, seed = 42)) +
    scale_fill_manual(values = sp_cols, labels = sp_labs) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6,
             label = rho_lab, parse = TRUE, size = 3.4) +
    labs(x = xlab, y = "Plasmid clusters per genome") + theme_pub
}
pa <- mk(cg, "n_plasmid_contig",
         "Plasmid replicons per genome (Complete Genomes, n = 47)",
         '"Spearman"~rho==0.987*","~italic(P)*" < 0.001"')
pb <- mk(m, "n_partition", "PFam32 partition loci per genome (n = 135)",
         '"Spearman"~rho==0.895*","~italic(P)*" < 0.001"')
fig <- (pa | pb) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave("plasmid_cluster_proxy_validation.pdf", fig, width = 10, height = 4.8)