# ============================================================
# 07. Cluster size-frequency trade-off (supplementary)
#     Spearman correlation between cluster representative length
#     and detection rate, at two scopes
# Input : data/Supplementary_Table_S1.xls
#         data/strain_plasmid_status_master_corrected.csv
#         data/plasmid_clusters_cluster_nochrom.tsv (483 plasmid clusters, chromosomes removed)
# Output: results/cluster_size_frequency_tradeoff.pdf
#         results/cluster_size_frequency.csv
# ============================================================
pkgs <- c("tidyverse", "readxl")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse); library(readxl)

dir.create("results", showWarnings = FALSE)

master <- read_csv("data/strain_plasmid_status_master_corrected.csv", show_col_types = FALSE) %>%
  mutate(Genome_ID = sub("_genomic$", "", Genome_ID))
plasmid_ids <- master %>% filter(has_plasmid_corrected == TRUE) %>% pull(Genome_ID)

s1 <- read_excel("data/Supplementary_Table_S1.xls") %>%
  mutate(Genome_ID = sub("_genomic$", "", Genome_ID),
         species = case_when(
           grepl("burgdorferi", species, ignore.case = TRUE) ~ "Borreliella burgdorferi",
           grepl("garinii",     species, ignore.case = TRUE) ~ "Borreliella garinii",
           grepl("afzelii",     species, ignore.case = TRUE) ~ "Borreliella afzelii",
           grepl("bavariensis", species, ignore.case = TRUE) ~ "Borreliella bavariensis",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(species)) %>% select(Genome_ID, species)

# --- cluster rep length is embedded in rep_id ("length=xxxxx") ---
clusters <- read_tsv("data/plasmid_clusters_cluster_nochrom.tsv",
                     col_names = c("rep_id", "member_id"), show_col_types = FALSE) %>%
  mutate(Genome_ID = sub("_genomic$", "", str_extract(member_id, "^[^|]+")),
         length_bp = as.numeric(str_extract(rep_id, "(?<=length=)\\d+")))
stopifnot(!any(is.na(clusters$length_bp)))
cat("Cluster rep length range:", min(clusters$length_bp), "-",
    max(clusters$length_bp), "bp (no chromosome-scale clusters)\n")

pav_len <- clusters %>%
  distinct(Genome_ID, rep_id, length_bp) %>%
  inner_join(s1, by = "Genome_ID") %>%
  filter(Genome_ID %in% plasmid_ids)

# --- scope 1: all 483 clusters, 135-genome denominator ---
scope_all <- pav_len %>%
  group_by(rep_id) %>%
  summarise(length_bp = first(length_bp),
            rate = n_distinct(Genome_ID) / 135 * 100, .groups = "drop") %>%
  mutate(scope = "All clusters (n=483, 135 genomes)")

# --- scope 2: clusters detected in B. burgdorferi, 91-genome denominator ---
bb_genomes <- pav_len %>% filter(species == "Borreliella burgdorferi") %>%
  distinct(Genome_ID) %>% pull()
stopifnot(length(bb_genomes) == 91)
scope_bb <- pav_len %>%
  filter(species == "Borreliella burgdorferi") %>%
  group_by(rep_id) %>%
  summarise(length_bp = first(length_bp),
            rate = n_distinct(Genome_ID) / 91 * 100, .groups = "drop") %>%
  mutate(scope = "B. burgdorferi clusters (91 genomes)")

tradeoff <- bind_rows(scope_all, scope_bb)

# --- Spearman correlations (both scopes) ---
cat("Spearman correlations (length vs detection rate):\n")
for (sc in unique(tradeoff$scope)) {
  dd <- tradeoff %>% filter(scope == sc)
  ct <- cor.test(dd$length_bp, dd$rate, method = "spearman", exact = FALSE)
  cat(sc, ": rho =", round(unname(ct$estimate), 3),
      " P =", format.pval(ct$p.value, digits = 3), " n =", nrow(dd), "\n")
}

write_csv(tradeoff %>% arrange(scope, desc(rate)), "results/cluster_size_frequency.csv")

# --- plot (two panels) ---
fig <- ggplot(tradeoff, aes(x = rate, y = length_bp / 1000)) +
  geom_point(size = 1.6, alpha = 0.55, color = "#1F618D", shape = 16) +
  geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", linewidth = 0.8) +
  facet_wrap(~scope, ncol = 2) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Detection rate (%)", y = "Cluster representative length (kb, log10)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        axis.ticks.length = unit(2.5, "pt"),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))
print(fig)
ggsave("results/cluster_size_frequency_tradeoff.pdf", fig, width = 10, height = 4.8, dpi = 300)
cat("Done: cluster size-frequency trade-off\n")