# ============================================================
# 01. Intra-species core gene SNP distance
# Input : data/Supplementary_Table_S1.xls
#         data/snp_distance_pairwise_with_ids.csv
#         data/strain_plasmid_status_master_corrected.csv (B. bavariensis bimodality check)
# Output: results/intra_species_snp_distance.pdf
#         results/intra_species_snp_distance_summary.csv
# ============================================================
pkgs <- c("tidyverse", "readxl", "rstatix", "multcompView")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse); library(readxl); library(rstatix); library(multcompView)

dir.create("results", showWarnings = FALSE)

sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_colors <- c("Borreliella burgdorferi" = "#1F618D",
               "Borreliella garinii"     = "#E74C3C",
               "Borreliella afzelii"     = "#1ABC9C",
               "Borreliella bavariensis" = "#5DADE2")

# --- species annotation (S1 as reference) ---
s1 <- read_excel("data/Supplementary_Table_S1.xls") %>%
  mutate(Genome_ID = sub("_genomic$", "", Genome_ID),
         species = case_when(
           grepl("burgdorferi", species, ignore.case = TRUE) ~ "Borreliella burgdorferi",
           grepl("garinii",     species, ignore.case = TRUE) ~ "Borreliella garinii",
           grepl("afzelii",     species, ignore.case = TRUE) ~ "Borreliella afzelii",
           grepl("bavariensis", species, ignore.case = TRUE) ~ "Borreliella bavariensis",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(species)) %>% select(Genome_ID, species)

master <- read_csv("data/strain_plasmid_status_master_corrected.csv", show_col_types = FALSE) %>%
  mutate(Genome_ID = sub("_genomic$", "", Genome_ID)) %>%
  select(Genome_ID, has_plasmid_corrected)

# --- pairwise distances: rebuild comparison from S1 species at both ends ---
pairs <- read_csv("data/snp_distance_pairwise_with_ids.csv", show_col_types = FALSE) %>%
  mutate(genome_1 = sub("_genomic$", "", genome_1),
         genome_2 = sub("_genomic$", "", genome_2)) %>%
  left_join(s1 %>% rename(sp1 = species), by = c("genome_1" = "Genome_ID")) %>%
  left_join(s1 %>% rename(sp2 = species), by = c("genome_2" = "Genome_ID")) %>%
  filter(!is.na(sp1), !is.na(sp2)) %>%
  mutate(comparison = ifelse(sp1 == sp2, "intra", "inter"))
stopifnot(nrow(pairs) == 28920)   # C(241, 2)

# --- genome-level analysis unit: median intra-species distance per genome ---
intra <- pairs %>% filter(comparison == "intra")
genome_stat <- bind_rows(
  intra %>% transmute(Genome_ID = genome_1, species = sp1, d = snp_distance),
  intra %>% transmute(Genome_ID = genome_2, species = sp2, d = snp_distance)) %>%
  group_by(Genome_ID, species) %>%
  summarise(median_d = median(d), .groups = "drop")

dist_sum <- genome_stat %>%
  group_by(species) %>%
  summarise(n_genomes = n(),
            median_pct = median(median_d) * 100,
            q1_pct  = quantile(median_d, 0.25) * 100,
            q3_pct  = quantile(median_d, 0.75) * 100,
            min_pct = min(median_d) * 100,
            max_pct = max(median_d) * 100, .groups = "drop")
print(dist_sum)

# --- B. bavariensis bimodality check ---
bav_check <- genome_stat %>%
  filter(species == "Borreliella bavariensis") %>%
  left_join(master, by = "Genome_ID") %>%
  mutate(group = ifelse(median_d < 0.005, "clonal", "divergent"))
print(bav_check %>% count(group, has_plasmid_corrected))

# --- inter-species pairwise medians ---
inter_combo <- pairs %>%
  filter(comparison == "inter") %>%
  mutate(pair = paste(pmin(sp1, sp2), pmax(sp1, sp2), sep = " vs ")) %>%
  group_by(pair) %>%
  summarise(n_pairs = n(), median_pct = median(snp_distance) * 100,
            q1_pct = quantile(snp_distance, 0.25) * 100,
            q3_pct = quantile(snp_distance, 0.75) * 100, .groups = "drop") %>%
  arrange(median_pct)
print(inter_combo)

# --- Welch ANOVA + Games-Howell + CLD (4 species, genome level) ---
welch_res <- genome_stat %>% welch_anova_test(median_d ~ species)
gh <- genome_stat %>% games_howell_test(median_d ~ species)
print(welch_res); print(gh)

p_mat <- matrix(1, 4, 4, dimnames = list(sp_levels, sp_levels))
for (i in seq_len(nrow(gh))) {
  a <- as.character(gh$group1[i]); b <- as.character(gh$group2[i])
  p_mat[a, b] <- gh$p.adj[i]; p_mat[b, a] <- gh$p.adj[i]
}
raw_cld <- multcompLetters(p_mat, threshold = 0.05)$Letters
med_order <- dist_sum %>% arrange(desc(median_pct)) %>% pull(species) %>% as.character()
uniq <- unique(unlist(strsplit(raw_cld[med_order], "")))
remap <- setNames(letters[seq_along(uniq)], uniq)
cld_letters <- sapply(raw_cld, function(x)
  paste(sort(remap[strsplit(x, "")[[1]]]), collapse = ""))
print(cld_letters)

# --- summary table (intra + inter) ---
snp_summary <- bind_rows(
  dist_sum %>% mutate(section = "intra_species_genome_level", pair = NA_character_, n_pairs = NA_integer_),
  inter_combo %>% mutate(section = "inter_species_pairwise",
                         n_genomes = NA_integer_, min_pct = NA_real_, max_pct = NA_real_)) %>%
  relocate(section)
write_csv(snp_summary, "results/intra_species_snp_distance_summary.csv")

# --- plot ---
genome_stat <- genome_stat %>% mutate(species = factor(species, levels = sp_levels))
cld_df <- genome_stat %>%
  group_by(species) %>%
  summarise(ypos = max(median_d) * 1.06, .groups = "drop") %>%
  mutate(letter = cld_letters[as.character(species)])

fig <- ggplot(genome_stat, aes(x = species, y = median_d, fill = species)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.6,
               color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 1.8, shape = 21,
              color = "black", stroke = 0.3, alpha = 0.7) +
  geom_text(data = cld_df, aes(x = species, y = ypos, label = letter),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = sp_colors) +
  scale_x_discrete(labels = c("Borreliella burgdorferi" = "B. burgdorferi",
                              "Borreliella garinii"     = "B. garinii",
                              "Borreliella afzelii"     = "B. afzelii",
                              "Borreliella bavariensis" = "B. bavariensis")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.5),
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Species", y = "Median intra-species core SNP distance") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        axis.ticks.length = unit(2.5, "pt"),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "italic", size = 11, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))
print(fig)
ggsave("results/intra_species_snp_distance.pdf", fig, width = 7, height = 5.5, dpi = 300)
cat("Done: intra-species core SNP distance\n")