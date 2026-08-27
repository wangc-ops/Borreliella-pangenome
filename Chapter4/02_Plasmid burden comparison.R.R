# ============================================================
# 02. Plasmid burden comparison (cluster count per genome)
# Input : data/Supplementary_Table_S1.xls
#         data/strain_plasmid_status_master_corrected.csv
#         data/plasmid_clusters_cluster_nochrom.tsv (483 plasmid clusters, chromosomes removed)
# Output: results/cluster_burden.pdf
#         results/plasmid_burden_summary.csv
# Note  : analyses conditional on 135 plasmid-carrying genomes;
#         B. afzelii (n=6) shown descriptively only (dagger), excluded from tests
# ============================================================
pkgs <- c("tidyverse", "readxl", "rstatix", "multcompView")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse); library(readxl); library(rstatix); library(multcompView)

dir.create("results", showWarnings = FALSE)

sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_test3  <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella bavariensis")
sp_colors <- c("Borreliella burgdorferi" = "#1F618D",
               "Borreliella garinii"     = "#E74C3C",
               "Borreliella afzelii"     = "#1ABC9C",
               "Borreliella bavariensis" = "#5DADE2")

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

clusters <- read_tsv("data/plasmid_clusters_cluster_nochrom.tsv",
                     col_names = c("rep_id", "member_id"), show_col_types = FALSE) %>%
  mutate(Genome_ID = sub("_genomic$", "", str_extract(member_id, "^[^|]+")))

# --- 135-genome x cluster PAV matrix ---
pav_sp <- clusters %>%
  distinct(Genome_ID, rep_id) %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = Genome_ID, names_from = rep_id,
              values_from = present, values_fill = 0) %>%
  inner_join(s1, by = "Genome_ID") %>%
  filter(Genome_ID %in% plasmid_ids)
stopifnot(nrow(pav_sp) == 135, ncol(pav_sp) - 2 == 483)

# --- burden = cluster count per genome ---
burden_df <- pav_sp %>%
  transmute(Genome_ID, species,
            burden = rowSums(across(-c(Genome_ID, species)))) %>%
  mutate(species = factor(species, levels = sp_levels))

burden_sum <- burden_df %>%
  group_by(species) %>%
  summarise(n = n(), median = median(burden),
            q1 = quantile(burden, 0.25), q3 = quantile(burden, 0.75), .groups = "drop")
print(burden_sum)

# --- Welch ANOVA + Games-Howell + CLD (bb/gar/bav only) ---
d <- burden_df %>% filter(species %in% sp_test3) %>%
  mutate(species = factor(as.character(species), levels = sp_test3))
welch_res <- welch_anova_test(d, burden ~ species)
gh <- games_howell_test(d, burden ~ species)
print(welch_res); print(gh)

p_mat <- matrix(1, 3, 3, dimnames = list(sp_test3, sp_test3))
for (i in seq_len(nrow(gh))) {
  a <- as.character(gh$group1[i]); b <- as.character(gh$group2[i])
  p_mat[a, b] <- gh$p.adj[i]; p_mat[b, a] <- gh$p.adj[i]
}
raw_cld <- multcompLetters(p_mat, threshold = 0.05)$Letters
med_order <- d %>% group_by(species) %>% summarise(m = median(burden), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(species) %>% as.character()
uniq <- unique(unlist(strsplit(raw_cld[med_order], "")))
remap <- setNames(letters[seq_along(uniq)], uniq)
cld_letters <- sapply(raw_cld, function(x)
  paste(sort(remap[strsplit(x, "")[[1]]]), collapse = ""))
print(cld_letters)

write_csv(burden_sum %>% mutate(cld = cld_letters[as.character(species)]),
          "results/plasmid_burden_summary.csv")

# --- plot ---
cld_df <- burden_df %>%
  filter(as.character(species) %in% names(cld_letters)) %>%
  group_by(species) %>%
  summarise(ypos = max(burden) + 0.06 * diff(range(burden_df$burden)), .groups = "drop") %>%
  mutate(letter = cld_letters[as.character(species)])

fig <- ggplot(burden_df, aes(x = species, y = burden, fill = species)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.6,
               color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 1.8, shape = 21,
              color = "black", stroke = 0.3, alpha = 0.7) +
  geom_text(data = cld_df, aes(x = species, y = ypos, label = letter),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = sp_colors) +
  scale_x_discrete(labels = c("Borreliella burgdorferi" = "B. burgdorferi",
                              "Borreliella garinii"     = "B. garinii",
                              "Borreliella afzelii"     = "B. afzelii †",
                              "Borreliella bavariensis" = "B. bavariensis")) +
  labs(x = "Species", y = "Plasmid cluster count per genome") +
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
ggsave("results/cluster_burden.pdf", fig, width = 7, height = 5.5, dpi = 300)
cat("Done: plasmid burden comparison\n")