# ============================================================
# 05. cp32-like plasmid family burden comparison
# Input : data/Supplementary_Table_S1.xls
#         data/strain_plasmid_status_master_corrected.csv
#         data/cp32_like_burden.csv (genome_id / cp32_count; explicit zeros, 241 rows)
# Output: results/cp32_burden.pdf
#         results/cp32_burden_summary.csv
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

# --- species from S1 (drop any species column in the burden file) ---
cp32_df <- read_csv("data/cp32_like_burden.csv", show_col_types = FALSE) %>%
  mutate(Genome_ID = sub("_genomic$", "", genome_id)) %>%
  select(-any_of("species")) %>%
  inner_join(s1, by = "Genome_ID") %>%
  filter(Genome_ID %in% plasmid_ids) %>%
  transmute(Genome_ID,
            species = factor(.data$species, levels = sp_levels),
            cp32_count)
stopifnot(nrow(cp32_df) == 135)

cp32_sum <- cp32_df %>%
  group_by(species) %>%
  summarise(n = n(), median = median(cp32_count),
            q1 = quantile(cp32_count, 0.25), q3 = quantile(cp32_count, 0.75), .groups = "drop")
print(cp32_sum)

# --- Welch ANOVA + Games-Howell + CLD (bb/gar/bav only) ---
d <- cp32_df %>% filter(species %in% sp_test3) %>%
  mutate(species = factor(as.character(species), levels = sp_test3))
welch_res <- welch_anova_test(d, cp32_count ~ species)
gh <- games_howell_test(d, cp32_count ~ species)
print(welch_res); print(gh)

p_mat <- matrix(1, 3, 3, dimnames = list(sp_test3, sp_test3))
for (i in seq_len(nrow(gh))) {
  a <- as.character(gh$group1[i]); b <- as.character(gh$group2[i])
  p_mat[a, b] <- gh$p.adj[i]; p_mat[b, a] <- gh$p.adj[i]
}
raw_cld <- multcompLetters(p_mat, threshold = 0.05)$Letters
med_order <- d %>% group_by(species) %>% summarise(m = median(cp32_count), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(species) %>% as.character()
uniq <- unique(unlist(strsplit(raw_cld[med_order], "")))
remap <- setNames(letters[seq_along(uniq)], uniq)
cld_letters <- sapply(raw_cld, function(x)
  paste(sort(remap[strsplit(x, "")[[1]]]), collapse = ""))
print(cld_letters)

write_csv(cp32_sum %>% mutate(cld = cld_letters[as.character(species)]),
          "results/cp32_burden_summary.csv")

# --- plot ---
cld_df <- cp32_df %>%
  filter(as.character(species) %in% names(cld_letters)) %>%
  group_by(species) %>%
  summarise(ypos = max(cp32_count) + 0.06 * diff(range(cp32_df$cp32_count)), .groups = "drop") %>%
  mutate(letter = cld_letters[as.character(species)])

fig <- ggplot(cp32_df, aes(x = species, y = cp32_count, fill = species)) +
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
  labs(x = "Species", y = "cp32-like plasmid copies per genome") +
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
ggsave("results/cp32_burden.pdf", fig, width = 7, height = 5.5, dpi = 300)
cat("Done: cp32-like plasmid family burden comparison\n")