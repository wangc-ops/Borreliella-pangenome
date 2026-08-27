# ============================================================
# Geographic region and host/vector composition of
# the four major pathogenic Borreliella species (n = 241),
# plus patchwork assembly of the full Fig. S1 (a/b/c).
#
# Input : Supplementary_Table_S1.xls (columns: species, geo_region,
#          host_group)
# Output: FigS1b_region.pdf, FigS1c_host.pdf, FigS1_combined.pdf
#
# Note: run after 11_Geographic_distribution_map.R so that p_a
#       exists in the environment (or source that script first).
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

s1 <- read_excel("Supplementary_Table_S1.xls")   # adjust path as needed

four <- c("Borreliella burgdorferi", "Borreliella garinii",
          "Borreliella afzelii", "Borreliella bavariensis")

## ---- Panel b: geographic region composition (geo_region, 100% coverage) ----
df_b <- s1 %>%
  filter(species %in% four) %>%
  mutate(species = factor(species, levels = four),
         geo_region = factor(geo_region,
                             levels = c("Europe", "North America",
                                        "East Asia", "Other_or_Unknown"))) %>%
  count(species, geo_region)

p_b <- ggplot(df_b, aes(species, n, fill = geo_region)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  scale_fill_manual(values = c("Europe"           = "#0072B2",
                               "North America"    = "#D55E00",
                               "East Asia"        = "#009E73",
                               "Other_or_Unknown" = "grey70"),
                    labels = c("Europe", "North America",
                               "East Asia", "Other/Unknown")) +
  scale_x_discrete(labels = c("B. burgdorferi", "B. garinii",
                              "B. afzelii", "B. bavariensis")) +
  labs(x = NULL, y = "Number of genomes", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "right")

## ---- Panel c: host/vector group composition (host_group, 100% coverage) ----
df_c <- s1 %>%
  filter(species %in% four) %>%
  mutate(species = factor(species, levels = four),
         host_group = factor(host_group,
                             levels = c("Tick", "Human", "Rodent",
                                        "Other_or_Unknown"))) %>%
  count(species, host_group)

p_c <- ggplot(df_c, aes(species, n, fill = host_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(n > 2, n, "")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  scale_fill_manual(values = c("Tick"             = "#E69F00",
                               "Human"            = "#56B4E9",
                               "Rodent"           = "#009E73",
                               "Other_or_Unknown" = "grey70"),
                    labels = c("Tick (vector)", "Human (patient)",
                               "Rodent", "Other/Unknown")) +
  scale_x_discrete(labels = c("B. burgdorferi", "B. garinii",
                              "B. afzelii", "B. bavariensis")) +
  labs(x = NULL, y = "Number of genomes", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "right")

print(p_b)
print(p_c)
ggsave("FigS1b_region.pdf", p_b, width = 6, height = 4.5)
ggsave("FigS1c_host.pdf",   p_c, width = 6, height = 4.5)

## ---- Full Fig. S1 assembly: map on top, panels b/c below ----
p_s1 <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "a")

print(p_s1)
ggsave("FigS1_combined.pdf", p_s1, width = 12, height = 9)
