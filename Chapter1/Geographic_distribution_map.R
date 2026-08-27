# ============================================================
# Geographic distribution of the 287 Borreliella genomes
#
# Input : Supplementary_Table_S1.xls (columns: species, geo_loc)
# Output: FigS1a_map.pdf
#
# Notes:
#   - 228 of 287 genomes carry valid country information (geo_loc);
#     the remaining 59 (placeholders such as "missing",
#     "not collected", "Unknown") are not shown (stated in caption).
#   - Precise coordinates (lat_lon) are available for only 68 genomes
#     and are strongly biased (only 3 for B. burgdorferi), so a
#     country-centroid bubble map is used instead; species sharing a
#     country are slightly jittered to avoid overplotting.
#   - set.seed(1) makes the jitter reproducible.
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

s1 <- read_excel("Supplementary_Table_S1.xls")   # adjust path as needed

# Four major pathogenic species analyzed in the main text
four <- c("Borreliella burgdorferi", "Borreliella garinii",
          "Borreliella afzelii", "Borreliella bavariensis")

# Country centroids (display-level approximations for all 17 countries
# present in the dataset)
centroids <- tribble(
  ~country,      ~lat,  ~lon,
  "Germany",     51.1,  10.4,  "Japan",       36.2, 138.3,
  "Canada",      56.1, -106.3, "USA",         39.8, -98.6,
  "Russia",      55.0,  90.0,  "Norway",      64.5,  17.0,
  "Portugal",    39.6,  -8.0,  "France",      46.6,   2.4,
  "Slovenia",    46.1,  14.9,  "Spain",       40.3,  -3.7,
  "South Korea", 36.4, 127.9,  "China",       35.0, 104.0,
  "Austria",     47.6,  14.1,  "Netherlands", 52.2,   5.3,
  "Chile",      -31.8, -71.0,  "Denmark",     56.0,  10.0,
  "Belarus",     53.5,  28.0
)

# Missing-value placeholders in geo_loc
bad <- c("missing", "not applicable", "Not Applicable", "Unknown",
         "not recorded", "not collected")

# Extract country name (geo_loc format: "Japan: island of Honshu ...")
df <- s1 %>%
  filter(!is.na(geo_loc), !geo_loc %in% bad) %>%
  mutate(country = str_trim(str_split_fixed(geo_loc, ":", 2)[, 1])) %>%
  inner_join(centroids, by = "country")

nrow(df)   # 228 expected

# Colors: four major species match the main-text palette;
# remaining species use a muted pastel set
pal4 <- c("Borreliella burgdorferi" = "#1F618D",
          "Borreliella garinii"     = "#E74C3C",
          "Borreliella afzelii"     = "#1ABC9C",
          "Borreliella bavariensis" = "#5DADE2")

other_sp <- sort(setdiff(unique(df$species), four))
pastel <- c("#B39DDB", "#FFCC80", "#DCE775", "#A1887F", "#F48FB1",
            "#FFE082", "#C5E1A5", "#D7CCC8", "#CE93D8", "#A5D6A7",
            "#BCAAA4", "#E6C78A")
pal_other <- setNames(pastel[seq_along(other_sp)], other_sp)

# Legend order: descending genome count (top four = major species)
sp_order <- s1 %>%
  count(species, name = "total") %>%
  arrange(desc(total), species) %>%
  filter(species %in% unique(df$species)) %>%
  pull(species)

# Abbreviated legend labels ("B. burgdorferi")
sp_labels <- setNames(paste0("B. ", sub("Borreliella ", "", sp_order)), sp_order)

# Aggregate by country x species; jitter species sharing a country
set.seed(1)
df_agg <- df %>%
  group_by(country, lon, lat, species) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(country) %>%
  mutate(lon_j = lon + runif(n(), -3, 3),
         lat_j = lat + runif(n(), -2, 2)) %>%
  ungroup()
df_agg$species <- factor(df_agg$species, levels = sp_order)

world <- map_data("world")

p_a <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey96", color = "grey85", linewidth = 0.15) +
  geom_point(data = df_agg %>% filter(!species %in% four),
             aes(lon_j, lat_j, color = species, size = n), alpha = 0.7) +
  geom_point(data = df_agg %>% filter(species %in% four),
             aes(lon_j, lat_j, color = species, size = n), alpha = 0.9) +
  scale_size_continuous(range = c(2, 9), name = "Genomes") +
  scale_color_manual(values = c(pal4, pal_other),
                     limits = sp_order,   # explicit legend order (ggplot2 >= 3.5)
                     labels = sp_labels) +
  coord_quickmap() +
  labs(color = NULL) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3), order = 1),
         size  = guide_legend(order = 2)) +
  theme_void() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.justification = "top",
        legend.text = element_text(face = "italic", size = 8))

print(p_a)
ggsave("FigS1a_map.pdf", p_a, width = 10.5, height = 5.5)