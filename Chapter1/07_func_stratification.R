# ============================================================
# Project: Borreliella pan-genome analysis
# Module: func_stratification
# Function: Functional category stratification (Core/Soft/Shell/Cloud)
# Input: ./data/gene_presence_absence_roary.csv
# Output: 
#   - ./output/pan_genome_func_annotation_rate.pdf
#   - ./output/pan_genome_func_category_stacked.pdf
#   - ./output/pan_genome_func_core_cloud_enrichment.pdf
#   - ./output/func_annotation_rate_summary.csv
#   - ./output/func_category_distribution.csv
#   - ./output/func_core_cloud_enrichment.csv
# Dependencies: data.table, ggplot2
# Note: No personal paths. Run from project root.
# ============================================================

PATH_ROARY <- "./data/gene_presence_absence_roary.csv"
PATH_OUT   <- "./output"

library(data.table)
library(ggplot2)

N_TOTAL <- 287
SAVE_OUTPUT <- TRUE

THRESH <- list(Core = c(1.0, 1.0), Soft = c(0.95, 1.0), 
               Shell = c(0.15, 0.95), Cloud = c(0, 0.15))

COL_RGB <- c("Core"  = rgb(217, 83, 79, maxColorValue = 255),
             "Soft"  = rgb(239, 172, 78, maxColorValue = 255),
             "Shell" = rgb(91, 155, 213, maxColorValue = 255),
             "Cloud" = rgb(160, 203, 232, maxColorValue = 255))

assign_strat <- function(freq) {
  ifelse(freq >= THRESH$Core[1]  & freq <= THRESH$Core[2],  "Core",
         ifelse(freq >= THRESH$Soft[1]  & freq <  THRESH$Soft[2],  "Soft",
                ifelse(freq >= THRESH$Shell[1] & freq <  THRESH$Shell[2], "Shell",
                       ifelse(freq >= THRESH$Cloud[1] & freq <  THRESH$Cloud[2], "Cloud", NA))))
}

classify_func <- function(annot) {
  annot <- tolower(as.character(annot))
  raw <- ifelse(grepl("transposase|integrase|recombinase|phage|prophage|plasmid|conjuga", annot), "Mobile Elements",
                ifelse(grepl("ospc|dbp|dbpa|dbpb|vls|vlsE|bbk|bfp|lipoprotein|outer membrane|surface|adhesin|invasin", annot), "Surface & Virulence",
                       ifelse(grepl("antigen|variable|recombination|resolvase|invertase|site-specific", annot), "Antigen Variation",
                              ifelse(grepl("ribosomal|rRNA|tRNA|translation|synthetase|polymerase|replication|helicase|topoisomerase|primase|transcription", annot), "Information Storage",
                                     ifelse(grepl("transporter|permease|abc|pts|binding-protein|import|uptake|efflux|tolC|pump", annot), "Transport Systems",
                                            ifelse(grepl("two-component|sensor|kinase|response regulator|sigma factor|transcriptional regulator|repressor|activator", annot), "Signal & Regulation",
                                                   ifelse(grepl("flagell|motility|chemotaxis|cell wall|peptidoglycan|penicillin-binding|murein|flaB|flaA|fli|flg", annot), "Cell Envelope & Motility",
                                                          ifelse(grepl("glycolysis|gluconeogenesis|sugar|carbohydrate|atp synthase|dehydrogenase|cytochrome|oxidoreductase", annot), "Energy & Carbohydrate",
                                                                 ifelse(grepl("amino acid|peptidase|protease|nucleotide|purine|pyrimidine", annot), "Amino Acid & Nucleotide",
                                                                        ifelse(grepl("coenzyme|cofactor|vitamin|biotin|lipid|fatty acid|phospholipid", annot), "Coenzyme & Lipid",
                                                                               ifelse(grepl("chaperone|clp|lon|ftsH|groEL|dnaK", annot), "Chaperones & Turnover",
                                                                                      ifelse(grepl("hypothetical|unknown|uncharacterized|putative|domain protein|conserved protein", annot), "Unknown", "Other"))))))))))))
  
  ifelse(raw %in% c("Antigen Variation", "Mobile Elements"), "Host Adaptation",
         ifelse(raw %in% c("Signal & Regulation", "Chaperones & Turnover", "Amino Acid & Nucleotide", "Coenzyme & Lipid"), "Metabolism & Regulation",
                ifelse(raw %in% c("Unknown", "Other"), "Unknown / Other", as.character(raw))))
}

# Main
roary <- fread(PATH_ROARY, header = TRUE)
roary[, freq := `No. isolates` / N_TOTAL]
roary[, strat := assign_strat(freq)]
roary[, strat := factor(strat, levels = c("Core","Soft","Shell","Cloud"))]
roary[, func_cat := sapply(Annotation, classify_func)]

annot_rate <- roary[, .(Total = .N, Annotated = sum(func_cat != "Unknown / Other"),
                        Rate = round(sum(func_cat != "Unknown / Other") / .N * 100, 2)), by = strat]

p_rate <- ggplot(annot_rate[!is.na(strat)], aes(x = strat, y = Rate, fill = strat)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = paste0(Annotated, "/", Total)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = COL_RGB) +
  scale_x_discrete(labels = c("Core" = "Core", "Soft" = "Soft core", 
                              "Shell" = "Shell", "Cloud" = "Cloud")) +
  labs(x = "Pan-genome category", y = "Functional annotation rate (%)") +
  ylim(0, max(annot_rate$Rate) * 1.15) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none", panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text.x = element_text(angle = 0, face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 12))

cat_dt <- roary[!is.na(strat), .N, by = .(strat, func_cat)]
cat_dt[, prop := N / sum(N), by = strat]
cat_levels <- c("Unknown / Other", "Host Adaptation", "Surface & Virulence", 
                "Metabolism & Regulation", "Transport Systems", "Cell Envelope & Motility",
                "Energy & Carbohydrate", "Information Storage")
cat_dt[, func_cat := factor(func_cat, levels = cat_levels)]

FUNC_COL <- c("Information Storage" = "#C75B5B", "Energy & Carbohydrate" = "#D4965C",
              "Cell Envelope & Motility" = "#6B9E75", "Transport Systems" = "#5B7E9E",
              "Metabolism & Regulation" = "#8E6B9E", "Surface & Virulence" = "#C27BA0",
              "Host Adaptation" = "#C9B85C", "Unknown / Other" = "#BBBBBB")

p_stack <- ggplot(cat_dt, aes(x = strat, y = prop, fill = func_cat)) +
  geom_bar(stat = "identity", width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = FUNC_COL) +
  scale_x_discrete(labels = c("Core" = "Core", "Soft" = "Soft core", 
                              "Shell" = "Shell", "Cloud" = "Cloud")) +
  labs(x = "Pan-genome category", y = "Proportion of genes", fill = "Functional category") +
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text.x = element_text(angle = 0, face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 12),
        legend.position = "right", legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"))

func_count <- roary[func_cat != "Unknown / Other" & !is.na(strat), .N, by = .(strat, func_cat)]
enrich <- func_count[, .(strat, func_cat, N)]
enrich[, total_strat := sum(N), by = strat]
enrich[, total_func := sum(N), by = func_cat]
enrich[, total_all := sum(N)]
enrich[, enrich_ratio := (N / total_strat) / (total_func / total_all)]

ec <- enrich[strat %in% c("Core", "Cloud")]
plot_ec <- ec
plot_ec[, func_lab := factor(func_cat, levels = rev(cat_levels[cat_levels != "Unknown / Other"]))]

p_bubble <- ggplot(plot_ec, aes(x = enrich_ratio, y = func_lab, size = N, colour = strat)) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = c("Core" = rgb(217, 83, 79, maxColorValue = 255),
                                 "Cloud" = rgb(160, 203, 232, maxColorValue = 255))) +
  scale_size_continuous(range = c(3, 10)) +
  labs(x = "Enrichment ratio", y = NULL, size = "Gene count", colour = "Stratification") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

print(p_rate); print(p_stack); print(p_bubble)

if (SAVE_OUTPUT) {
  dir.create(PATH_OUT, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(PATH_OUT, "pan_genome_func_annotation_rate.pdf"), p_rate, width = 5.5, height = 4.5, dpi = 300)
  ggsave(file.path(PATH_OUT, "pan_genome_func_category_stacked.pdf"), p_stack, width = 8, height = 5, dpi = 300)
  ggsave(file.path(PATH_OUT, "pan_genome_func_core_cloud_enrichment.pdf"), p_bubble, width = 7, height = 5, dpi = 300)
  fwrite(annot_rate, file.path(PATH_OUT, "func_annotation_rate_summary.csv"))
  fwrite(cat_dt, file.path(PATH_OUT, "func_category_distribution.csv"))
  fwrite(dcast(ec, func_cat ~ strat, value.var = "enrich_ratio", fill = 0), 
         file.path(PATH_OUT, "func_core_cloud_enrichment.csv"))
}