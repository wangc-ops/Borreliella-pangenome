# ============================================================
# Project: Borreliella pan-genome analysis
# Module: MGE_stacked_bar (CORRECTED)
# Function: MGE distribution 100% stacked bar + Table S3
# Input: ./data/gene_presence_absence_roary.csv
# Output: ./output/MGE_stacked_bar.pdf, ./output/TableS3_MGE_distribution.csv
# Dependencies: data.table, ggplot2
# Note: No personal paths. Run from project root.
##
# ============================================================

PATH_ROARY <- "./data/gene_presence_absence_roary.csv"
PATH_OUT   <- "./output"

library(data.table)
library(ggplot2)

N_TOTAL <- 287
SAVE_OUTPUT <- TRUE

THRESH <- list(Core = c(1.0, 1.0), Soft = c(0.95, 1.0), Shell = c(0.15, 0.95), Cloud = c(0, 0.15))

assign_strat <- function(freq) {
  ifelse(freq >= THRESH$Core[1]  & freq <= THRESH$Core[2],  "Core",
         ifelse(freq >= THRESH$Soft[1]  & freq <  THRESH$Soft[2],  "Soft",
                ifelse(freq >= THRESH$Shell[1] & freq <  THRESH$Shell[2], "Shell",
                       ifelse(freq >= THRESH$Cloud[1] & freq <  THRESH$Cloud[2], "Cloud", NA))))
}

classify_mge <- function(annot) {
  annot <- tolower(as.character(annot))
  ifelse(
    grepl("vls|ospc|dbp|bbk|cp32|cp26|lp54|lp28", annot), "Borrelia_plasmid_virulence",
    ifelse(grepl("\\btra[a-z]\\b|\\btra[0-9]\\b|\\btrb[a-z]\\b|\\btrb[0-9]\\b|\\bmob[a-z]\\b|\\bmob[0-9]\\b|relaxase|type iv secretion|conjugat", annot), "Conjugation",
           ifelse(grepl("phage|prophage|tail|capsid|lysin|holin|portal|terminase", annot), "Phage",
                  ifelse(grepl("transposase|transposon|insertion sequence", annot), "Transposase",
                         ifelse(grepl("integrase|recombinase|resolvase|invertase|tyrosine recombinase", annot), "Integrase_recombinase",
                                ifelse(grepl("plasmid|partition|par\\s|rep\\s|replication protein", annot), "Plasmid_maintenance",
                                       ifelse(grepl("reverse transcriptase|group ii intron|crispr|\\bcas[0-9]?\\b", annot), "Other_MGE", "Non-MGE"))))))
  )
}

# Main
roary <- fread(PATH_ROARY, header = TRUE)
roary[, freq := `No. isolates` / N_TOTAL]
roary[, strat := assign_strat(freq)]
roary[, strat := factor(strat, levels = c("Core","Soft","Shell","Cloud"))]
roary[, mge_type := classify_mge(Annotation)]

cat("=== MGE classification summary (full pan-genome, n =", nrow(roary), ") ===\n")
print(roary[, .N, by = mge_type][order(-N)])

# Table S3 (full pan-genome)
mge_tab <- roary[, .N, by = .(strat, mge_type)]
mge_tab[, prop := round(N / sum(N) * 100, 2), by = strat]

# Merge for plot
mge_tab[, mge_merge := fcase(
  mge_type == "Non-MGE", "Non-MGE",
  mge_type == "Plasmid_maintenance", "Plasmid maintenance",
  mge_type == "Phage", "Phage",
  mge_type == "Borrelia_plasmid_virulence", "Borrelia virulence",
  default = "Other MGE"
)]

mge_merge <- mge_tab[, .(N = sum(N)), by = .(strat, mge_merge)]
mge_merge[, prop := N / sum(N), by = strat]
mge_merge[, strat := factor(strat, levels = c("Core", "Soft", "Shell", "Cloud"))]
mge_merge[, mge_merge := factor(mge_merge, levels = c("Non-MGE", "Plasmid maintenance", "Phage", "Borrelia virulence", "Other MGE"))]

MGE_COL <- c("Non-MGE" = "#CCCCCC", "Plasmid maintenance" = "#5B9BD5",
             "Phage" = "#984EA3", "Borrelia virulence" = "#D9534F", "Other MGE" = "#F0AD4E")

p_stack <- ggplot(mge_merge, aes(x = strat, y = prop, fill = mge_merge)) +
  geom_bar(stat = "identity", width = 0.7, colour = "white", size = 0.3) +
  geom_text(aes(label = ifelse(prop > 0.03, paste0(round(prop*100, 1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = MGE_COL) +
  scale_x_discrete(labels = c("Core" = "Core", "Soft" = "Soft core", "Shell" = "Shell", "Cloud" = "Cloud")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pan-genome category", y = "Proportion", fill = "MGE category") +
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(angle = 0, face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 12),
        legend.position = "right", legend.title = element_text(face = "bold"))

print(p_stack)

if (SAVE_OUTPUT) {
  dir.create(PATH_OUT, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(PATH_OUT, "MGE_stacked_bar.pdf"), p_stack, width = 6.5, height = 4.5, dpi = 300)
  
  s3_cast <- dcast(mge_tab, mge_type ~ strat, value.var = c("N", "prop"), fill = 0)
  s3_final <- data.table(
    MGE_type = s3_cast$mge_type,
    Core = paste0(s3_cast$N_Core, "/", s3_cast$prop_Core),
    `Soft core` = paste0(s3_cast$N_Soft, "/", s3_cast$prop_Soft),
    Shell = paste0(s3_cast$N_Shell, "/", s3_cast$prop_Shell),
    Cloud = paste0(s3_cast$N_Cloud, "/", s3_cast$prop_Cloud)
  )
  fwrite(s3_final, file.path(PATH_OUT, "TableS3_MGE_distribution.csv"))
  cat("Table S3 saved.\n")
  print(s3_final)
}

cat("\nDone.\n")