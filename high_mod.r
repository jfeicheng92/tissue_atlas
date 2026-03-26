library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(cowplot)
library(reshape2)
library(stringr)
library(ggbreak) 
library(patchwork)
.libPaths(c("/gpfs3/users/ludwig/cfo155/R/x86_64-pc-linux-gnu-library/4.3" , .libPaths()))
.libPaths()
options(bitmapType='cairo-png')
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls")
####Constitutive meth regions ####
prefix <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500"
# meth_mean <- fread(paste0("dmr_call_grail/",prefix,".mean.txt")) %>% as.data.frame()
# anno <- fread("resource/hg38_ws1000.s500.feature.bed", col.names = c("chr","start","end","feature"))
# anno$V1 <- paste0(anno$chr,"_", anno$start+1, "_", anno$end)
# 
# meth_mean$median_all_hmC<- apply(meth_mean %>% select(contains("_hmC")), 1, median)
# meth_mean$median_all_umC <- apply(meth_mean %>% select(contains("_umC")), 1, median)
# meth_mean$median_all_mC <- apply(meth_mean %>% select(contains("_mC")), 1, median)
# 
# bg_feature <-  anno %>%
#   pull(feature) %>%                            # Extract feature column
#   strsplit(",") %>%                            # Split by comma
#   unlist() %>%                                 # Flatten to vector
#   table() %>%                                  # Frequency table
#   as.data.frame() %>%                          # Convert to data frame for easier manipulation
#   rename(bg_Count = Freq) %>%                     # Rename columns
#   arrange(desc(bg_Count)) %>%                     # Sort by frequency
#   mutate(bg_Percentage = round(bg_Count / nrow(anno), 4)) # Calculate percentage
# colnames(bg_feature)[1] <- c("feature")
# 
# 
# sel_n <- round(nrow(meth_mean)*0.005)  
# meth_cols <- colnames(meth_mean)[grepl("mean|median", colnames(meth_mean))]
# all_feature_tables <- list()
# for (col_name in meth_cols) {
#   # Order by current column, select top sel_n
#   tmp <- meth_mean %>%
#     arrange(desc(.data[[col_name]])) %>%
#     select(V1, all_of(col_name)) %>%
#     head(sel_n)
#   
#   print(paste(col_name, tail(tmp, 1)))
#   
#   # Merge with annotation and calculate feature distribution
#   tmp_tbl <- merge(tmp, anno, by = "V1") %>%
#     pull(feature) %>%
#     strsplit(",") %>%
#     unlist() %>%
#     table() %>%
#     as.data.frame() %>%
#     rename(Count = Freq) %>%
#     arrange(desc(Count)) %>%
#     mutate(
#       Percentage = round(Count / sel_n, 4),
#       Methylation = col_name
#     )
#   
#   # Store in list
#   all_feature_tables[[col_name]] <- tmp_tbl
# }
# 
# 
# all_feature_tables <- do.call(rbind, all_feature_tables)
# colnames(all_feature_tables)[1] <- c("feature")
# fwrite(all_feature_tables, "figs/fig2/top_mod_feature_sta.csv", sep = "\t", quote = FALSE)
# fwrite(bg_feature, "figs/fig2/all_feature_sta.csv", sep = "\t", quote = FALSE)

all_feature_tables <- fread("figs/fig2/top_mod_feature_sta.csv")
bg_feature <- fread("figs/fig2/all_feature_sta.csv")
# select feature for plotting
all_features <- c("CGIshelves", "CGIshore", "cpgIsland", 
                  "GRCh38.Regulatory_Build.CTCF_binding_site", "GRCh38.Regulatory_Build.enhancer", "GRCh38.Regulatory_Build.open_chromatin_region", "GRCh38.Regulatory_Build.promoter", "GRCh38.Regulatory_Build.promoter_flanking_region", "GRCh38.Regulatory_Build.TF_binding_site", 
                  "hg38.repeatmasker", "hg38.repeatmasker.DNA_hAT", "hg38.repeatmasker.DNA_TcMar", "hg38.repeatmasker.LINE_CR1", "hg38.repeatmasker.LINE_L1", "hg38.repeatmasker.LINE_L2", "hg38.repeatmasker.LINE_RTE", "hg38.repeatmasker.LTR_ERV", "hg38.repeatmasker.Retroposon_SVA", "hg38.repeatmasker.Satellite_centromeric", "hg38.repeatmasker.SINE_Alu", "hg38.repeatmasker.SINE_MIR", 
                  "hg38lift_genome_100_segments/Acet", "hg38lift_genome_100_segments/BivProm", "hg38lift_genome_100_segments/DNase", "hg38lift_genome_100_segments/EnhA", "hg38lift_genome_100_segments/EnhWk", "hg38lift_genome_100_segments/GapArtf", "hg38lift_genome_100_segments/HET", "hg38lift_genome_100_segments/hg38lift_genome_100_segments.Enh", "hg38lift_genome_100_segments/hg38lift_genome_100_segments.Prom", "hg38lift_genome_100_segments/PromF", "hg38lift_genome_100_segments/Quies", "hg38lift_genome_100_segments/ReprPC", "hg38lift_genome_100_segments/TSS", "hg38lift_genome_100_segments/Tx", "hg38lift_genome_100_segments/TxEnh", "hg38lift_genome_100_segments/TxEx", "hg38lift_genome_100_segments/TxWk", "hg38lift_genome_100_segments/znf", 
                  "Homo_sapiens.GRCh38.regulatory_features.v114.CTCF_binding_site", "Homo_sapiens.GRCh38.regulatory_features.v114.enhancer", "Homo_sapiens.GRCh38.regulatory_features.v114.open_chromatin_region", "Homo_sapiens.GRCh38.regulatory_features.v114.promoter", 
                  "lincRNA", 
                  "MANE.GRCh38.v1.0.refseq_genomic.exon", "MANE.GRCh38.v1.0.refseq_genomic.gene", "MANE.GRCh38.v1.0.refseq_genomic.gene.down1k", "MANE.GRCh38.v1.0.refseq_genomic.gene.up1k", "MANE.GRCh38.v1.0.refseq_genomic.gene.up2k", "MANE.GRCh38.v1.0.refseq_genomic.intron",
                  "PMD_coordinates_hg38.commonPMD")
select_features <- c("Homo_sapiens.GRCh38.regulatory_features.v114.enhancer", "Homo_sapiens.GRCh38.regulatory_features.v114.promoter", "Homo_sapiens.GRCh38.regulatory_features.v114.CTCF_binding_site",  "Homo_sapiens.GRCh38.regulatory_features.v114.open_chromatin_region", 
                     "cpgIsland","CGIshore", "CGIshelves",  "MANE.GRCh38.v1.0.refseq_genomic.gene",
                     "hg38.repeatmasker.DNA_hAT", "hg38.repeatmasker.DNA_TcMar", "hg38.repeatmasker.LINE_CR1", "hg38.repeatmasker.LINE_L1", "hg38.repeatmasker.LINE_L2", "hg38.repeatmasker.LINE_RTE", "hg38.repeatmasker.LTR_ERV", "hg38.repeatmasker.SINE_Alu", "hg38.repeatmasker.SINE_MIR"
)

rownames(all_feature_tables) <- NULL

features_keep <- c(
  "Homo_sapiens.GRCh38.regulatory_features.v114.enhancer",
  "Homo_sapiens.GRCh38.regulatory_features.v114.promoter",
  "Homo_sapiens.GRCh38.regulatory_features.v114.CTCF_binding_site",
  "Homo_sapiens.GRCh38.regulatory_features.v114.open_chromatin_region",
  "cpgIsland", "CGIshore", "CGIshelves",
  "MANE.GRCh38.v1.0.refseq_genomic.gene",
  "LINE", "SINE", "LTR", "DNA"
)
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "Liver-Cirrhosis","Pancreas-Pancreatitis",
                  "Brain-Tumour", "Breast-Tumour", "Colon-Tumour", "Kidney-Tumour", "Liver-Tumour", "Lung-Tumour", "Ovary-Tumour", "Pancreas-Tumour", "Prostate-Tumour","Stomach-Tumour",
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")

# Build category, then sum Count and Percentage within each category
summary_tbl <- all_feature_tables %>%
  mutate(
    category = case_when(
      str_detect(feature, "^hg38\\.repeatmasker\\.LINE_") ~ "LINE",
      str_detect(feature, "^hg38\\.repeatmasker\\.SINE_") ~ "SINE",
      str_detect(feature, "^hg38\\.repeatmasker\\.LTR_")  ~ "LTR",
      str_detect(feature, "^hg38\\.repeatmasker\\.DNA_")  ~ "DNA",
      TRUE ~ feature
    )
  ) %>%
  filter(category %in% features_keep) %>%
  group_by(Methylation, category) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    Percentage = sum(Percentage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(category = factor(category, levels = features_keep)) %>%
  arrange(Methylation, category)

bg_tbl <- bg_feature %>%
  mutate(
    category = case_when(
      str_detect(feature, "^hg38\\.repeatmasker\\.LINE_") ~ "LINE",
      str_detect(feature, "^hg38\\.repeatmasker\\.SINE_") ~ "SINE",
      str_detect(feature, "^hg38\\.repeatmasker\\.LTR_")  ~ "LTR",
      str_detect(feature, "^hg38\\.repeatmasker\\.DNA_")  ~ "DNA",
      TRUE ~ feature
    )
  ) %>%
  filter(category %in% features_keep) %>%
  group_by(category) %>%
  summarise(
    Count = sum(bg_Count , na.rm = TRUE),
    Percentage = sum(bg_Percentage, na.rm = TRUE),
    .groups = "drop"
  )

summary_tbl <- summary_tbl %>%
  left_join(bg_tbl, by = c("category"))
summary_tbl$odds <- summary_tbl$Percentage.x / summary_tbl$Percentage.y
summary_tbl$Methylation <- gsub("CD34-erythroblasts","Erythroid-precursors",summary_tbl$Methylation)
summary_tbl$Methylation <- gsub("CD34-megakaryocytes","Megakaryocytes",summary_tbl$Methylation)
summary_tbl <-summary_tbl %>%
  as.data.frame() %>%
  separate(
    Methylation,
    into = c("stat", "tissue", "mark"),             
    sep = "_",
    extra = "merge",                                  
    fill  = "right"
  )%>%
  mutate(
    tissue=factor(tissue,levels=c(tissue_order,"all")),
    category=factor(category, levels=rev(features_keep))
  )%>%
  mutate(
    type = "solid",
    type = ifelse(grepl("Tumour", tissue), "tumour", type),
    type = ifelse(grepl("-Pancreatitis|-Cirrhosis", tissue), "pre", type),
    type = ifelse(tissue %in% c(
      "CD4-T-cells","CD8-T-cells","NK-cells","B-cells",
      "Neutrophils","Eosinophils","Monocytes",
      "Erythroid-precursors","Megakaryocytes"
    ), "blood", type),
    type = factor(type, levels = c("blood","pre","tumour","solid")) # optional order
  )




base_plot <- function(mark_label) {
  summary_tbl %>%
    filter(mark == mark_label, type %in% c("blood","solid")) %>%
    mutate(category = factor(category, levels = rev(features_keep))) %>%
    ggplot(aes(x = odds, y = category)) +
    geom_boxplot(aes(y=category),outlier.shape = NA, width = 0.6) +
    geom_jitter(
      aes(color = type),
      size = 0.6, alpha = 0.6
    ) +
    scale_color_manual(values = c(blood = "#6BAED6", solid = "#E34A33")) +
    theme_light() +
    theme(
      legend.position = "bottom",
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# Common x-scale (limits + breaks) and the same break for all three
common_x <- scale_x_continuous(
  limits = c(0, 60),
  breaks  = c(0, 2, 4, 6, 15, 37.5, 60),
  expand  = expansion(mult = c(0, 0.02))
)
xbreak <- scale_x_break(c(6, 12), scales = 0.2)

p1 <- base_plot("hmC") + ggtitle("hmC enrich score") + common_x + xbreak
p2 <- base_plot("mC")  + ggtitle("mC enrich score")  + common_x + xbreak
p3 <- base_plot("umC") + ggtitle("umC enrich score") + common_x + xbreak

# Equal widths; one shared legend
pdf("figs/fig2/top_mod_feature_enrich.pdf", width = 10, height = 6)
(p3 | p1 | p2) +
  plot_layout(widths = c(2, 1, 0.2)) &
  theme(legend.position = "bottom")
p1+theme(
  legend.position = "bottom",
  axis.text.y  = element_text()
)
dev.off()

