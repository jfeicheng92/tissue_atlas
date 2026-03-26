library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(cowplot)
library(reshape2)
library(stringr)
library(patchwork)
options(bitmapType='cairo-png')
.libPaths(c(.libPaths(),"/gpfs3/users/ludwig/cfo155/R/x86_64-pc-linux-gnu-library/4.3"))
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls")

#### Load methylation data ####
prefix <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500"
meth <- fread(paste0(prefix, ".bed"), header=TRUE) %>% as.data.frame()

tissue_order <- c( "Brain", "Brain.Tumour","Breast","Breast.Tumour","Heart", "Esophagus", "Colon","Colon.Tumour",
                   "Kidney", "Kidney.Tumour","Liver","Liver.Tumour","Liver.Cirrhosis", "Lung","Lung.Tumour","Ovary", "Ovary.Tumour",
                   "Pancreas", "Pancreas.Tumour","Pancreas.Pancreatitis","Prostate","Prostate.Tumour","Stomach","Stomach.Tumour",
                   "Spleen", "CD4.T.cells", "CD8.T.cells", "Neutrophils", "NK.cells", "B.cells", "Eosinophils", "Monocytes",
                   "CD34.erythroblasts", "CD34.megakaryocytes") %>% gsub("\\.","-",.)

mods <- c("mC", "hmC","umC")

#### Mean methy per tissue ####
for (tissue in tissue_order) {
  for (mod in mods) {
    pattern <- paste0(tissue, "_", mod)
    col_name <- paste0("mean_", tissue, "_", mod)
    print(col_name)

    if(sum(grepl(pattern, names(meth)))>1){
      meth[[col_name]] <- rowMeans(meth[, grepl(pattern, names(meth))], na.rm = TRUE)
    }
    else
      meth[[col_name]] <- meth[, grepl(pattern, names(meth))]
  }
}


rownames(meth) <- paste(meth$chr,meth$start+1,meth$end,sep="_")
meth <- meth %>%
  .[complete.cases(.), ]

meth %>%
  select(contains("mean")) %>%
  write.table(.,paste0(prefix,".mean.txt"),sep="\t", row.names = TRUE)

#### Correlation ####
prefix <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500"
meth_mean <- fread(paste0("dmr_call_grail/",prefix,".mean.txt")) %>% as.data.frame()
anno <- fread("resource/hg38_ws1000.s500.feature.bed", col.names = c("chr","start","end","feature"))
anno$pos <- paste0(anno$chr,"_",anno$start+1, "_", anno$end)
features_keep <- c(
  "Homo_sapiens.GRCh38.regulatory_features.v114.enhancer",
  "Homo_sapiens.GRCh38.regulatory_features.v114.promoter",
  "Homo_sapiens.GRCh38.regulatory_features.v114.CTCF_binding_site",
  "Homo_sapiens.GRCh38.regulatory_features.v114.open_chromatin_region",
  "cpgIsland", "CGIshore", "CGIshelves",
  "MANE.GRCh38.v1.0.refseq_genomic.gene",
  "LINE", "SINE", "LTR", "DNA"
)
results <- list()
for (sel_feature in features_keep) {
  meth_mean_sel <- anno %>%
    filter(grepl(sel_feature, feature)) %>%
    select(pos) %>%
    left_join(meth_mean, by = c("pos" = "V1")) %>%
    pivot_longer(
      cols = starts_with("mean_"),
      names_to = c("tissue", "mark"),
      names_pattern = "mean_(.*)_(mC|hmC|umC)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = mark, values_from = value) %>%
    group_by(tissue)
  
  meth_mean_sel_sta <- meth_mean_sel %>%
    summarise(
      cor_mC_hmC  = cor(mC, hmC,  use = "pairwise.complete.obs"),
      cor_mC_umC  = cor(mC, umC,  use = "pairwise.complete.obs"),
      cor_hmC_umC = cor(hmC, umC, use = "pairwise.complete.obs"),
      mean_mC  = mean(mC,  na.rm = TRUE),
      mean_hmC = mean(hmC, na.rm = TRUE),
      mean_umC = mean(umC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(feature = sel_feature)
  
  results[[sel_feature]] <- meth_mean_sel_sta
}


meth_mean_sta <- bind_rows(results)
write.table(meth_mean_sta, paste0("figs/fig3/",prefix,".summary.txt"),sep="\t", quote = FALSE)


meth_mean_sta <- fread(paste0("figs/fig3/",prefix,".summary.txt"))
print(dim(meth_mean_sta))
meth_mean_sta$tissue <- gsub("CD34-erythroblasts","Erythroid-precursors",meth_mean_sta$tissue)
meth_mean_sta$tissue <- gsub("CD34-megakaryocytes","Megakaryocytes",meth_mean_sta$tissue)

meth_mean_sta <- meth_mean_sta %>% filter(feature %in% features_keep)

meth_mean_sta <- meth_mean_sta %>%
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
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "Liver-Cirrhosis","Pancreas-Pancreatitis",
                  "Brain-Tumour", "Breast-Tumour", "Colon-Tumour", "Kidney-Tumour", "Liver-Tumour", "Lung-Tumour", "Ovary-Tumour", "Pancreas-Tumour", "Prostate-Tumour","Stomach-Tumour",
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")

meth_mean_sta$tissue <- factor(meth_mean_sta$tissue, levels = tissue_order)
meth_mean_sta <- meth_mean_sta %>% mutate(tissue = factor(tissue, levels = tissue_order))

base_plot <- function(cor_type) {
  stopifnot(cor_type %in% names(meth_mean_sta))
  meth_mean_sta %>%
    filter(type %in% c("solid","blood"),
           feature %in% features_keep) %>%
    mutate(feature = factor(feature, levels = rev(features_keep))) %>%
    ggplot(aes(x = .data[[cor_type]], y = feature)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
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
    ) +
    xlim(-1,1) +
    geom_vline(xintercept = 0)
}

# If you want the three pairwise correlations:
p1 <- base_plot("cor_mC_umC") + ggtitle("mC umC correlation")
p2 <- base_plot("cor_mC_hmC") + ggtitle("mC hmC correlation")
p3 <- base_plot("cor_hmC_umC") + ggtitle("hmC umC correlation")
# Equal widths; one shared legend
pdf("figs/fig3/feature_cor.pdf", width = 10, height = 6)
(p1 | p2 | p3) +
  plot_layout(widths = c(1,1,1)) 
p1+theme(
  axis.text.y  = element_text()
)
dev.off()

meth_mean_sta %>%
  select(cor_mC_umC) %>%
  apply(2,mean)

meth_mean_sta %>%
  group_by(feature) %>%
  summarise(
    across(contains("cor"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    .groups = "drop"
  )
# # A tibble: 12 × 4
# feature                                                            mean_cor_mC_hmC mean_cor_mC_umC mean_cor_hmC_umC
# <chr>                                                                        <dbl>           <dbl>            <dbl>
# 1 CGIshelves                                                                 -0.0895          -0.985         -0.0514 
# 2 CGIshore                                                                    0.209           -0.994         -0.296  
# 3 DNA                                                                        -0.219           -0.976          0.0397 
# 4 Homo_sapiens.GRCh38.regulatory_features.v114.CTCF_binding_site              0.232           -0.994         -0.319  
# 5 Homo_sapiens.GRCh38.regulatory_features.v114.enhancer                       0.0481          -0.987         -0.182  
# 6 Homo_sapiens.GRCh38.regulatory_features.v114.open_chromatin_region          0.0200          -0.991         -0.125  
# 7 Homo_sapiens.GRCh38.regulatory_features.v114.promoter                       0.464           -0.995         -0.529  
# 8 LINE                                                                       -0.201           -0.976          0.0238 
# 9 LTR                                                                        -0.194           -0.977          0.0204 
# 10 MANE.GRCh38.v1.0.refseq_genomic.gene                                       -0.151           -0.982         -0.00554
# 11 SINE                                                                       -0.228           -0.977          0.0540 
# 12 cpgIsland                                                                   0.395           -0.997         -0.449 
#### dmr quant parameters ####
bg_tg_q <- "bgQ0.05-0.05.tgQ0.25"
#### hyper vs. hypo DMR ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/") ## using quantiles of groups instead of individual samples
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")

# --- inputs ---
dmr_files <- c(
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"))
selN <- 1000

# optional: set your desired tissue order
# tissue_order <- c("Brain","Liver","Lung",...) 

mark_from_path <- function(p) {
  if (grepl("\\.hmC\\.", p)) "5hmC"
  else if (grepl("\\.mC\\.", p)) "5mC"
  else if (grepl("\\.umC\\.", p)) "5umC"
  else "mark"
}

read_pair <- function(hypo_path) {
  hyper_path <- sub("hypo_dmrs", "hyper_dmrs", hypo_path)
  mark <- mark_from_path(hypo_path)
  
  dt_hypo <- fread(hypo_path, showProgress = FALSE)
  dt_hyper <- fread(hyper_path, showProgress = FALSE)
  
  # keep only needed columns if they exist
  keep <- c("chr","start","end","delta_quants","Hyper_or_Hypo","selected_tissue")
  dt <- rbindlist(list(dt_hypo[, ..keep], dt_hyper[, ..keep]), use.names = TRUE, fill = TRUE)
  dt[, mark := mark]
  dt
}

dmr_all <- purrr::map_dfr(dmr_files, read_pair) %>%
  mutate(
    Hyper_or_Hypo = factor(Hyper_or_Hypo, levels = c("Hypo","Hyper"))
  )

# choose top selN per tissue *within each mark*
plot_df <- dmr_all %>%
  group_by(mark, selected_tissue) %>%
  arrange(desc(delta_quants), .by_group = TRUE) %>%
  slice_head(n = selN) %>%
  count(mark, selected_tissue, Hyper_or_Hypo, name = "n") %>%
  ungroup()

plot_df <- plot_df[plot_df$selected_tissue %in% tissue_order,]
plot_df$selected_tissue <- factor(plot_df$selected_tissue, levels = rev(tissue_order))
plot_df$mark <- factor(plot_df$mark, levels=c("5mC","5umC","5hmC"))

# stacked bars, faceted by mark
p <- ggplot(plot_df, aes(y = selected_tissue, x = n, fill = Hyper_or_Hypo)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ mark, nrow = 1, scales = "free_x") +
  labs(
    y = "Tissue",
    x = paste0("Number of top ", selN, " DMRs"),
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("#0571B0","#CA0020")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(1, "lines")
  )
p
ggsave("../figs/fig3/topDMR_hyper_hypo.pdf",p, width = 6, height = 6)


#### DMR heatmap ####
setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls")
dir.create("figs/fig3/heatmap")
plot_heatmap <- function(dmr_file, all_meth_file, tissue_order, dmr_type ="healthy", selN=100, width = 12,height = 7,prefix){
  ## select top N dmr for each tissue ##
  dmr <- fread(dmr_file) %>% 
    group_by(selected_tissue) %>%
    arrange(desc(delta_quants)) %>%
    slice_head(n = selN) %>% 
    select("chr","start","end","selected_tissue") %>%
    as.data.frame()
  dmr$tissue <- as.factor(dmr$selected_tissue)
  
  ## load meth ##
  all_meth <- fread(all_meth_file)
  colnames(all_meth) <- gsub("CD34-erythroblasts","Erythroid-precursors",colnames(all_meth))
  colnames(all_meth) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(all_meth))
  meth <- merge(dmr,all_meth,by=c("chr","start","end"))
  
  if(dmr_type == "healthy"){
    tissue_order <- grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,invert = TRUE, value = TRUE,
                         grep("Tumour|Cirrhosis|Pancreatitis",tissue_order, value = TRUE))
    meth <- meth[grep("Tumour|Cirrhosis|Pancreatitis",meth$selected_tissue,invert = TRUE),]
    meth$tissue <- factor(meth$tissue,levels=grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,invert = TRUE, value = TRUE))
  }else if(dmr_type == "tumour"){
    tissue_order <- c(grep("Tumour|Cirrhosis|Pancreatitis",tissue_order, value = TRUE),
                      grep("Tumour|Cirrhosis|Pancreatitis",tissue_order, invert = TRUE, value = TRUE))
    meth <- meth[grep("Tumour|Cirrhosis|Pancreatitis",meth$selected_tissue),]
    meth$tissue <- factor(meth$tissue,levels=grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,value = TRUE))
  }else if(dmr_type=="all"){
    tissue_order <- c(grep("Tumour|Cirrhosis|Pancreatitis",tissue_order, value = TRUE),
                      grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,invert = TRUE, value = TRUE))
    meth$tissue <- factor(meth$tissue,levels=tissue_order)
  }
    
  anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(meth$tissue))))
  names(anno_colors) <- unique(meth$tissue)
  anno_colors <- list(type=anno_colors)
  anno_row <- data.frame(type=meth$tissue); rownames(anno_row)<- rownames(meth)
  
  meth <- meth[order(meth$tissue),]
  
  
  
  
  smp_order <- data.frame(smp=colnames(meth)[-c(1:5)])
  smp_order$tissue <- factor(gsub("_hmC|_umC|_mC","",smp_order$smp) %>% gsub(".*_","",.), levels=tissue_order)
  smp_order <- smp_order[order(smp_order$tissue),]
  smp_order <- smp_order[complete.cases(smp_order),]
  
  print(dim(meth))
  
  # extract data subsets
  dmr_umC <- meth[meth$tissue %in% smp_order$tissue,
                  c(1:5, (smp_order[grep("_umC", smp_order$smp),] %>% rownames() %>% as.numeric())+5)]
  
  dmr_mC  <- meth[meth$tissue %in% smp_order$tissue,
                  c(1:5, (smp_order[grep("_mC", smp_order$smp),] %>% rownames() %>% as.numeric())+5)]
  
  dmr_hmC <- meth[meth$tissue %in% smp_order$tissue,
                  c(1:5, (smp_order[grep("_hmC", smp_order$smp),] %>% rownames() %>% as.numeric())+5)]
  
  
  # --- Heatmaps ---
  
  # umC 
  pheatmap(dmr_umC[,-c(1:5)], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
           filename = paste0(prefix,".umC.label.png"), main=paste0(prefix,".umC"),
           width = width, height = height)
  
  pheatmap(dmr_umC[,-c(1:5)], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
           width = width, height = height,
           annotation_colors=anno_colors, annotation_row = anno_row, 
           labels_col= colnames(dmr_umC)[-c(1:5)] %>% gsub("_hmC|_umC|_mC","",.) %>% gsub(".*_","",.), 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           filename = paste0(prefix,".umC.png"))
  
  
  # mC
  pheatmap(dmr_mC[,-c(1:5)], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
           filename = paste0(prefix,".mC.label.png"), main=paste0(prefix,".mC"),
           width = width, height = height)
  
  pheatmap(dmr_mC[,-c(1:5)], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
           width = width, height = height,
           annotation_colors=anno_colors, annotation_row = anno_row, 
           labels_col= colnames(dmr_mC)[-c(1:5)] %>% gsub("_hmC|_umC|_mC","",.) %>% gsub(".*_","",.), 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           filename = paste0(prefix,".mC.png"))
  
  
  # hmC
  pheatmap(dmr_hmC[,-c(1:5)], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
           filename = paste0(prefix,".hmC.label.png"), main=paste0(prefix,".hmC"),
           width = width, height = height)
  
  pheatmap(dmr_hmC[,-c(1:5)], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
           width = width, height = height,
           annotation_colors=anno_colors, annotation_row = anno_row, 
           labels_col= colnames(dmr_hmC)[-c(1:5)] %>% gsub("_hmC|_umC|_mC","",.) %>% gsub(".*_","",.), 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           filename = paste0(prefix,".hmC.png"))
  
  
}



tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Monocytes", "Eosinophils", "Erythroid-precursors", "Megakaryocytes")


dmr_files <- c(
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt")
  )
selN <- 200
for(dmr_file in dmr_files){
  plot_heatmap(dmr_file = paste0("dmr_call_new1/",dmr_file),
               all_meth_file = "dmr_call_new1/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.bed",
               tissue_order=tissue_order,
               selN=selN, width = 10, height = 7, dmr_type = "healthy",
               prefix = paste0("figs/fig3/heatmap/",gsub("txt$",selN, dmr_file))
               
  )
}

tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes",
                  "Brain-Tumour", "Breast-Tumour", "Colon-Tumour", "Kidney-Tumour", "Liver-Tumour", "Lung-Tumour", "Ovary-Tumour", "Pancreas-Tumour", "Prostate-Tumour","Stomach-Tumour")
for(dmr_file in dmr_files){
  plot_heatmap(dmr_file = paste0("dmr_call_new1/",dmr_file),
               all_meth_file = "dmr_call_new1/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.bed",
               tissue_order=tissue_order,
               selN=selN, width = 12, height = 7, dmr_type = "all",
               prefix = paste0("figs/fig3/heatmap/",gsub("txt$",selN, dmr_file))
               
  )
}

#### DMR Feature enrichment ####
# run bash script /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/code_tissue_methylome_atlas/ws_1kb_dmr.sh
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/")
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes",
                  "Brain-Tumour", "Breast-Tumour", "Colon-Tumour", "Kidney-Tumour", "Liver-Tumour", "Lung-Tumour", "Ovary-Tumour", "Pancreas-Tumour", "Prostate-Tumour","Stomach-Tumour")
all_feature_file <- c("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmr_call_new1/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.feature.txt")

dmr_files <- c(
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"))

feature_enrich_files <- c("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.feature.summary.txt",
                          gsub(".txt",".feature.summary.txt", dmr_files))

dat <- NULL  # initialize empty data.table

for (file in feature_enrich_files) {
  tmp <- fread(file)
  dat <- rbind(dat, tmp)
}


dat <- as.data.frame(dat)
colnames(dat) <- c("region","feature","num","pct")
sel_df <- dat %>%
  select(region, feature, pct) %>%
  pivot_wider(names_from = feature, values_from = pct) %>% as.data.frame()
rownames(sel_df)<- sel_df$region

sel_df$region <- factor(gsub(".*_","",sel_df$region),levels=c("windows",tissue_order))
sel_df <- sel_df[order(sel_df$region),]
sel_df$region <- NULL

normalized_df <- sel_df
for (col in colnames(sel_df)) {
  normalized_df[[col]] <- sel_df[[col]] / sel_df[[col]][1]
}
select_features <- c("GRCh38.Regulatory_Build.CTCF_binding_site", "GRCh38.Regulatory_Build.enhancer", "GRCh38.Regulatory_Build.open_chromatin_region", "GRCh38.Regulatory_Build.promoter", "GRCh38.Regulatory_Build.promoter_flanking_region", "GRCh38.Regulatory_Build.TF_binding_site",
                     "lincRNA",  "MANE.GRCh38.v1.0.refseq_genomic.gene", "PMD_coordinates_hg38.commonPMD",
                     "hg38.repeatmasker.DNA_hAT", "hg38.repeatmasker.DNA_TcMar", "hg38.repeatmasker.LINE_CR1", "hg38.repeatmasker.LINE_L1", "hg38.repeatmasker.LINE_L2", "hg38.repeatmasker.LINE_RTE", "hg38.repeatmasker.LTR_ERV", "hg38.repeatmasker.SINE_Alu", "hg38.repeatmasker.SINE_MIR"
)
normalized_df_sel <- normalized_df %>%
  select(all_of(select_features))
rownames(normalized_df_sel) <- gsub("bg.*groups.","",rownames(normalized_df_sel))
my_colors <- colorRampPalette(c("white", "red"))(100)  # adjust colors as needed
breaks <- seq(1, 3, length.out = length(my_colors) + 1)

pdf("../figs/fig3/1kb_top2000_dmr_feature_distribution.pdf", width = 10, height = 8)
pheatmap(normalized_df_sel[grepl("^hmC.*Hyper", rownames(normalized_df_sel)),], main="hmC_hyper", cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6, color = my_colors,breaks = breaks)
pheatmap(normalized_df_sel[grepl("^umC.*Hyper", rownames(normalized_df_sel)),], main="umC_hyper", cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6, color = my_colors,breaks = breaks)
pheatmap(normalized_df_sel[grepl("^mC.*Hypo", rownames(normalized_df_sel)),], main="mC_hypo",cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6, color = my_colors,breaks = breaks)
pheatmap(normalized_df_sel[grepl("^hmC.*Hypo", rownames(normalized_df_sel)),], main="hmC_hypo", cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6, color = my_colors,breaks = breaks)
pheatmap(normalized_df_sel[grepl("^umC.*Hypo", rownames(normalized_df_sel)),], main="umC_hypo", cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6, color = my_colors,breaks = breaks)
pheatmap(normalized_df_sel[grepl("^mC.*Hyper", rownames(normalized_df_sel)),], main="mC_hyper", cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6, color = my_colors,breaks = breaks)
dev.off()

#### Venn diagram ####
library(VennDiagram)
library(grid)
library(eulerr)
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/")
selN<-1000
select_tissues <- c(
  "Brain","Breast","Heart","Kidney","Liver","Lung","Ovary","Pancreas",
  "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells","CD8-T-cells",
  "NK-cells","B-cells","Neutrophils","Eosinophils","Monocytes",
  "Erythroid-precursors","Megakaryocytes"
)


hmC_hyper_dmr <- fread("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.all.txt") %>%
  filter(selected_tissue %in% select_tissues) %>%
  as_tibble() %>% select(chr, start, end, selected_tissue,delta_quants) %>%
  group_by(selected_tissue) %>%
  arrange(desc(delta_quants)) %>%
  slice_head(n = selN) %>%
  select(chr,start,end, selected_tissue) 


umC_hyper_dmr <- fread("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.all.txt") %>%
  filter(selected_tissue %in% select_tissues) %>%
  as_tibble() %>% select(chr, start, end, selected_tissue,delta_quants) %>%
  group_by(selected_tissue) %>%
  arrange(desc(delta_quants)) %>%
  slice_head(n = selN) %>%
  select(chr,start,end, selected_tissue) 

mC_hypo_dmr <- fread("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.bgQ0.05-0.05.tgQ0.25.hypo_dmrs.bg_quant_modegroups.all.txt") %>%
  filter(selected_tissue %in% select_tissues) %>%
  as_tibble() %>% select(chr, start, end, selected_tissue,delta_quants) %>%
  group_by(selected_tissue) %>%
  arrange(desc(delta_quants)) %>%
  slice_head(n = selN) %>%
  select(chr,start,end, selected_tissue) 

# helper: filter tissues and make an ID per interval
to_ids <- function(df) {
  df %>%
    filter(selected_tissue %in% select_tissues) %>%
    transmute(id = paste0(chr, ":", start, "-", end)) %>%
    pull(id) %>%
    unique()
}

dmr_list <- list(
  `5hmC hyper` = to_ids(hmC_hyper_dmr),
  `5umC hyper` = to_ids(umC_hyper_dmr),
  `5mC hypo`   = to_ids(mC_hypo_dmr)
)


pdf("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig3/DMR_overlap.venn.pdf", width = 6, height =6)
plot(euler(dmr_list, shape = "ellipse"), fills = list(
  fill  = c("#C03B2F", "#373c95", "#76b7b3"),  # one colour per set
  edges = list(col = "#ebebeb", lwd = 2, lty = 1, alhpa=0.1),
  alpha = 0.6
),quantities = TRUE)
dev.off()


#### ranksum test ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/")

dmrp <- fread("dmr_rank/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.wilcox_results.csv")

# Standardize tissue labels in column names
colnames(dmrp) <- gsub("CD34-erythroblasts", "Erythroid-precursors", colnames(dmrp))
colnames(dmrp) <- gsub("CD34-megakaryocytes", "Megakaryocytes", colnames(dmrp))

select_tissues <- c(
  "Brain","Breast","Heart","Kidney","Liver","Lung","Ovary","Pancreas",
  "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells","CD8-T-cells",
  "NK-cells","B-cells","Neutrophils","Eosinophils","Monocytes",
  "Erythroid-precursors","Megakaryocytes"
)

# Binarize wilcoxon p-values: <= 0.05 -> 1 else 0
dmrp <- dmrp %>%
  mutate(
    across(
      .cols = matches("_(greater|less)$"),
      .fns  = ~ ifelse(. <= 0.05, 1, 0),
      .names = "{.col}_bin"
    )
  )

  

# --- Build per-library bin-sums ---
hmC_hyper_dmr_binsums <- hmC_hyper_dmr %>%
  inner_join(dmrp %>% select(chr, start, end, contains("bin")),
             by = c("chr","start","end")) %>%
  group_by(selected_tissue) %>%
  summarise(across(contains("bin"), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

umC_hyper_dmr_binsums <- umC_hyper_dmr %>%
  inner_join(dmrp %>% select(chr, start, end, contains("bin")),
             by = c("chr","start","end")) %>%
  group_by(selected_tissue) %>%
  summarise(across(contains("bin"), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

mC_hypo_dmr_binsums <- mC_hypo_dmr %>%
  inner_join(dmrp %>% select(chr, start, end, contains("bin")),
             by = c("chr","start","end")) %>%
  group_by(selected_tissue) %>%
  summarise(across(contains("bin"), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

# Denominators: number of DMR regions per tissue for each set
denoms <- list(
  hmC_hyper = hmC_hyper_dmr %>% count(selected_tissue, name = "selN"),
  umC_hyper = umC_hyper_dmr %>% count(selected_tissue, name = "selN"),
  mC_hypo   = mC_hypo_dmr   %>% count(selected_tissue, name = "selN")
)

# Bundle what we need for the loop
plot_items <- list(
  list(name = "hmC_hyper", binsums = hmC_hyper_dmr_binsums, denom = denoms$hmC_hyper,
       keep_types = c("hmC_greater","mC_less","umC_greater")),
  list(name = "umC_hyper", binsums = umC_hyper_dmr_binsums, denom = denoms$umC_hyper,
       keep_types = c("umC_greater","mC_less","hmC_greater")),
  list(name = "mC_hypo",   binsums = mC_hypo_dmr_binsums,   denom = denoms$mC_hypo,
       keep_types = c("mC_less","umC_greater","hmC_greater"))
)
sel_res <- data.frame()

pdf("figs/fig3/dmr_ranksum_comb.pdf", height = 6, width = 8)
for (it in plot_items) {
  dmr_binsums <- it$binsums
  denom_tbl   <- it$denom
  keep_types  <- it$keep_types
  
  res <- dmr_binsums %>%
    melt(id.vars = "selected_tissue") %>%
    extract(
      variable,
      into  = c("ranksum_tissue", "mark"),
      regex = "^(.*?)_wilcox_([[:alnum:]]+_[[:alnum:]]+)(?:_bin)?$",
      remove = FALSE
    ) %>%
    filter(
      selected_tissue == ranksum_tissue,
      mark %in% keep_types
    ) %>%
    # attach denominators (number of DMRs for that tissue in this set)
    left_join(denom_tbl, by = "selected_tissue") %>%
    mutate(prop = value / selN) %>%
    filter(selected_tissue %in% select_tissues) %>%
    mutate(selected_tissue = factor(selected_tissue, levels = select_tissues)) %>%
    arrange(selected_tissue) %>%
    mutate(
      type = "solid",
      type = ifelse(ranksum_tissue %in% c(
        "CD4-T-cells","CD8-T-cells","NK-cells","B-cells",
        "Neutrophils","Eosinophils","Monocytes",
        "Erythroid-precursors","Megakaryocytes"
      ), "blood", type),
      type = factor(type, levels = c("blood","solid")) # optional order
    ) %>%
    mutate(
      selected_tissue = factor(selected_tissue, levels=select_tissues)
    ) %>%
    mutate(
      dmr_type=it$name
    )
  
  sel_res <- rbind(sel_res,res)
  p <- ggplot(res, aes(x = selected_tissue, y = prop)) +
    geom_bar(stat = "identity") +
    ggtitle(it$name)+
    facet_grid(mark ~ .) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Proportion significant (value / #DMRs)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
  
  print(p)
  
}
dev.off()


base_plot <- function(plot_df) {
  plot_df %>%
  ggplot(aes(x = mark, y = prop)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(
      aes(color = type),
      position = position_jitter(width = 0.1, height = 0.15, seed = 1),
      size = 0.6, alpha = 0.6
    ) +
    facet_wrap(~ dmr_type, nrow = 1) +               # <- correct formula
    scale_color_manual(values = c(blood = "#6BAED6", solid = "#E34A33")) +  # <- add '+'
    theme_classic() +
    theme(legend.position = "bottom") +
    ylim(0,1)
}
pdf("figs/fig3/dmr_ranksum_comb.boxplot.pdf", height = 3, width = 5)
(sel_res %>% filter(dmr_type == "mC_hypo"  & mark != "mC_less") %>%base_plot |
    sel_res %>% filter(dmr_type == "umC_hyper"  & mark != "umC_greater") %>%base_plot |
    sel_res %>% filter(dmr_type == "hmC_hyper"  & mark != "hmC_greater") %>%base_plot 
) %>%  plot_layout(widths = c(1,1,1)) 
dev.off()


sel_res %>%
  group_by(type,dmr_type,mark) %>%
  summarise(
    across(contains("prop"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    .groups = "drop"
  )

# # A tibble: 18 × 4
# type  dmr_type  mark        mean_prop
# <fct> <chr>     <chr>           <dbl>
#   1 blood hmC_hyper hmC_greater     0.992
# 2 blood hmC_hyper mC_less         0.773
# 3 blood hmC_hyper umC_greater     0.493 #49.3%
# 4 blood mC_hypo   hmC_greater     0.332
# 5 blood mC_hypo   mC_less         1    
# 6 blood mC_hypo   umC_greater     0.999
# 7 blood umC_hyper hmC_greater     0.167
# 8 blood umC_hyper mC_less         1    
# 9 blood umC_hyper umC_greater     1    
# 10 solid hmC_hyper hmC_greater     0.869
# 11 solid hmC_hyper mC_less         0.479
# 12 solid hmC_hyper umC_greater     0.163 #83.7%
# 13 solid mC_hypo   hmC_greater     0.405
# 14 solid mC_hypo   mC_less         0.957
# 15 solid mC_hypo   umC_greater     0.950
# 16 solid umC_hyper hmC_greater     0.236
# 17 solid umC_hyper mC_less         0.947
# 18 solid umC_hyper umC_greater     0.956

#### check DMR length by ranksum test ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/")

dmrp <- fread("dmr_rank/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.wilcox_results.csv")

# Standardize tissue labels in column names
colnames(dmrp) <- gsub("CD34-erythroblasts", "Erythroid-precursors", colnames(dmrp))
colnames(dmrp) <- gsub("CD34-megakaryocytes", "Megakaryocytes", colnames(dmrp))

select_tissues <- c(
  "Brain","Breast","Heart","Kidney","Liver","Lung","Ovary","Pancreas",
  "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells","CD8-T-cells",
  "NK-cells","B-cells","Neutrophils","Eosinophils","Monocytes",
  "Erythroid-precursors","Megakaryocytes"
)

# Binarize wilcoxon p-values: <= 0.05 -> 1 else 0
dmrp <- dmrp %>%
  mutate(
    across(
      .cols = matches("_(greater|less)$"),
      .fns  = ~ ifelse(. <= 0.05, 1, 0),
      .names = "{.col}_bin"
    )
  )
dmr_sta <- dmrp %>%
  select(contains("bin")) %>%
  apply(2,sum) %>%
  as.data.frame()
dmr_sta$id <- rownames(dmr_sta)
dmr_sta <- dmr_sta %>%
  separate(
    id,
    into = c("tissue", "test", "mark", "direction", "bin"),
    sep = "_",
    remove = TRUE
  ) 

dmr_sta %>%
  filter(tissue %in% select_tissues) %>%
  mutate(tissue = factor(tissue, levels = select_tissues)) %>%
  ggplot(aes(x = tissue, y = ., fill = direction)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ mark, nrow = 3, scales = "free_x") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(1, "lines")
  )

library(GenomicRanges) 

select_tissue <- "Liver"
threshold <- 0.05

col_hmC <- paste0(select_tissue, "_wilcox_hmC_greater")
col_umC <- paste0(select_tissue, "_wilcox_umC_greater")
print(col_hmC)
print(col_umC)
get_length_freq <- function(dmrp, col_name, label) {
  
  dmrp_sel <- dmrp %>%
    select(chr,start,end,col_name) %>%
    filter(.data[[col_name]] < threshold)
  n_dmr <- nrow(dmrp_sel)
  gr <- makeGRangesFromDataFrame(
    dmrp_sel,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = FALSE
  )
  
  merged_gr <- reduce(gr)
  
  ggplot(data.frame(length = width(merged_gr)),
         aes(x = "", y = length)) +
    geom_violin() +
    scale_y_continuous(trans = "log10") +   # highly recommended
    labs(y = "Merged interval length (bp)", x = col_name) +
    theme_minimal() +
    ylim(0,100000)
}

p1 <- get_length_freq(dmrp, col_hmC, "hmC_greater")
p2 <- get_length_freq(dmrp, col_umC, "umC_greater")

cowplot::plot_grid(p1,p2)










#### DMR heatmap for genebody ####
setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls")
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Monocytes", "Eosinophils", "Erythroid-precursors", "Megakaryocytes")


dmr_files <- c(
  paste0("all_sample.merged.mlml.mincov10_common.groupby.MANE.GRCh38.v1.0.refseq_genomic.gene.hmC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.MANE.GRCh38.v1.0.refseq_genomic.gene.mC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.MANE.GRCh38.v1.0.refseq_genomic.gene.umC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.MANE.GRCh38.v1.0.refseq_genomic.gene.hmC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.MANE.GRCh38.v1.0.refseq_genomic.gene.mC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.MANE.GRCh38.v1.0.refseq_genomic.gene.umC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt")
)
selN <- 20
for(dmr_file in dmr_files){
  plot_heatmap(dmr_file = paste0("dmr_call_genebody/",dmr_file),
               all_meth_file = "dmr_call_genebody/all_sample.merged.mlml.mincov10_common.groupby.MANE.GRCh38.v1.0.refseq_genomic.gene.bed",
               tissue_order=tissue_order,
               selN=selN, width = 10, height = 7, dmr_type = "healthy",
               prefix = paste0("dmr_call_genebody/",gsub("txt$",selN, dmr_file))
               
  )
}