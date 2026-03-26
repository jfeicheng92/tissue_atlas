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
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls")

#### dmr quant parameters ####
bg_tg_q <- "bgQ0.05-0.05.tgQ0.25"
#### hyper vs. hypo DMR ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/") ## using quantiles of groups instead of individual samples
tissue_order <- c("Brain-Tumour","Breast-Tumour","Colon-Tumour","Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Stomach-Tumour")


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
ggsave("../figs/fig6/topDMR_hyper_hypo.tumour.pdf",p, width = 6, height = 5)


#### DMR heatmap ####
setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls")
dir.create("figs/fig6/heatmap")
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

dmr_files <- c(
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt")
)
selN <- 200


tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Monocytes", "Eosinophils", "Erythroid-precursors", "Megakaryocytes")
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
                  "Brain-Tumour", "Breast-Tumour",  "Kidney-Tumour", "Liver-Tumour", "Lung-Tumour", "Ovary-Tumour", "Pancreas-Tumour", "Prostate-Tumour","Colon-Tumour","Stomach-Tumour","Liver-Cirrhosis","Pancreas-Pancreatitis")
for(dmr_file in dmr_files){
  # plot_heatmap(dmr_file = paste0("dmr_call_new1/",dmr_file),
  #              all_meth_file = "dmr_call_new1/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.bed",
  #              tissue_order=tissue_order,
  #              selN=selN, width = 12, height = 7, dmr_type = "all",
  #              prefix = paste0("figs/fig6/heatmap/",gsub("txt$",selN, dmr_file))
  #              
  # )
  plot_heatmap(dmr_file = paste0("dmr_call_new1/",dmr_file),
               all_meth_file = "dmr_call_new1/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.bed",
               tissue_order=tissue_order,
               selN=selN, width = 12, height = 7, dmr_type = "tumour",
               prefix = paste0("figs/fig6/heatmap/tumour.",gsub("txt$",selN, dmr_file))
               
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


