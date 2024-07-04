library(data.table)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(dendextend)
library(ggplot2)
library(reshape2)
library(tidyr)
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/plots/poster")
source("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/code_tissue_methylome_atlas/rename_sample.r")

#### overall methylation level ####
dat <- read.table("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/meth/sample_merged_common_cpg.mincov3.all_meth.sta",col.names = c("smp","meth"))
dat[grep("umC",dat$smp),2] <- 100- dat[grep("umC",dat$smp),2]
dat$tissue <- lapply(dat$smp, function(x){unlist(strsplit(x,"_"))[2]}) %>% unlist()
dat$tissue <- gsub("CD34-erythroblasts","Erythroid-precursors",dat$tissue)
dat$tissue <- gsub("CD34-megakaryocytes","Megakaryocytes",dat$tissue)
dat$tumour <-"no"; dat$tumour[grep("Tumour",dat$tissue)] <- "yes";dat$tumour[grep("-Pancreatitis|-Cirrhosis",dat$smp)] <- "pre"
dat$tissue <- gsub("-Tumour|-Pancreatitis|-Cirrhosis","",dat$tissue)
dat$type <- lapply(dat$smp, function(x){unlist(strsplit(x,"_"))[3]}) %>% unlist()
dat$type[grep("umC",dat$type)] <- "mC+hmC"
dat$type <- factor(dat$type,levels=c("mC","hmC","mC+hmC"));dat<- dat[order(dat$type),]
# tissue_order <- c( "Heart", "Brain", "Brain-Tumour","Breast","Breast-Tumour","Esophagus", "Colon","Colon-Tumour",
#                    "Kidney", "Kidney-Tumour","Liver","Liver-Tumour","Liver-Cirrhosis", "Lung","Lung-Tumour","Ovary", "Ovary-Tumour",
#                    "Pancreas", "Pancreas-Tumour","Pancreas-Pancreatitis","Prostate","Prostate-Tumour","Stomach","Stomach-Tumour",
#                    "Spleen", "CD4-T-cells", "CD8-T-cells", "Neutrophils", "NK-cells", "B-cells", "Eosinophils", "Monocytes", 
#                    "CD34-erythroid-precursors", "CD34-megakaryocytes")
tissue_order <- c( "Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                   "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
dat$tissue <- factor(dat$tissue,levels=tissue_order)
dat <- dat[order(dat$tissue),]
p <- ggplot(dat,aes(x=tissue,y=meth, group=tumour, color=tumour)) + 
  geom_point(aes(shape=tumour), position=position_dodge(width=0.85)) +
  facet_grid(type~.,scale="free_y") +
  scale_shape_manual(values=c(16, 15, 8)) +
  scale_color_manual(values=c("#473183","#829fd5","#e8908f"))+
  theme_bw()+
  theme(legend.position = NULL)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color="black")) +
  ylab("genome methylation %") + xlab("")+
  stat_summary(fun = "mean", 
               geom = "errorbar", color="black",
               aes(ymax = ..y.., ymin = ..y..), 
               position = position_dodge(width = 0.85), 
               width = 0.5)

ggsave("sample_merged_common_cpg.mincov3.all_meth.sta.pdf",p,width = 10,height = 6 )

p <- ggplot(dat,aes(x=tissue,y=meth, group=tumour, color=tumour)) + 
  geom_point(aes(shape=tumour), position=position_dodge(width=0.85)) +
  facet_grid(type~.,scale="free_y") +
  scale_shape_manual(values=c(16, 15, 8)) +
  scale_color_manual(values=c("#473183","#829fd5","#e8908f"))+
  theme_bw()+
  theme(legend.position = NULL)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color="black")) +
  ylab("genome methylation %") + xlab("")+
  stat_summary(fun = "mean", colour = "black", size = 1.8,
               geom = "text", aes(label = round(after_stat(y), 1)),
               position = position_dodge(width=0.85), vjust = 1.5) +
  stat_summary(fun = "mean", 
               geom = "errorbar", color="black",
               aes(ymax = ..y.., ymin = ..y..), 
               position = position_dodge(width = 0.85), 
               width = 0.5)

ggsave("sample_merged_common_cpg.mincov3.all_meth.sta.num.pdf",p,width = 10,height = 6 )


plot_dend <- function(infile,normal_only=FALSE,use_cor=FALSE){
  dat <- fread(infile) %>% as.data.frame()
  colnames(dat) <- gsub("-","_",colnames(dat))
  dat <- rename_columns(dat)
  if(normal_only==TRUE){
    dat <- dat[,grep("Tumour|Cirrhosis|Pancreatitis",colnames(dat),invert = TRUE)]
    outfix <- gsub(".csv",".healthy",paste0(basename(infile),".pdf"))
  }else{
    outfix <- gsub(".csv","",paste0(basename(infile),".pdf"))
  }
  colnames(dat) <- gsub("CD34-erythroblasts","Erythroid-precursors",colnames(dat))
  colnames(dat) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(dat))
  rownames(dat) <- dat[,1]; dat[,1] <- NULL
  if(use_cor==TRUE){
    dend <- dat %>% cor() %>%
      dist %>% hclust %>% as.dendrogram
    outfix <- gsub(".pdf",".cor.pdf",outfix)
  }else{
    dend <- dat %>% t %>%
      dist %>% hclust %>% as.dendrogram
  }
    
  num_colors <- unique(gsub("[0-9]$","",labels(dend))) %>% length()
  color_scheme <- data.frame(smp=labels(dend),tissue=gsub("[0-9]$","",labels(dend))) %>%
    plyr::join(data.frame(tissue=unique(gsub("[0-9]$","",colnames(dat)[-1])),color=colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)))
  dend <- color_labels(dend, col = color_scheme$color)

  
  # num_colors <- unique(gsub(".*_","",labels(dend))) %>% length()
  # color_scheme <- data.frame(smp=labels(dend),tissue=gsub(".*_","",labels(dend))) %>%
  #   plyr::join(data.frame(tissue=unique(gsub(".*_","",colnames(dat)[-1])),color=colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)))
  # dend <- color_labels(dend, col = color_scheme$color)
  # labels(dend) <- gsub(".*_","",labels(dend))
  
  pdf(outfix,width = 10, height = 2)
  par(mar=c(10, 1, 1, 2) + 0.1, cex=0.6)
  plot(dend, edge.root=TRUE, horiz=FALSE, axes=FALSE)
  dev.off()
}

plot_dend("/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/CAPS_samples.all_samples.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv")
plot_dend("/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/TAPSbeta_samples.all_samples.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv")
plot_dend("/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/CAPS_samples.all_samples.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv", normal_only = TRUE)
plot_dend("/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/TAPSbeta_samples.all_samples.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv", normal_only = TRUE)
plot_dend("/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/CAPS_samples.all_samples.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv", normal_only = TRUE,use_cor=TRUE)
plot_dend("/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/TAPSbeta_samples.all_samples.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv", normal_only = TRUE,use_cor=TRUE)


#### Umap & tSNE ####
library(ggrepel)
library(cowplot)
library(Rtsne)
library(uwot)
options(ggrepel.max.overlaps = Inf)
plot_umap <- function(dat, smp_status, outputlabel,  nneighbors=seq(2,10,1), width=10, height=10, min_dist=0.01, learning_rate=0.05, spread=0.1, scale="none"){
  umap_plot <- data.frame()
  for(n in nneighbors){
    umap_out <- umap(as.matrix(t(dat)), n_neighbors = n, init = "pca", min_dist=min_dist, learning_rate=learning_rate)
    umap_plot <- rbind(umap_plot,
                       data.frame(x = umap_out[,1], 
                                  y = umap_out[,2], 
                                  col = smp_status,smp=colnames(dat),
                                  n_neighbour=n))
  }
  p1 <- ggplot(umap_plot,aes(x,y,color=col, label = smp, shape=col)) + 
    geom_point() + 
    geom_text_repel(data=umap_plot, aes(label=smp), size=2) +
    scale_shape_manual(values = rep(seq(1:14),3)[1:length(unique(smp_status))]) +
    facet_wrap(~ n_neighbour, scales = "free", ncol = 3) +
    ggtitle(outputlabel) + 
    theme_classic() +
    theme(plot.title = element_text(size = 8),
          legend.position = "bottom") 
  
  p2 <- ggplot(umap_plot,aes(x,y,color=col, label = smp, shape=col)) +
    geom_point() +
    scale_shape_manual(values = rep(seq(1:14),3)[1:length(unique(smp_status))]) +
    facet_wrap(~ n_neighbour, scales = "free", ncol = 3) +
    ggtitle(outputlabel) +
    theme_classic() +
    theme(plot.title = element_text(size = 8),
          legend.position = "bottom")
  p3 <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
  ggsave(paste0(outputlabel, ".pdf"), p3, width=width*2, height = height)
}




plot_tsne <- function(dat, smp_status,outputlabel, perplexity=seq(2,10,1), width=10, height=10,min_dist=0.01, learning_rate=0.05, spread=0.1, scale="none"){
  tsne_plot <- tibble()
  for(i in perplexity){
    tsne_out <- Rtsne(as.matrix(t(dat)), perplexity = i) # Run TSNE
    tsne_plot <- rbind(tsne_plot,
                       data.frame(x = tsne_out$Y[,1], 
                                  y = tsne_out$Y[,2], 
                                  col = smp_status,
                                  smp=colnames(dat),
                                  perplexity=i)
    )
  }
  
  p1 <- ggplot(tsne_plot,aes(x,y,color=col, label = smp, shape=col)) + 
    geom_point() + 
    geom_text_repel(data=tsne_plot, aes(label=smp), size=2) +
    scale_shape_manual(values = rep(seq(1:14),3)[1:length(unique(smp_status))]) +
    facet_wrap(~ perplexity, scales = "free", ncol = 3) +
    ggtitle(outputlabel) + 
    theme_classic() +
    theme(plot.title = element_text(size = 8),
          legend.position = "bottom") 
  
  p2 <- ggplot(tsne_plot,aes(x,y,color=col, label = smp, shape=col)) +
    geom_point() +
    scale_shape_manual(values = rep(seq(1:14),3)[1:length(unique(smp_status))]) +
    facet_wrap(~ perplexity, scales = "free", ncol = 3) +
    ggtitle(outputlabel) +
    theme_classic() +
    theme(plot.title = element_text(size = 8),
          legend.position = "bottom")
  p3 <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
  ggsave(paste0(outputlabel, ".pdf"), p3, width=width*2, height = height)
}
tissue_order <- c("Brain-Tumour","Breast-Tumour", "Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Colon-Tumour","Stomach-Tumour",
                  "Brain","Breast", "Heart",
                  "Kidney","Liver", "Lung","Ovary", 
                  "Pancreas", "Prostate","Colon","Stomach","Esophagus",
                  "Spleen", "CD4-T-cells", "CD8-T-cells", "Neutrophils", "NK-cells", "B-cells", "Eosinophils", "Monocytes", 
                  "Erythroid-precursors", "Megakaryocytes","Liver-Cirrhosis","Pancreas-Pancreatitis")
caps_file <-"/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/CAPS_samples.all_samples.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv"
caps <- fread(caps_file)
colnames(caps) <- gsub("CD34-erythroblasts","Erythroid-precursors",colnames(caps)); colnames(caps) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(caps))

smp_status <- factor(lapply(colnames(caps)[-1],function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist(),levels=tissue_order)
plot_umap(dat=caps[,-1], smp_status = smp_status,outputlabel = paste0(basename(caps_file),"umap"))
plot_tsne(dat=caps[,-1], smp_status = smp_status,outputlabel = paste0(basename(caps_file),"tsne"))

selcaps <- caps %>% select(-contains("-Tumour"), -contains("-Pancreatitis"), -contains("-Cirrhosis"))
smp_status <- factor(lapply(colnames(selcaps)[-1],function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist(),levels=grep("-Tumour|-Pancreatitis|-Cirrhosis",tissue_order,invert=TRUE,value=TRUE))
plot_umap(dat=selcaps[,-1], smp_status = smp_status,outputlabel = paste0(basename(caps_file),"umap.healthy"))
plot_tsne(dat=selcaps[,-1], smp_status = smp_status,outputlabel = paste0(basename(caps_file),"tsne.healthy"))



tapsb_file <-"/well/ludwig/users/dyp502/tissue_atlas_v3/dmr/variance/TAPSbeta_samples.all_samples.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Top10_variance_blocks.csv"
tapsb <- fread(tapsb_file)
colnames(tapsb) <- gsub("CD34-erythroblasts","Erythroid-precursors",colnames(tapsb)); colnames(tapsb) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(tapsb))

smp_status <- factor(lapply(colnames(tapsb)[-1],function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist(),levels=tissue_order)
plot_umap(dat=tapsb[,-1], smp_status = smp_status,outputlabel = paste0(basename(tapsb_file),"umap"))
plot_tsne(dat=tapsb[,-1], smp_status = smp_status,outputlabel = paste0(basename(tapsb_file),"tsne"))

seltapsb <- tapsb %>% select(-contains("-Tumour"), -contains("-Pancreatitis"), -contains("-Cirrhosis"))
smp_status <- factor(lapply(colnames(seltapsb)[-1],function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist(),levels=grep("-Tumour|-Pancreatitis|-Cirrhosis",tissue_order,invert=TRUE,value=TRUE))
plot_umap(dat=seltapsb[,-1], smp_status = smp_status,outputlabel = paste0(basename(tapsb_file),"umap.healthy"))
plot_tsne(dat=seltapsb[,-1], smp_status = smp_status,outputlabel = paste0(basename(tapsb_file),"tsne.healthy"))

