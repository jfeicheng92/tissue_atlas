#### Umap & tSNE ####
library(ggrepel)
library(cowplot)
library(Rtsne)
library(uwot)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(cowplot)
library(reshape2)
library(eulerr)
library(stringr)
library(matrixStats)
library(dendextend)
options(bitmapType='cairo-png')
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/")
source("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/code_tissue_methylome_atlas/rename_sample.r")
#### Load methylation data ####
prefix <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500"
meth <- fread(paste0("dmr_call_grail/",prefix, ".bed"), header=TRUE) %>% as.data.frame()


for(mod in c("_mC","_umC","_hmC")){
  # matrix of the mC columns (samples)
  mC_mat <- meth %>% select(ends_with(mod)) %>% as.matrix()
  
  # row-wise variance across samples (needs ≥2 non-NA per row)
  valid_rows <- rowSums(!is.na(mC_mat)) >= 2
  mC_var <- rep(NA_real_, nrow(mC_mat))
  mC_var[valid_rows] <- rowVars(mC_mat[valid_rows, , drop = FALSE], na.rm = TRUE)
  
  # add to data
  meth$mC_var <- mC_var
  
  # choose top 5‰ (0.5%) most variable bins
  k <- max(1, ceiling(nrow(meth) * 5/1000))              # count-based selection
  idx_top <- order(meth$mC_var, decreasing = TRUE)[seq_len(k)]
  meth_top <- meth[idx_top, ]
  
  # (alternative) quantile-based threshold (may return a few more if ties)
  thr <- quantile(meth$mC_var, probs = 1 - 5/100, na.rm = TRUE)
  meth_top_q <- meth %>% filter(mC_var >= thr) %>% select(chr,start,end,ends_with(mod))
  write.csv(meth_top_q, paste0("figs/fig2/",prefix,mod,"top_var.txt"),quote = FALSE, row.names = FALSE)
}







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
mod_files <-c("figs/fig2/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500_hmCtop_var.txt",
              "figs/fig2/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500_umCtop_var.txt",
              "figs/fig2/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500_mCtop_var.txt")
for(mod_file in mod_files){
  mod <- fread(mod_file)
  colnames(mod) <- gsub("CD34-erythroblasts","Erythroid-precursors",colnames(mod)); 
  colnames(mod) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(mod))
  
  smp_status <- factor(lapply(colnames(mod)[-c(1:3)],function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist(),levels=tissue_order)
  plot_umap(dat=mod[,-c(1:3)], smp_status = smp_status,outputlabel = paste0(mod_file,"umap"))
  plot_tsne(dat=mod[,-c(1:3)], smp_status = smp_status,outputlabel = paste0(mod_file,"tsne"))
  
  sel_mod <- mod %>% select(-contains("-Tumour"), -contains("-Pancreatitis"), -contains("-Cirrhosis"))
  smp_status <- factor(lapply(colnames(sel_mod)[-c(1:3)],function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist(),levels=grep("-Tumour|-Pancreatitis|-Cirrhosis",tissue_order,invert=TRUE,value=TRUE))
  plot_umap(dat=sel_mod[,-c(1:3)], smp_status = smp_status,outputlabel = paste0(mod_file,"umap.healthy"))
  plot_tsne(dat=sel_mod[,-c(1:3)], smp_status = smp_status,outputlabel = paste0(mod_file,"tsne.healthy"))
  
}



plot_dend <- function(infile,normal_only=FALSE,use_cor=FALSE){
  dat <- fread(infile) %>% as.data.frame()
  colnames(dat) <- gsub("-","_",colnames(dat)) %>% gsub("_hmC|_mC|_umC","",.)
  dat <- rename_columns(dat)
  if(normal_only==TRUE){
    dat <- dat[,grep("Tumour|Cirrhosis|Pancreatitis",colnames(dat),invert = TRUE)]
    outfix <- gsub(".csv",".healthy",paste0(infile,".pdf"))
  }else{
    outfix <- gsub(".csv","",paste0(infile,".pdf"))
  }
  colnames(dat) <- gsub("CD34-erythroblasts","Erythroid-precursors",colnames(dat))
  colnames(dat) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(dat))
  rownames(dat) <- paste(dat$chr, dat$start, dat$end,sep="_"); dat[,1:3] <- NULL
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
  
  pdf(outfix,width = 10, height = 2)
  par(mar=c(10, 1, 1, 2) + 0.1, cex=0.6)
  plot(dend, edge.root=TRUE, horiz=FALSE, axes=FALSE)
  dev.off()
}

for(mod_file in mod_files){
  plot_dend(mod_file)
  plot_dend(mod_file, normal_only = TRUE)
  plot_dend(mod_file, normal_only = TRUE,use_cor=TRUE)
 }

