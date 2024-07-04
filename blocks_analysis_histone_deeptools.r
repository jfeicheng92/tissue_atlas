### deeptools_plot
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr15/histone_normal/")
options(bitmapType='cairo-png')
.libPaths(c(.libPaths(),
            "/gpfs3/users/ludwig/cfo155/R/x86_64-pc-linux-gnu-library/4.2",
            "/gpfs3/users/ludwig/cfo155/R/x86_64-pc-linux-gnu-library/4.3"))

library(rGREAT)
library(GenomicRanges)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(preprocessCore)
library(cowplot)
library(ggrepel)
library(forecast)
options(ggrepel.max.overlaps = Inf)
tissue_pairs <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/histone/tissue_pairs.txt", header=FALSE)
tissue_pairs <- tissue_pairs[!is.na(tissue_pairs$V2),]


plot_dmb_histone <- function(rname, fname, oname, dmb, tissue_pairs, seltissue="None"){
  mat <- fread(fname,skip = 1)
  tmp <-  readLines(fname, n = 1) %>% strsplit('],"') %>% unlist()
  
  params <- data.frame(par=lapply(gsub('\\{|\\}|\\"|\\[|\\]',"", tmp), function(x)unlist(strsplit(x,":"))[[1]]) %>% unlist,
                       value=lapply(gsub('\\{|\\}|\\"|\\[|\\]',"", tmp), function(x)unlist(strsplit(x,":"))[[2]]) %>% unlist)
  n <- params[params$par=='sample_boundaries',]$value %>% strsplit(",",3) %>% unlist %>% as.numeric() %>% `[`(2)
  samples <- gsub(",_","_",params[params$par=='sample_labels',]$value) %>% strsplit(",",3) %>% unlist %>% as.character()
  region <- fread(rname)
  region$region <- paste0(region$V4,"_",region$V5)
  region$V5 <- paste0(region$V1,":",region$V2,"-",region$V3)
  
  mat <- merge(region%>%select(region,V5),mat,by.x=c("V5"),by.y="V4")
  mat_agg <- data.frame(matrix(ncol = 6, nrow =0))
  pdf(oname, width=24, height = 4)
  for(i in tissue_pairs$V1){
    mat_mean <- mat[,-c(1,3:7)] %>%
      filter(region==paste0(dmb,"_",i)) %>%
      summarize(across(starts_with("V"), mean, na.rm = TRUE))
    colnames(mat_mean) <-  c(paste0(rep(samples,each=n),":",rep(seq(1,n),length(samples))))
    mat_mean <- pivot_longer(mat_mean, cols = contains(":"), names_to = "idx", values_to = "level")
    
    mat_mean$pos <- as.numeric(gsub(".*:","",mat_mean$idx))
    mat_mean$label <- gsub(":.*","",mat_mean$idx)
    mat_mean$cell <- gsub(".*\\.","",mat_mean$label)
    if(all(seltissue == "None")){
      mat_mean <- mat_mean
    }else{
      mat_mean <- mat_mean[grepl(paste(tissue_pairs$V2[tissue_pairs$V1%in%unique(c(seltissue,i))], collapse = "|"), mat_mean$cell),]
    }
    
    mat_mean_sel <- mat_mean %>%
      select(level,pos,cell) %>%
      group_by(cell, pos) %>%
      summarize(mean_value = mean(level))
    mat_mean_sel$target <- "2"
    mat_mean_sel$target[grep(tissue_pairs$V2[tissue_pairs$V1==i],mat_mean_sel$cell)] <- "1"
    
    anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(mat_mean_sel$cell))))
    p1 <- ggplot(mat_mean_sel, aes(x=pos,y=mean_value,color=cell)) + 
      geom_line(aes(linetype=target)) +
      theme_bw() + 
      ggtitle(paste0(i," ",gsub(".gz","",fname))) +
      scale_x_continuous(breaks = c(0,150,300),labels = c("-15K","0","15K"))+
      scale_color_manual(values=anno_colors) + ylab(paste0(gsub(".*DMB.|.gz","",fname),"\n raw"))
    mat_mean_sel$seltissue <- i
    mat_mean_sel$value <- "raw"
    ## smoothing average
    ma_order<-10
    mat_mean_sel_avg <- mat_mean_sel %>% 
      select(cell,pos,mean_value) %>%
      pivot_wider(names_from = cell, values_from = mean_value) %>%
      apply(2,function(x)c(x[1:(ma_order/2)],
                           ma(x,order=ma_order)[-c(1:(ma_order/2),(length(x)-(ma_order/2)):length(x))],
                           x[(length(x)-(ma_order/2)):length(x)])) %>%
      as.data.frame() %>%
      pivot_longer(cols = -pos, names_to = "cell", values_to = "mean_value")
    mat_mean_sel_avg$target <- "2"
    mat_mean_sel_avg$target[grep(tissue_pairs$V2[tissue_pairs$V1==i],mat_mean_sel_avg$cell)] <- "1"
    p2 <- ggplot(mat_mean_sel_avg, aes(x=pos,y=mean_value,color=cell)) + 
      geom_line(aes(linetype=target)) +
      theme_bw() + 
      ggtitle(paste0(i," ",gsub(".gz","",fname))) +
      scale_x_continuous(breaks = c(0,150,300),labels = c("-15K","0","15K"))+
      scale_color_manual(values=anno_colors) + ylab(paste0(gsub(".*DMB.|.gz","",fname),"\n smooth average"))
    mat_mean_sel_avg$seltissue <- i
    mat_mean_sel_avg$value <- "average"
    ## smoothing average
    mat_mean_sel_avg_nor <- mat_mean_sel %>% 
      select(cell,pos,mean_value) %>%
      pivot_wider(names_from = cell, values_from = mean_value)
    nor_factor <- apply(mat_mean_sel_avg_nor[,-1],2,median)
    for (j in names(nor_factor)) {
      mat_mean_sel_avg_nor <- mat_mean_sel_avg_nor %>%
        mutate(!!j := .data[[j]] / nor_factor[[j]])
    }
    
    mat_mean_sel_avg_nor <- apply(mat_mean_sel_avg_nor, 2,function(x)c(x[1:(ma_order/2)],
                                                                       ma(x,order=ma_order)[-c(1:(ma_order/2),(length(x)-(ma_order/2)):length(x))],
                                                                       x[(length(x)-(ma_order/2)):length(x)])) %>%
      as.data.frame() %>%
      pivot_longer(cols = -pos, names_to = "cell", values_to = "mean_value")
    mat_mean_sel_avg_nor$target <- "2"
    mat_mean_sel_avg_nor$target[grep(tissue_pairs$V2[tissue_pairs$V1==i],mat_mean_sel_avg_nor$cell)] <- "1"
    p3 <- ggplot(mat_mean_sel_avg_nor, aes(x=pos,y=mean_value,color=cell)) + 
      geom_line(aes(linetype=target)) +
      theme_bw() + 
      ggtitle(paste0(i," ",gsub(".gz","",fname))) +
      scale_x_continuous(breaks = c(0,150,300),labels = c("-15K","0","15K"))+
      scale_color_manual(values=anno_colors) + ylab(paste0(gsub(".*DMB.|.gz","",fname),"\n normalized by median"))
    mat_mean_sel_avg_nor$seltissue <- i
    mat_mean_sel_avg_nor$value <- "nor"
    print(cowplot::plot_grid(p1,p2,p3,nrow=1,rel_widths = c(1,1,1)))
    mat_agg <- rbind(mat_agg, mat_mean_sel_avg, mat_mean_sel_avg_nor)
    
  }
  dev.off()
  write.table(mat_agg,gsub(".pdf",".agg.txt",oname),sep="\t",quote = FALSE, row.names = FALSE)
  p4 <- mat_agg %>%
    filter(target==1) %>%
    filter(value=="average") %>%
    ggplot(aes(x=pos,y=mean_value,color=cell)) + 
    geom_line(aes(linetype=target)) +
    theme_bw() + 
    ggtitle(paste0(i," ",gsub(".gz","",fname))) +
    scale_x_continuous(breaks = c(0,150,300),labels = c("-15K","0","15K"))+
    ylab(dmb_histone)
  mat_agg_wide <- mat_agg %>%
    filter(target==1) %>%
    filter(value=="average") %>%
    select(pos,cell,mean_value) %>%
    pivot_wider(names_from = pos,values_from = mean_value)
  
  p5 <- cbind(mat_agg_wide[,1],mat_agg_wide[,-1] %>% apply(., 1, scale) %>% t()) %>%
    melt(id.vars=c("cell")) %>%
    ggplot(aes(x = cell, y = variable)) +
    geom_tile(aes(fill = value)) +
    coord_flip() +
    theme(legend.position = "right") +
    xlab("Tissue") + ylab("") +
    theme(axis.text=element_text(size=10)) +
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)) +
    scale_y_discrete(breaks=c(1,150,300),labels = c("-15K","0","15K")) +
    theme_minimal()
  pdf(gsub(".pdf",".agg.pdf",oname), width=20, height = 4)
  print(cowplot::plot_grid(p4,p5,nrow=1,rel_widths = c(1,1,1)))
  dev.off()
}

for(dmb_histone in list.files(pattern = "TAPSbeta_Hypo_DMB.histone_.*.gz")){
  plot_dmb_histone(rname="TAPSbeta_Hypo_DMB.bed",fname=dmb_histone,
                   oname=gsub(".gz",".pdf",dmb_histone),dmb="Hypo",tissue_pairs=tissue_pairs)
  
    
}


for(dmb_histone in list.files(pattern = "TAPSbeta_Hypo_DMB.histone_.*.gz")){
  plot_dmb_histone(rname="TAPSbeta_Hypo_DMB.bed",fname=dmb_histone,
                   oname=gsub(".gz",".pdf",dmb_histone),dmb="Hypo",tissue_pairs=tissue_pairs)
  
  
}



tissue_pairs <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/histone/tissue_pairs.txt", header=FALSE, col.names = c("tissue","cell"))
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
tissue_pairs$tissue <- factor(tissue_pairs$tissue, levels=tissue_order)
tissue_pairs <- tissue_pairs[order(tissue_pairs$tissue),]


for(fname in list.files(pattern=".*agg.txt")){
  mat_agg <- fread(fname)
  anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(mat_agg$cell))))
  cell_order <- sapply(tissue_pairs$cell , function(p) {grep(patt=p, unique(mat_agg$cell))} ) %>% unlist() %>% unique(mat_agg$cell)[.] 
  
  p4 <- mat_agg %>%
    filter(target==1) %>%
    filter(value=="average") %>%
    ggplot(aes(x=pos,y=mean_value,color=cell)) + 
    geom_line() +
    theme_bw() + 
    ggtitle(gsub(".txt","",fname)) +
    scale_x_continuous(breaks = c(0,150,300),labels = c("-15K","0","15K")) +
    ylab(gsub(".txt","",fname)) +
    scale_color_manual(values=anno_colors)
  mat_agg_wide <- mat_agg %>%
    filter(target==1) %>%
    filter(value=="average") %>%
    select(pos,cell,mean_value) %>%
    pivot_wider(names_from = pos,values_from = mean_value)
  mat_agg_wide$cell <- factor(mat_agg_wide$cell,rev(cell_order))
  mat_agg_wide <- mat_agg_wide[order(mat_agg_wide$cell),]
  
  
  p5 <- cbind(mat_agg_wide[,1],mat_agg_wide[,-1] %>% apply(., 1, scale) %>% t()) %>%
    melt(id.vars=c("cell")) %>%
    ggplot(aes(x = cell, y = variable)) +
    geom_tile(aes(fill = value)) +
    coord_flip() +
    theme(legend.position = "right") +
    xlab("Tissue") + ylab("") +
    theme(axis.text=element_text(size=10)) +
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)) +
    scale_y_discrete(breaks=c(1,150,300),labels = c("-15K","0","15K")) +
    theme_minimal()
  
  pdf(gsub(".txt",".reorder.pdf",fname), width=12, height = 4)
  print(cowplot::plot_grid(p4,p5,nrow=1,rel_widths = c(1,0.6)))
  dev.off()
  
}





tissue_pairs <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/histone/tissue_pairs.txt", header=FALSE, col.names = c("tissue","cell"))
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroblasts", "Megakaryocytes")
tissue_pairs$tissue <- factor(tissue_pairs$tissue, levels=tissue_order)
tissue_pairs <- tissue_pairs[order(tissue_pairs$tissue),]


for(fname in list.files(pattern=".*agg.txt")){
  mat_agg <- fread(fname)
  anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(mat_agg$cell))))
  cell_order <- sapply(tissue_pairs$cell , function(p) {grep(patt=p, unique(mat_agg$cell))} ) %>% unlist() %>% unique(mat_agg$cell)[.]
  mat_agg$cell <- factor(mat_agg$cell, levels=cell_order)
 p <- mat_agg%>%
   filter(!seltissue %in% c("Erythroblasts","Megakaryocytes")) %>%
   filter(value=="average") %>%
   ggplot(aes(x=pos,y=mean_value,color=cell)) + 
   geom_line(aes(linetype=as.factor(target))) +
   theme_classic() + 
   facet_wrap( ~ seltissue, nrow = 5, scales="free")+
   scale_x_continuous(breaks = c(0,150,300),labels = c("-15K","0","15K"))+
   scale_color_manual(values=anno_colors)  +
   theme(legend.position = "bottom")
 ggsave(gsub(".txt",".all_tissue.pdf",fname),p,width = 8, height = 10)

}


