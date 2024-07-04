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
library(ggseqlogo)
#### convert motif to meme format ####
setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/homer/motifs")
for(i in list.files(pattern=c(".motif"))){
  motif <- fread(i,skip = 1) %>% t() %>% as.data.frame()
  #(A/C/G/T)
  colnames(motif) <- gsub("V","",colnames(motif))
  motif <- cbind(seq=c('"A"','"C"','"G"','"T"'),motif)
  colnames(motif)[1] <- paste0('"',gsub(".motif","",i),'"')
  write.table(motif,gsub(".motif",".taipale.fmt",i),quote = FALSE,sep="\t",row.names =FALSE)
}

#### plot motif by cluster ####
setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/homer/motifs")
motif_cluster <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/homer/motifs/all_motif_tomtom.all.cluster.name.txt",col.names=c("id","motif","cluster"))
pfms <- list()
for(i in motif_cluster$id){
  motif <- fread(paste0(i,".motif"),skip = 1) %>% t() %>% as.matrix()
  #(A/C/G/T)
  rownames(motif) <- c("A","C","G","T")
  pfms[[paste0(i, ":cluster ", motif_cluster$cluster[motif_cluster$id==i])]] <- motif
}
# pdf("motif_logo.cluster.pdf",width = 30, height = 30)
# print(ggseqlogo(pfms,nrow=30))
# dev.off()


pfms <- list()
for(i in motif_cluster$id){
  motif <- fread(paste0(i,".motif"),skip = 1) %>% t() %>% as.matrix()
  #(A/C/G/T)
  rownames(motif) <- c("A","C","G","T")
  pfms[[i]] <- motif
}


#### Plot enriched motif ####
WORKDIR<-"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr15/motif/"
setwd(WORKDIR)
# allmotif <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr9/motif/all_motif.txt")
allmotif <- fread("all_motif.txt")
motif_cluster <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/homer/motifs/all_motif_tomtom.all.cluster.name.txt",col.names=c("id","motif","cluster"))
allmotif <- merge(allmotif,motif_cluster,by.x=c("Motif Name"),by.y = "motif")

# allmotif %>% filter(`q-value (Benjamini)`<0.001) %>% merge(., motif_cluster,by.x=c("Motif Name"),by.y=c("motif")) %>%
#   write.table(.,"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr9/motif/all_motif.q0.001.txt",sep = "\t",quote = FALSE, row.names = FALSE)



for(exp in c("TAPSbeta_Hypo_","CAPS_Hyper_")){
  motif <- allmotif[grep(exp,allmotif$bed),]
  motif$odds <- gsub("%","",motif$`% of Target Sequences with Motif`) %>% as.numeric()/
    gsub("%","",motif$`% of Background Sequences with Motif`) %>% as.numeric()
  ### choose enriched motif based on q-value
  # enriched_motif <- motif %>% filter(`q-value (Benjamini)`<0.0001) %>% 
  #   select(`Motif Name`) %>% unique() %>% unlist() %>% as.character()
  enriched_motif <- motif %>% filter(`P-value`<0.01) %>% 
    select(`Motif Name`) %>% unique() %>% unlist() %>% as.character()

  
  motif_wide <- motif %>% select(bed,`Motif Name`,`Log P-value`) %>%
    pivot_wider(names_from = `Motif Name`,
                values_from = `Log P-value`,values_fn=mean ) %>%
    select(c("bed",enriched_motif))
  
  ### for each cluster choose the top enriched motif as representatives.
  motif_score <- apply(motif_wide[,-1],2,min) %>% as.data.frame()%>% merge(motif_cluster,by.x=0,by.y="motif") 
  colnames(motif_score)[c(1:2)] <- c("motif",y="score")
  motif_score$logp <- -motif_score$score
  motif_score_cluster <- motif_score %>% group_by(cluster) %>%
    top_n(1, logp)
  sel_motif_wide <- motif_wide %>% select(c("bed",motif_score_cluster$motif )) %>% as.data.frame()
  rownames(sel_motif_wide) <- sel_motif_wide$bed;sel_motif_wide$bed<- NULL
  colnames(sel_motif_wide) <- data.frame(motif=c(colnames(sel_motif_wide))) %>% merge(., motif_cluster) %>% select(id) %>% unlist()
  
  ### reorder tissue and motif
  tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                    "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroblasts", "Megakaryocytes")
  sel_motif_wide <- sel_motif_wide[order(factor(rownames(sel_motif_wide),levels=paste0(exp,tissue_order))),]
  
                                                                      
                                                                                           
                               
  
  if(exp=="TAPSbeta_Hypo_"){
    # motif_order <- colnames(sel_motif_wide)[apply(sel_motif_wide,1,function(x)which(x< -10)) %>% unlist() %>% as.numeric() %>% unique()]
    # motif_order <- c("sox21","ap2gamma","egr1","hoxb13","klf5","p63","otx2","hif1b","fos","mef2b","pax5","tead1","hnf4","nkx2.1","nr5a2","e2a","hnf6","rbpj1","snai2","ar-half","fox-ebox","gata3.ir3","gata3.ir4","hoxc9","nf1-half","gata6","erg","lef1","runx1","stat3.il23","bzip-irf","cebp-ap1","cebp-cebp","cebp","myb","tbr1","ebf2","ets-ebox","oct11","pu1-irf8","gfi1b","gata-scl")
     motif_order <- c("sox3", "ap2gamma", "mef2d",  "tead", "hnf6", "mafF", "sf1", "rbpj1-ebox","hoxc9","ar-half","nf1-half","tata","hoxd13", "p53",
                      "zfx","lef1","znf264","sp1","fosl2","runx2","tbr1","ebf2","ascl1","pax5","e47","oct11","pu1","cebp-ap1","cebp","pu1-irf8","myb","gata3.ir3","gata3.ir4","gata6","gata-scl","irf4","stat5")
  }else{
     motif_order <- colnames(sel_motif_wide)[apply(sel_motif_wide,1,function(x)which(x< -10)) %>% unlist() %>% as.numeric() %>% unique()]
  }
  sel_motif_wide <- sel_motif_wide %>% select(motif_order)

  
  sel_motif_long <- motif %>% 
    filter(motif$`Motif Name` %in% motif_score_cluster$motif)
  sel_motif_long$id <- factor(sel_motif_long$id, levels=motif_order)
  sel_motif_long <- sel_motif_long[order(sel_motif_long$id),]
  sel_motif_long$bed <-factor(sel_motif_long$bed,levels=paste0(exp,rev(tissue_order)))
  sel_motif_long <- sel_motif_long[order(sel_motif_long$bed),]
  
  ### Plots
  outfix <- paste0(WORKDIR,exp)
  pdf(paste0(outfix,"enriched_motif_logP_heatmap.pdf"),height = 6,width = 15)
  breaksList1 <- seq(-20,0,1)
  print(pheatmap(sel_motif_wide,
                 color = rev(colorRampPalette(brewer.pal(n = 7, name = "Blues"))(length(breaksList1))),
                 breaks = breaksList1,
                 cluster_rows = FALSE, cluster_cols = FALSE))
  breaksList2 <- seq(0,8,1)
  print(sel_motif_long %>%
    ggplot(aes(x=bed, y = id, color =odds , size = -`Log P-value`)) + 
    geom_point() + 
    scale_color_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(length(breaksList2)), 
                          limits = c(0,max(breaksList2)), oob = scales::squish, name = 'odds') +
    cowplot::theme_cowplot() + 
    coord_flip()+
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') +
    theme(axis.ticks = element_blank(),
          legend.position = "bottom"))
  dev.off()
  pdf(paste0(outfix,"enriched_motif.pdf"),height = 15,width = 4)
  p1 <- ggseqlogo(pfms[colnames(sel_motif_wide)],ncol=1) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 8, margin = margin()))
  p2 <- ggseqlogo(pfms[colnames(sel_motif_wide)],ncol=1) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 0, margin = margin()))
  print(cowplot::plot_grid(p1,p2,nrow=1))

  dev.off()
}
