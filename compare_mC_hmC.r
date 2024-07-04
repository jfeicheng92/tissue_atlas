library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(stringr)
library(tidyr)
.libPaths(c(.libPaths(),"/gpfs3/users/ludwig/cfo155/R/x86_64-pc-linux-gnu-library/4.2"))
# setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr9/genelist")
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr15/genelist")
result_list <- data.frame()
caps_pattern <- "capsHyper_.*_genelist$"
taps_pattern <- "tapsbetaHypo_.*_genelist$"
caps_files <- grep("Tumour|Pancreas-Pancreatitis|Liver-Cirrhosis",list.files(pattern = caps_pattern),invert = TRUE,value=TRUE)
taps_files <- grep("Tumour|Pancreas-Pancreatitis|Liver-Cirrhosis",list.files(pattern = taps_pattern),invert = TRUE,value=TRUE)

tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroblasts", "Megakaryocytes")


idx <- order(factor(lapply(caps_files, function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist,levels=tissue_order))

for (i in caps_files[idx]) {
  for (j in taps_files[idx]) {
    # Read data from files
    x1 <- read.table(i, header = FALSE)
    x2 <- read.table(j, header = FALSE)
    
    # Calculate overlap and lengths
    overlap <- length(intersect(x1$V1, x2$V1))
    nx1 <- length(x1$V1)
    nx2 <- length(x2$V1)
    
    # Create a 2x2 contingency table
    contingency_table <- matrix(c(
      overlap, nx1 - overlap,
      nx2 - overlap, 19062 + overlap - nx1 - nx2
    ), nrow = 2, dimnames = list(
      list1 = c("hmC", "nhmC"),
      list2 = c("mC", "nmC")
    ))
    
    # Perform Fisher's exact test
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    
    # Extract relevant statistics
    p_value <- fisher_result$p.value
    odds_ratio <- fisher_result$estimate %>% round(2)
    
    # Store results in a list
    result_list <- c(
      result_list,
      list(
        gsub("_genelist","",i),
        gsub("_genelist","",j),
        overlap,
        nx1,
        nx2,
        odds_ratio,
        p_value
      )
    )
  }
}

# Combine the results into a data frame
enrich <-  data.frame(matrix(do.call(rbind, result_list)%>%unlist,
                         nrow=length(result_list)/7,byrow=TRUE))

# Print the results
# colnames(result_df) <- c("caps_Hyper","tapsbeta_Hypo","noverlap",
#                          "num_caps_Hyper","num_tapsbeta_Hypo",
#                          "odds_ratio","pvalue")
colnames(enrich)  <- c("list1","list2","overlap","nlist1","nlist2", "odds_ratio","pvalue")

odds <- reshape(enrich %>% select("list1","list2","odds_ratio"), 
        idvar = "list1", timevar = "list2", direction = "wide") %>% as.data.frame()
odds[,-1] <- apply(odds[-1], 2, as.numeric) %>% as.data.frame()
rownames(odds) <- odds[,1]; odds[,1] <- NULL;
colnames(odds) <- gsub("odds_ratio.","",colnames(odds))


pvalue <- reshape(enrich %>% select("list1","list2","pvalue"), 
                idvar = "list1", timevar = "list2", direction = "wide") %>% as.data.frame()
pvalue[,-1] <- apply(pvalue[-1], 2, as.numeric) %>% as.data.frame()
rownames(pvalue) <- pvalue[,1]; pvalue[,1] <- NULL;
colnames(pvalue) <- gsub("pvalue.","",colnames(pvalue))



enrich$anno <- paste0(enrich$overlap,"\n",enrich$nlist1,"\n",enrich$nlist2)
enrich_anno <- reshape(enrich[,c(1,2,8)], idvar = "list1", timevar = "list2", direction = "wide") 
rownames(enrich_anno) <- enrich_anno$list1
colnames(enrich_anno) <- gsub("ratio.","",colnames(enrich_anno))
enrich_anno[,1]<- NULL
enrich_anno_sel <- enrich_anno[grep("capsHyper",rownames(enrich_anno)),grep("tapsbetaHypo",colnames(enrich_anno))] 

pdf("compare_mC_hmC_DMB_gene_overlap.pdf", width = 12, height = 12)
rownames(odds) <- gsub("capsHyper_","Hyper DhMB ",rownames(odds))
colnames(odds) <- gsub("tapsbetaHypo_","Hypo DMB ",colnames(odds))

pheatmap(odds, cluster_rows = FALSE, cluster_cols = FALSE, 
         color =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         breaks = seq(0,10,10/100),display_numbers = enrich_anno_sel, main = "color by odds; display:num_overlap(num_DMB-genelist)")
rownames(pvalue) <- gsub("capsHyper_","Hyper DhMB ",rownames(pvalue))
colnames(pvalue) <- gsub("tapsbetaHypo_","Hypo DMB ",colnames(pvalue))
pheatmap(-log(pvalue,10), cluster_rows = FALSE, cluster_cols = FALSE, 
         color =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         breaks = seq(0,10,10/100),display_numbers = enrich_anno_sel, main = "color by -log10 P-value; display:num_overlap(num_DMB-genelist)")
dev.off()

### compare tss distance
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/blocks_tissue/genelist")
result_list <- data.frame()
caps_pattern <- "caps_Hyper_.*_dis$"
taps_files <- list.files(pattern = taps_pattern)
for (i in caps_files) {
    # Read data from files
    x1 <- read.table(i, header = FALSE);x1$V3 <- c("caps_Hyper")
    x2 <- gsub("caps_Hyper","tapsbeta_Hypo",i) %>% read.table(header = FALSE);x2$V3 <- c("tapsbeta_Hypo")
    summary(dis$V2.x)
    summary(dis$V2.y)
    # Calculate overlap and lengths
    dis <- merge(x1,x2,by="V1")
      
    ggplot(dis,aes(x=V3,y=V4)) + geom_bar()
}


### compare mC and hmC in DMB
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/mC_vs_hmC")
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "CD34-erythroblasts", "CD34-megakaryocytes")
tumour_order <- c("Brain-Tumour", "Breast-Tumour", "Colon-Tumour", "Kidney-Tumour", "Liver-Tumour", "Lung-Tumour", "Ovary-Tumour", "Pancreas-Tumour", "Prostate-Tumour","Stomach-Tumour")
tissue_order <- c(tissue_order, tumour_order)


source("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/code_tissue_methylome_atlas/rename_sample.r")

plot_meth <- function(infile, dmr_bed, dmr, tissue_order, tissue="normal",selN=50) {
  bed <- fread(dmr_bed)
  bed$tissue <- paste(bed$selected_tissue, bed$Hyper_or_Hypo,sep="_")
  bed <- bed %>%
    filter(Hyper_or_Hypo==dmr) %>%
    group_by(tissue) %>%
    slice_head(n = selN) %>% as.data.frame() 

  dat1 <- fread(infile) %>% as.data.frame()
  dat1$chr_start_end <- paste(dat1$chr,dat1$start,dat1$end,sep="_")
  dat1 <- merge(dat1,bed%>%select(chr_start_end),by=c("chr_start_end"))
  dat1$chr_start_end  <- NULL
  
  dat1[,grep("umC",colnames(dat1))] <- 100-dat1[,grep("umC",colnames(dat1))]

  colnames(dat1) <- gsub("umC","amC",colnames(dat1))
  # colnames(dat1) <- gsub("CD34-erythroblasts","Erythroblasts",colnames(dat1)); 
  # colnames(dat1) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(dat1))
  dat1 <- dat1[grep(dmr,dat1$tissue),]
  dat1$tissue <- as.factor(gsub(paste0("_",dmr),"",dat1$tissue))
  
  if(tissue=="normal"){
    dat1 <- dat1 %>% 
      filter(!grepl('Tumour|Cirrhosis|Pancreatitis', tissue)) 
    dat1$tissue <- factor(dat1$tissue,levels=grep("Tumour|Cirrhosis|Pancreatitis",gsub("CD34-megakaryocytes","Megakaryocytes",tissue_order) %>% gsub("CD34-erythroblasts","Erythroblasts",.),invert = TRUE, value = TRUE))
  }else if(tissue=="tumour"){
    dat1 <- dat1 %>% filter(grepl('Tumour|Cirrhosis|Pancreatitis', tissue))
    dat1$tissue <- factor(dat1$tissue,levels=grep("Tumour|Cirrhosis|Pancreatitis",gsub("CD34-megakaryocytes","Megakaryocytes",tissue_order) %>% gsub("CD34-erythroblasts","Erythroblasts",.), value = TRUE))
  }else{
    dat1$tissue <- factor(dat1$tissue,levels=gsub("CD34-megakaryocytes","Megakaryocytes",tissue_order) %>% gsub("CD34-erythroblasts","Erythroblasts",.))
  }
  
  dat1 <- dat1[order(dat1$tissue),]
  smp_order <- data.frame(smp=colnames(dat1)[-c(1:4)])
  smp_order$tissue <- factor(gsub(".*_","",gsub("_mC|_hmC|_amC","",smp_order$smp)), levels=tissue_order)
  smp_order <- smp_order[order(smp_order$tissue),]
  dat1 <- dat1[,c(1:4,rownames(smp_order)%>% as.numeric()+4)]
  colnames(dat1) <- gsub("-","_",colnames(dat1))
  tmp1 <- dat1[,grep("_hmC",colnames(dat1))];colnames(tmp1) <- gsub("_hmC","",colnames(tmp1)); tmp1 <- rename_columns(tmp1); colnames(tmp1) <- paste0(colnames(tmp1),"_5hmC")
  tmp2 <- dat1[,grep("_mC",colnames(dat1))];colnames(tmp2) <- gsub("_mC","",colnames(tmp2)); tmp2 <- rename_columns(tmp2); colnames(tmp2) <- paste0(colnames(tmp2),"_5mC")
  tmp3 <- dat1[,grep("_amC",colnames(dat1))];colnames(tmp3) <- gsub("_amC","",colnames(tmp3)); tmp3 <- rename_columns(tmp3); colnames(tmp3) <- paste0(colnames(tmp3),"_5mC+5hmC")
  
  dat1 <- cbind(dat1[1:4],tmp1,tmp2,tmp3)
  dat1 <- dat1 %>% 
    select(-contains("Tumour"), -contains("Cirrhosis"), -contains("Pancreatitis"))
  pheatmap(dat1[,c(grep("_5hmC",colnames(dat1)),grep("_5mC",colnames(dat1)),grep("_5mC+5hmC",colnames(dat1)) )],
           cluster_rows = FALSE,cluster_cols = FALSE, show_rownames = FALSE,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           filename = paste0(dmr,gsub(".txt",".png",infile)), width = 20, height = 6, fontsize_col = 5)
  pheatmap(dat1[,c(grep("_5hmC",colnames(dat1)),grep("_5mC",colnames(dat1)),grep("_5mC+5hmC",colnames(dat1)) )],
           cluster_rows = FALSE,cluster_cols = FALSE, show_rownames = FALSE,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           filename = paste0(dmr,gsub(".txt",".pdf",infile)), width = 20, height = 6, fontsize_col = 5)
}

plot_meth(infile="dmr15_TAPSbeta.Top300_block_DMBs.all.meth.txt",
          dmr_bed="dmr15_TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
          dmr="Hypo",
          tissue_order = tissue_order,tissue="normal", selN = 50)
plot_meth(infile="dmr15_CAPS.Top300_block_DMBs.all.meth.txt",
          dmr_bed="dmr15_CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
          dmr="Hyper",
          tissue_order = tissue_order,tissue="normal", selN = 50)


#### Plot mC / hmC DMB around tissue specific genes ####
library(ggplot2)
library(forecast)
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/mC_vs_hmC")
dat <- read.table("all_gene_updown10k.vs.dmb.txt", header=TRUE, sep="\t")
tissue_pairs <- read.table("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/PanglaoDB_immune_GTEx_Tissue.pairs.txt",header=FALSE, sep="\t",col.names = c("expr_tissue","dmb_tissue"))
tissue_pairs <- tissue_pairs[grep("Erythroblasts|not_specific", tissue_pairs$dmb_tissue,invert = TRUE),]
res <- data.frame()
ma_order<-6
for(sel_tissue_expr in tissue_pairs$expr_tissue){
  sel_tissue_dmb <- tissue_pairs$dmb_tissue[tissue_pairs$expr_tissue==sel_tissue_expr]
  sta <- merge(dat %>% filter(str_detect(expr_tissue, sel_tissue_expr)) %>% filter(str_detect(dmb_tissue,sel_tissue_dmb)) %>% filter(Hyper_CAPS==1) %>% 
                 select(gene,idx) %>% unique() %>% select(idx) %>% table() %>% as.data.frame(),
               dat %>% filter(str_detect(expr_tissue, sel_tissue_expr)) %>% filter(str_detect(dmb_tissue,sel_tissue_dmb)) %>% filter(Hypo_TAPSbeta==1) %>% 
                 select(gene,idx) %>% unique() %>% select(idx) %>% table() %>% as.data.frame(),by="idx", all.x = TRUE, all.y=TRUE) %>%
    merge(dat %>% filter(str_detect(expr_tissue, "not_specific")) %>% filter(str_detect(dmb_tissue,sel_tissue_dmb)) %>% filter(Hyper_CAPS==1) %>% 
            select(gene,idx) %>% unique() %>% select(idx) %>% table() %>% as.data.frame(),by="idx", all.x = TRUE, all.y=TRUE ) %>%
    merge(dat %>% filter(str_detect(expr_tissue, "not_specific")) %>% filter(str_detect(dmb_tissue,sel_tissue_dmb)) %>% filter(Hypo_TAPSbeta==1) %>% 
            select(gene,idx) %>% unique() %>% select(idx) %>% table() %>% as.data.frame(),by="idx", all.x = TRUE, all.y=TRUE )
  features <- c("Hyper_CAPS","Hypo_TAPSbeta",
             paste0("not specific: ",c("Hyper_CAPS","Hypo_TAPSbeta")))
  colnames(sta) <- c("idx",paste0("num:",features))
  sta <- sta[order(sta$idx),]

  n1 <- dat %>% filter(str_detect(expr_tissue, sel_tissue_expr)) %>% select(gene) %>% unique() %>% nrow()
  n2 <- dat %>% filter(str_detect(expr_tissue, "not_specific")) %>% select(gene) %>% unique() %>% nrow()
  sta <- cbind(sta,sta[,2]/n1, sta[,3]/n1, sta[,4]/n2, sta[,5]/n2) %>% as.data.frame()
  colnames(sta)[6:9] <- gsub("num:","pct:",colnames(sta)[c(2:5)])
  sta[is.na(sta)]<-0
  ## use smooth average
  sta[,10:13] <- apply(sta[,6:9],2,function(x)c(x[1:(ma_order/2)],
                                                ma(x,order=ma_order)[-c(1:(ma_order/2),(length(x)-(ma_order/2)):length(x))],
                                                x[(length(x)-(ma_order/2)):length(x)]))
  colnames(sta)[10:13] <- gsub("num:","pct_ma:",colnames(sta)[c(2:5)])
  sta$tissue <- sel_tissue_expr
  sta$n_ts_gene <- n1
  res <- rbind(res,sta)
}
dim(res)

# by tissue type
res_long <- res %>% 
  select(contains("pct_ma") | contains("tissue") | contains("idx")) %>%
  pivot_longer(cols = contains("pct"),
               names_to = "gene_dmb",
               values_to = "pct")
seltissue <- c("Liver","Kidney","Ovary")
seltissue <- unique(res_long$tissue)
p1 <- res_long[res_long$tissue %in% seltissue,] %>% filter(str_detect("pct_ma:Hyper_CAPS",gene_dmb)) %>%
  ggplot(aes(x=idx,y=pct,group=tissue,color=tissue))+
  geom_line() + 
  theme_bw() + ylab("Hyper_CAPS")+
  scale_x_discrete(breaks = c(1,41,80,120),labels = c("-10K","TSS","TES","10K"))+
  scale_color_manual(values=colorRampPalette(brewer.pal(8, "Dark2"))(nrow(tissue_pairs)))


p2 <- res_long[res_long$tissue %in% seltissue,] %>% filter(str_detect("pct_ma:Hypo_TAPSbeta",gene_dmb)) %>%
  ggplot(aes(x=idx,y=pct,group=tissue,color=tissue))+
  geom_line() + 
  theme_bw() + ylab("Hypo_TAPSbeta")+
  scale_x_discrete(breaks = c(1,41,80,120),labels = c("-10K","TSS","TES","10K"))+
  scale_color_manual(values=colorRampPalette(brewer.pal(8, "Dark2"))(nrow(tissue_pairs)))
pdf("all_gene_updown10k.vs.dmb.alltissue.pdf",width = 13, height = 5)
cowplot::plot_grid(p1,p2,nrow=1)
dev.off()



res_wide <- res_long %>%
  pivot_wider(
    id_cols = c("tissue", "gene_dmb"), 
    names_from = "idx",
    values_from = "pct"
  ) %>% as.data.frame()
res_wide[is.na(res_wide)] <- 0
hmC_hyper <- res_wide %>% filter(str_detect(gene_dmb,"pct_ma:Hyper_CAPS")) %>% as.data.frame()
rownames(hmC_hyper)<- hmC_hyper$tissue
pheatmap(hmC_hyper[,-c(1:2)],cluster_rows = FALSE, cluster_cols = FALSE, scale = "row")


mC_hypo <- res_wide %>% filter(str_detect(gene_dmb,"pct_ma:Hypo_TAPSbeta")) %>% as.data.frame()
rownames(mC_hypo)<- mC_hypo$tissue
pheatmap(mC_hypo[,-c(1:2)],cluster_rows = FALSE, cluster_cols = FALSE, scale = "row")

###  combine all tissue specific gene 
res_wide <- res %>% 
  select(contains("num:Hyper_CAPS") | contains("num:Hypo_TAPSbeta") | contains("idx")|contains("tissue")|contains("n_ts_gene")) %>%
  pivot_wider(
    id_cols = c("tissue", "n_ts_gene"), 
    names_from = "idx",
    values_from = c("num:Hyper_CAPS","num:Hypo_TAPSbeta")
  )

res_wide_sum <- as.data.frame(apply(res_wide[,-c(1:2)],2,sum)/sum(res_wide[,2]))
colnames(res_wide_sum) <- "pct"
res_wide_sum$idx <- gsub(".*_","",rownames(res_wide_sum)) %>% as.character() %>% as.numeric()
res_wide_sum$type <- gsub("_[0-9]*","",rownames(res_wide_sum))
res_wide_sum <- res_wide_sum[order(res_wide_sum$type,res_wide_sum$idx),]
res_wide_sum$ma_pct <- c(ma(res_wide_sum$pct[res_wide_sum$type=="num:HyperCAPS"],order=6),
                         ma(res_wide_sum$pct[res_wide_sum$type=="num:HypoTAPSbeta"],order=6))
p <-  ggplot(res_wide_sum,aes(x=idx,y=ma_pct,color=type))+
  geom_line() + 
  theme_bw() + 
  facet_wrap(~type,nrow=2)+
  scale_x_continuous(breaks = c(1,41,80,120),labels = c("-10K","TSS","TES","10K"))+
  scale_color_manual(values=brewer.pal(2, "Dark2")) +
  ylab("# TSG with DMB / # TSG (TSG:Tissue specific gene)\nall_TSG")
ggsave("all_gene_updown10k.vs.dmb.combine_tissue.pdf",width = 5, height = 6)

immune_cell <- c("Eosinophils","Megakaryocytes","Monocytes","Neutrophils","B Cells","NK Cells","T Cytotoxic Cells","T Helper Cells")
res_wide_sel <- res_wide[!res_wide$tissue %in% c("Brain"),]
res_wide_sum <- as.data.frame(apply(res_wide_sel[,-c(1:2)],2,sum)/sum(res_wide_sel[,2]))
colnames(res_wide_sum) <- "pct"
res_wide_sum$idx <- gsub(".*_","",rownames(res_wide_sum)) %>% as.character() %>% as.numeric()
res_wide_sum$type <- gsub("_[0-9]*","",rownames(res_wide_sum))
res_wide_sum <- res_wide_sum[order(res_wide_sum$type,res_wide_sum$idx),]
res_wide_sum$ma_pct <- c(ma(res_wide_sum$pct[res_wide_sum$type=="num:HyperCAPS"],order=6),
                         ma(res_wide_sum$pct[res_wide_sum$type=="num:HypoTAPSbeta"],order=6))
res_wide_sum$type <- gsub("num:","",res_wide_sum$type)
p <-  ggplot(res_wide_sum,aes(x=idx,y=ma_pct,color=type))+
  geom_line() + 
  theme_minimal() +
  facet_wrap(~type,nrow=1)+
  scale_x_continuous(breaks = c(1,41,80,120),labels = c("-10K","TSS","TES","10K"))+
  scale_color_manual(values=brewer.pal(3, "Dark2")[c(1:2)]) +
  ylab("# TSG with DMB / # TSG \n(TSG:Tissue specific gene) no_brain") +
  theme(legend.position = "none")+
  xlab("")
ggsave("all_gene_updown10k.vs.dmb.combine_tissue_nobrain.pdf",p,width = 7, height = 3)
ggsave("all_gene_updown10k.vs.dmb.combine_tissue_nobrain.png",p,width = 7, height = 3)

res_wide_sel <- res_wide[!res_wide$tissue %in% c(immune_cell,"Brain"),]
res_wide_sum <- as.data.frame(apply(res_wide_sel[,-c(1:2)],2,sum)/sum(res_wide_sel[,2]))
colnames(res_wide_sum) <- "pct"
res_wide_sum$idx <- gsub(".*_","",rownames(res_wide_sum)) %>% as.character() %>% as.numeric()
res_wide_sum$type <- gsub("_[0-9]*","",rownames(res_wide_sum))
res_wide_sum <- res_wide_sum[order(res_wide_sum$type,res_wide_sum$idx),]
res_wide_sum$ma_pct <- c(ma(res_wide_sum$pct[res_wide_sum$type=="num:HyperCAPS"],order=6),
                         ma(res_wide_sum$pct[res_wide_sum$type=="num:HypoTAPSbeta"],order=6))
p <-  ggplot(res_wide_sum,aes(x=idx,y=ma_pct,color=type))+
  geom_line() + 
  theme_bw() + 
  facet_wrap(~type,nrow=2)+
  scale_x_continuous(breaks = c(1,41,80,120),labels = c("-10K","TSS","TES","10K"))+
  scale_color_manual(values=brewer.pal(2, "Dark2")) +
  ylab("# TSG with DMB / # TSG (TSG:Tissue specific gene)\nno_immune_cell_brain")
ggsave("all_gene_updown10k.vs.dmb.combine_tissue_noimmunecellbrain.pdf",width = 5, height = 6)

res_wide_sel <- res_wide[res_wide$tissue %in% immune_cell,]
res_wide_sum <- as.data.frame(apply(res_wide_sel[,-c(1:2)],2,sum)/sum(res_wide_sel[,2]))
colnames(res_wide_sum) <- "pct"
res_wide_sum$idx <- gsub(".*_","",rownames(res_wide_sum)) %>% as.character() %>% as.numeric()
res_wide_sum$type <- gsub("_[0-9]*","",rownames(res_wide_sum))
res_wide_sum <- res_wide_sum[order(res_wide_sum$type,res_wide_sum$idx),]
res_wide_sum$ma_pct <- c(ma(res_wide_sum$pct[res_wide_sum$type=="num:HyperCAPS"],order=6),
                         ma(res_wide_sum$pct[res_wide_sum$type=="num:HypoTAPSbeta"],order=6))
p <-  ggplot(res_wide_sum,aes(x=idx,y=ma_pct,color=type))+
  geom_line() + 
  theme_bw() + 
  facet_wrap(~type,nrow=2)+
  scale_x_continuous(breaks = c(1,41,80,120),labels = c("-10K","TSS","TES","10K"))+
  scale_color_manual(values=brewer.pal(2, "Dark2")) +
  ylab("# TSG with DMB / # TSG (TSG:Tissue specific gene)\nimmune_cell_only")
ggsave("all_gene_updown10k.vs.dmb.combine_tissue_immunecell.pdf",width = 5, height = 6)