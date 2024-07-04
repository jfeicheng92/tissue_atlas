setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/data_for_jingfei")
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
library(ggsci)
library(ggbeeswarm)


#### correlation between meth. and expression ####
tissue_specific <- fread("tissue_specific_pearson_corrs.csv") 
colnames(tissue_specific)[-1] <- paste0(colnames(tissue_specific)[-1], ".tissue_specific")
tissue_specific_mat <- tissue_specific %>% 
  dplyr::rename(tissue=V1,
         mC_promoter=MANE_genes.nearest_promoter.200bp_padding.CpG_only.TAPSbeta.tissue_specific,
         hmC_promoter=MANE_genes.nearest_promoter.200bp_padding.CpG_only.CAPS.tissue_specific,
         mC_enhancer=MANE_genes.nearest_enhancer.200bp_padding.CpG_only.TAPSbeta.tissue_specific,
         hmC_enhancer=MANE_genes.nearest_enhancer.200bp_padding.CpG_only.CAPS.tissue_specific,
         mC_genebody=MANE.GRCh38.v1.0.genes.TAPSbeta.tissue_specific,
         hmC_genebody=MANE.GRCh38.v1.0.genes.CAPS.tissue_specific) %>%
  melt(id.vars="tissue")
tissue_specific_mat$exp <- gsub(".*_","",tissue_specific_mat$variable)
tissue_specific_mat$feature <- gsub("_.*","",tissue_specific_mat$variable)
tissue_specific_mat$exp <- factor(tissue_specific_mat$exp, levels=c("genebody","enhancer","promoter"))
p <-  tissue_specific_mat%>%
  ggplot(aes(x=exp,y=value,fill=feature)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  geom_beeswarm(dodge.width=0.9, alpha=0.5, size=0.8)+
  theme_classic() + geom_hline(yintercept = 0, linetype="longdash") +
  theme(legend.position = "bottom")+
  ylab("Pearson's correlation with gene expression") + xlab("features") +
  ylim(-0.7,0.7) +
  scale_fill_manual(values=brewer.pal(4,"RdBu")[c(1,4)])
ggsave("pearson_correlation_with_geneexpre.pdf", p,width = 3, height = 3)

for(infile in c("housekeeping_pearson_corrs.csv","all_gene_pearson_corrs.csv", "tissue_specific_pearson_corrs.csv") ){
  tag <- gsub("_person_corrs.csv","",infile)
  dat <- fread(infile) 
  outfix <- gsub(".csv",".pdf",infile)
  dat_mat <- dat %>% 
    dplyr::rename(tissue=V1,
                  mC_promoter=MANE_genes.nearest_promoter.200bp_padding.CpG_only.TAPSbeta,
                  hmC_promoter=MANE_genes.nearest_promoter.200bp_padding.CpG_only.CAPS,
                  mC_enhancer=MANE_genes.nearest_enhancer.200bp_padding.CpG_only.TAPSbeta,
                  hmC_enhancer=MANE_genes.nearest_enhancer.200bp_padding.CpG_only.CAPS,
                  mC_genebody=MANE.GRCh38.v1.0.genes.TAPSbeta,
                  hmC_genebody=MANE.GRCh38.v1.0.genes.CAPS) %>%
    melt(id.vars="tissue")
  dat_mat$exp <- gsub(".*_","",dat_mat$variable)
  dat_mat$feature <- gsub("_.*","",dat_mat$variable)
  dat_mat$exp <- factor(dat_mat$exp, levels=c("genebody","enhancer","promoter"))
  p <-  dat_mat%>%
    ggplot(aes(x=exp,y=value,fill=feature)) +
    geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
    geom_beeswarm(dodge.width=0.9, alpha=0.5, size=0.8)+
    theme_classic() + geom_hline(yintercept = 0, linetype="longdash") +
    theme(legend.position = "bottom")+
    ylab("Pearson's correlation with gene expression") + xlab("features") +
    ylim(-0.7,0.7) +
    scale_fill_manual(values=brewer.pal(4,"RdBu")[c(1,4)])
  ggsave(outfix, p,width = 3, height = 3)
}




#### correlation between predicted expression and expression ####
dat <- fread("model_32.LOO_all_model_predictions.csv") %>% as.data.frame()
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")[c(2:20)]
dat$tissue_x <- factor(dat$tissue_x, levels = tissue_order)
dat <- dat[order(dat$tissue_x), ]
# by tissue type 
anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(dat$tissue_x))))[c(2:20)]
p <- dat[!grepl("Breast",dat$tissue_x),] %>% 
  ggplot(aes(x=logTPM, y=model_preds,color=tissue_x)) + geom_point(alpha = 1/10) +
  theme_classic() +
  facet_wrap( ~ tissue_x, nrow=5) +
  scale_color_manual(values=anno_colors)

ggsave("model_32.LOO_all_model_predictions_all_tissue.other.pdf",p,width = 10,height = 10)
p <- dat[grepl("Breast",dat$tissue_x),] %>% 
  ggplot(aes(x=logTPM, y=model_preds,color=tissue_x)) + geom_point(alpha = 1/10) +
  theme_classic() +
  facet_wrap( ~ tissue_x, nrow=1) +
  scale_color_manual(values=anno_colors) +
  theme(legend.position = "none")

ggsave("model_32.LOO_all_model_predictions_sel_tissue.pdf",p,width = 2,height = 2)

# correlations <- dat %>%
#   group_by(tissue_x) %>%
#   summarize(correlation = cor(logTPM,model_preds))
# cor(dat$logTPM, dat$model_preds)
# # on selected gene
# p <- dat[dat$core_gene_x =="CD300A",] %>%
#   ggplot(aes(x=logTPM, y=model_preds,color=tissue_x)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values=anno_colors)
# p
# ggsave("model_32.LOO_all_model_predictions.CD300A.pdf",p,width = 4,height = 3)


normal <- fread("model_32.healthy_tissue.LOO_scores.Pearson_Correlation.csv")
normal <- normal %>% melt(id.vars="tissue_x") %>% 
  filter(tissue_x==variable)
normal$tissue_x <- factor(normal$tissue_x, levels = tissue_order)
normal <- normal[order(normal$tissue_x),]
p <- normal  %>%
  ggplot(aes(x=tissue_x, y=value, color=tissue_x)) +
  geom_point(stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None") +
  ylim(0,1) +
  scale_color_manual(values=anno_colors)
ggsave("model_32.healthy_tissue.LOO_scores.Pearson_Correlation.pdf", p,width = 6, height = 2.5)


#### meth input for predictions ####
dat <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/meth_around_gene/MANE.GRCh38.v1.0.refseq_genomic.gene.flank20000.nwin20.all_meth.closest.1000.txt")
prefix <- "MANE.GRCh38.v1.0.refseq_genomic.gene.flank20000.nwin10.all_meth.closest.1000"
selgene <- "CD300A"
seltissue <- "C505016_Brain"
outfix <- paste(prefix,seltissue,selgene,"pdf",sep=".")
dat %>% 
  filter(gene == selgene) %>%
  select(contains(seltissue)|contains("gene")|contains("idx")) %>% 
  select(ends_with("C")) %>%
  t() %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,filename = outfix, width = 5, height = 0.5, border_color = NA)
seltissue <- "UKVAC-049-3_CD4-T-cells"
outfix <- paste(prefix,seltissue,selgene,"pdf",sep=".")
dat %>% 
  filter(gene == selgene) %>%
  select(contains(seltissue)|contains("gene")|contains("idx")) %>% 
  select(ends_with("C")) %>%
  t() %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,filename = outfix, width = 5, height = 0.5, border_color = NA)
seltissue <- "CD563257_Liver"
outfix <- paste(prefix,seltissue,selgene,"pdf",sep=".")
dat %>% 
  filter(gene == selgene) %>%
  select(contains(seltissue)|contains("gene")|contains("idx")) %>% 
  select(ends_with("C")) %>%
  t() %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,filename = outfix, width = 5, height = 0.5, border_color = NA)


dat %>% 
  filter(gene == selgene) %>%
  select(ends_with("_mC")|contains("gene")|contains("idx")) %>% 
  select(ends_with("C")) %>%
  t() %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,border_color = NA)
dat %>% 
  filter(gene == selgene) %>%
  select(ends_with("_hmC")|contains("gene")|contains("idx")) %>% 
  select(ends_with("C")) %>%
  t() %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,border_color = NA)







#### feature importance ####
hmC <- fread("model_32.feature_importance.hmC.csv", header = TRUE) %>% as.data.frame(); rownames(hmC) <- hmC$V1
mC <- fread("model_32.feature_importance.mC.csv", header = TRUE) %>% as.data.frame(); rownames(mC) <- mC$V1
umC <- fread("model_32.feature_importance.umC.csv", header = TRUE) %>% as.data.frame(); rownames(umC) <- umC$V1
hmC <- hmC[order(hmC$V1),]
mC <- mC[order(mC$V1),]

p1 <- mC %>%
  melt(id.vars=c("V1")) %>%
  ggplot(aes(x = V1, y = variable)) +
  geom_tile(aes(fill = value)) +
  coord_flip() +
  theme(legend.position = "right") +
  xlab("Tissue") + ylab("mC") +
  theme(axis.text=element_text(size=10)) +
  scale_fill_gradientn(limits = c(0, 0.1), colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10)) +
  scale_y_discrete(breaks=c(0,19,39,59),labels = c("-20K","TSS","TES","20K")) +
  theme_minimal()

p2 <- hmC %>%
  melt(id.vars=c("V1")) %>%
  ggplot(aes(x = V1, y = variable)) +
  geom_tile(aes(fill = value)) +
  coord_flip() +
  theme(legend.position = "right") +
  xlab("Tissue") + ylab("hmC") +
  theme(axis.text=element_text(size=10)) +
  scale_fill_gradientn(limits = c(0, 0.1), colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10)) +
  scale_y_discrete(breaks=c(0,19,39,59),labels = c("-20K","TSS","TES","20K")) +
  theme_minimal()

pdf("model_32.feature_importance.pdf",width=10, height = 5)
cowplot::plot_grid(p1,p2,nrow=2)
dev.off()


#### Pearson's correlation for cancer ####
cancer <- fread("model_32.cancer_GE.Z-score_Pearson_Correlation_df2.csv") %>% as.data.frame()
rownames(cancer) <- cancer[,1]
pheatmap(cancer[,-1],cluster_rows = FALSE,cluster_cols = FALSE,filename = "model_32.cancer_GE.Pearson_Correlation_df2.pdf")

