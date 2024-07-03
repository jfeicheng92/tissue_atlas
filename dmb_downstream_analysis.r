```{bash}
cp -r /gpfs3/well/ludwig/users/dyp502/tissue_atlas_v3/dmr14 raw
cp -r /gpfs3/well/ludwig/users/dyp502/tissue_atlas_v3/dmr15 raw
# workingdir: /users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr14
module load BEDTools
tapsb_dmb=TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
caps_dmb=CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
gsize=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.fai
ln -s raw/$tapsb_dmb
ln -s raw/$caps_dmb
for i in Hyper Hypo
do
  n=`head -1 $tapsb_dmb|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep -P "Hyper_or_Hypo|selected_tissue" |cut -f1|tr '\n' ',' |sed 's/,$//g'`
  cut -f1,$n $tapsb_dmb|grep $i |grep -P "Tumour"|sed 's/_/\t/g' |grep -v start |bedtools sort -i - -faidx <(grep -w chr[0-9,X,Y,M]* $gsize) >${tapsb_dmb/.csv/}.${i}.tumour.bed
  cut -f1,$n $tapsb_dmb|grep $i |grep -vP "Tumour|Cirrhosis|Pancreatitis"|sed 's/_/\t/g' |grep -v start|bedtools sort -i - -faidx <(grep -w chr[0-9,X,Y,M]* $gsize) >${tapsb_dmb/.csv/}.${i}.normal.bed
  n=`head -1 $tapsb_dmb|sed 's/\t/\n/g'|awk '{print NR"\t"$0}'|grep -P "Hyper_or_Hypo|selected_tissue" |cut -f1|tr '\n' ',' |sed 's/,$//g'`
  cut -f1,$n $caps_dmb |grep $i |grep -P "Tumour"|sed 's/_/\t/g' |grep -v start |bedtools sort -i - -faidx <(grep -w chr[0-9,X,Y,M]* $gsize) >${caps_dmb/.csv/}.${i}.tumour.bed
  cut -f1,$n $caps_dmb |grep $i |grep -vP "Tumour|Cirrhosis|Pancreatitis"|sed 's/_/\t/g' |grep -v start|bedtools sort -i - -faidx <(grep -w chr[0-9,X,Y,M]* $gsize) >${caps_dmb/.csv/}.${i}.normal.bed
done

```


#### 0. Heatmap ####

library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(pheatmap)


# WORKDIR <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr14/"
WORKDIR <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr15/"
setwd(WORKDIR)
dir.create("heatmap")
plot_heatmap <- function(infile, dmr, tissue_order, tumour ="no", selN=50, width = 12,height = 7,prefix){
  dat1 <- fread(infile) %>% as.data.frame()
  colnames(dat1) <- gsub("CD34-erythroblasts","Erythroid-precursors",colnames(dat1)); colnames(dat1) <- gsub("CD34-megakaryocytes","Megakaryocytes",colnames(dat1))
  dat1$tissue <- as.factor(dat1$selected_tissue)
  dat1$tissue <- gsub("Erythroblasts","Erythroid-precursors",dat1$tissue)
  
  dat1 <- dat1[grep(dmr,dat1$Hyper_or_Hypo),]
  if(tumour == "no"){
    tissue_order <- c(grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,invert = TRUE, value = TRUE), grep("Tumour|Cirrhosis|Pancreatitis",tissue_order, value = TRUE))
    dat1 <- dat1[grep("Tumour|Cirrhosis|Pancreatitis",dat1$selected_tissue,invert = TRUE),]
    dat1$tissue <- factor(dat1$tissue,levels=grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,invert = TRUE, value = TRUE))
  }else{
    tissue_order <- c(grep("Tumour|Cirrhosis|Pancreatitis",tissue_order, value = TRUE),grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,invert = TRUE, value = TRUE))
    dat1 <- dat1[grep("Tumour|Cirrhosis|Pancreatitis",dat1$selected_tissue),]
    dat1$tissue <- factor(dat1$tissue,levels=grep("Tumour|Cirrhosis|Pancreatitis",tissue_order,value = TRUE))
  }
  anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(dat1$tissue))))
  names(anno_colors) <- unique(dat1$tissue)
  anno_colors <- list(type=anno_colors)
  
  dat1 <- dat1[order(dat1$tissue),]
  colnames(dat1) <- gsub("_TAPSbeta_RATE|_CAPS_RATE","",colnames(dat1))
  dat1 <- dat1[order(dat1$tissue),]
  all_dat <- dat1
  dat1 <- all_dat %>%
    group_by(tissue) %>%
    arrange(desc(delta_quants)) %>%
    slice_head(n = selN) %>% as.data.frame()
  
  smp_order <- data.frame(smp=colnames(dat1)[c(2:grep("delta",colnames(dat1))[1]-1)])
  smp_order$tissue <- factor(gsub(".*_","",smp_order$smp), levels=tissue_order)
  smp_order <- smp_order[order(smp_order$tissue),]
  smp_order <- smp_order[complete.cases(smp_order),]
  anno_row <- data.frame(type=dat1$tissue); rownames(anno_row)<- rownames(dat1)
  dat1 <- dat1[dat1$tissue %in% smp_order$tissue,c(1,rownames(smp_order) %>% as.numeric())]
  
  
  
  pheatmap(dat1[,-c(1)],cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
           filename = paste0(prefix,".label.png"))
  pheatmap(dat1[,-c(1)]*100,cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
           annotation_colors=anno_colors, 
           annotation_row = anno_row,  width = width,height = height,labels_col=gsub(".*_","",colnames(dat1)[-1]), 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           filename = paste0(prefix,".png"))
  pheatmap(dat1[,-c(1)]*100,cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
           annotation_colors=anno_colors, 
           annotation_row = anno_row, width = width,height = height ,labels_col=gsub(".*_","",colnames(dat1)[-1]), 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           filename = paste0(prefix,".pdf"))
}

#### healthy on healthy dmb
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
plot_heatmap(infile="TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hypo",tumour="no",tissue_order=tissue_order,prefix = "heatmap/TAPSbeta.all_samples.Healthy.Hypo")
plot_heatmap(infile="CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hyper",tumour="no",tissue_order=tissue_order,prefix="heatmap/CAPS.all_samples.Healthy.Hyper")
plot_heatmap(infile="TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hyper",tumour="no",tissue_order=tissue_order,prefix = "heatmap/TAPSbeta.all_samples.Healthy.Hyper")
plot_heatmap(infile="CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hypo",tumour="no",tissue_order=tissue_order,prefix="heatmap/CAPS.all_samples.Healthy.Hypo")

#### tumour on healthy dmb
tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes",
                  "Brain-Tumour","Breast-Tumour", "Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Colon-Tumour","Stomach-Tumour",
                  "Liver-Cirrhosis","Pancreas-Pancreatitis")
plot_heatmap(infile="TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hypo",tumour="no", tissue_order=tissue_order,prefix = "heatmap/TAPSbeta.all_samples.all.Hypo")
plot_heatmap(infile="CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hyper",tumour="no",tissue_order=tissue_order,prefix="heatmap/CAPS.all_samples.all.Hyper")
plot_heatmap(infile="TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hyper",tumour="no",tissue_order=tissue_order,prefix = "heatmap/TAPSbeta.all_samples.all.Hyper")
plot_heatmap(infile="CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hypo",tumour="no",tissue_order=tissue_order,prefix="heatmap/CAPS.all_samples.all.Hypo")

#### tumour + healthy on tumour dmb
tissue_order <- c("Brain-Tumour","Breast-Tumour", "Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Colon-Tumour","Stomach-Tumour",
                  "Liver-Cirrhosis","Pancreas-Pancreatitis",
                  "Brain","Breast", "Heart",
                  "Kidney","Liver", "Lung","Ovary", 
                  "Pancreas", "Prostate","Colon","Stomach","Esophagus",
                  "Spleen", "CD4-T-cells", "CD8-T-cells", "Neutrophils", "NK-cells", "B-cells", "Eosinophils", "Monocytes", 
                  "Erythroid-precursors", "Megakaryocytes")
plot_heatmap(infile="TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hypo",tumour="yes",tissue_order=tissue_order,width = 15, prefix = "heatmap/TAPSbeta.all_samples.Tumour.Hypo")
plot_heatmap(infile="CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hyper",tumour="yes",tissue_order=tissue_order,width = 15, prefix="heatmap/CAPS.all_samples.Tumour.Hyper")
plot_heatmap(infile="TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hyper",tumour="yes",tissue_order=tissue_order,width = 15, prefix = "heatmap/TAPSbeta.all_samples.Tumour.Hyper")
plot_heatmap(infile="CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
             dmr="Hypo",tumour="yes",tissue_order=tissue_order,width = 15, prefix="heatmap/CAPS.all_samples.Tumour.Hypo")

# tissue_order <- c("Brain-Tumour","Breast-Tumour", "Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Colon-Tumour","Stomach-Tumour",
#                   "Brain","Breast", "Heart",
#                   "Kidney","Liver", "Lung","Ovary",
#                   "Pancreas", "Prostate","Colon","Stomach","Esophagus",
#                   "Spleen", "CD4-T-cells", "CD8-T-cells", "Neutrophils", "NK-cells", "B-cells", "Eosinophils", "Monocytes",
#                   "Erythroid-precursors", "Megakaryocytes","Liver-Cirrhosis","Pancreas-Pancreatitis")
# plot_heatmap(infile="TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
#              dmr="Hypo",tumour="yes",tissue_order=tissue_order,width = 15, prefix = "heatmap/TAPSbeta.all_samples.Tumour_pre.Hypo")
# plot_heatmap(infile="CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
#              dmr="Hyper",tumour="yes",tissue_order=tissue_order,width = 15, prefix="heatmap/CAPS.all_samples.Tumour_pre.Hyper")
# 


# plot_heatmap(infile="raw/TAPSbeta.healthy_high_cov_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
#              dmr="Hypo",tumour="no",tissue_order=tissue_order,prefix = "heatmap/TAPSbeta.all_samples.justHealthy.Hypo")
# plot_heatmap(infile="raw/CAPS.healthy_high_cov_samples.Tissue_Group_DMBs.top300.CAPS_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
#              dmr="Hyper",tumour="no",tissue_order=tissue_order,prefix="heatmap/CAPS.all_samples.justHealthy.Hyper")
# 




#### 1. Feature enrichment ####

```{bash}
tapsb_dmb=TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
caps_dmb=CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
gsize=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.fai

mkdir feature_enrichment
bed=(${tapsb_dmb/.csv/}.Hyper.tumour.bed ${caps_dmb/.csv/}.Hyper.tumour.bed ${tapsb_dmb/.csv/}.Hyper.normal.bed ${caps_dmb/.csv/}.Hyper.normal.bed ${tapsb_dmb/.csv/}.Hypo.tumour.bed  ${caps_dmb/.csv/}.Hypo.tumour.bed  ${tapsb_dmb/.csv/}.Hypo.normal.bed  ${caps_dmb/.csv/}.Hypo.normal.bed ${tapsb_dmb/.csv/}.tissue_specific_gene.bed ${caps_dmb/.csv/}.tissue_specific_gene.bed )

# shuffle using genome as background
for j in ${bed[@]}
do
    exp=`echo $j |cut -d '_' -f1`
    temp=feature_enrichment/${j/.bed/}.shuffle.genome.bed
    touch $temp
    for seed in `seq 1 1 100`
    do
    bedtools shuffle -i $j -g <(grep -w chr[0-9]* $gsize) -seed $seed |\
        awk -v s=$seed 'BEGIN{OFS="\t"}{print $0,s}'
    done >$temp
    awk -v s=0 'BEGIN{OFS="\t"}{print $0,s}' $j >>$temp
done


feature=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/all_feature.bed
echo -e "bed\ttest\tnum\tfeature" >feature_enrichment/all_feature.txt
for j in ${bed[@]}
do
    temp=feature_enrichment/${j/.bed/}.shuffle.genome.bed
    out=feature_enrichment/${j/.bed/}.shuffle.feature
    bedtools intersect \
    -a <(cat $feature|grep -w chr[0-9]*|bedtools sort -i - -g $gsize) \
    -b <(cat $temp|grep -w chr[0-9]*|bedtools sort -i - -g $gsize) -wa -wb -g <(grep -w chr[0-9]* $gsize ) -sorted |
    tee ${out}.bed|
    cut -f4,10|sort |uniq -c |sed 's/^ \+//g;s/ /\t/g'|sort -k3,3|\
    join -1 1 <(seq 0 1 100|sort -k1,1) -2 3 - -t$'\t' -a1|sort -k1,1n |\
    awk -v n=$j 'BEGIN{OFS="\t"}{print n,$0}'|\
    cat <(echo -e "bed\ttest\tnum\tfeature" ) - |\
    tee ${out}.txt|awk '$1!="bed"' >> feature_enrichment/all_feature.txt
done
```


dmb_dir <- paste0(WORKDIR,"feature_enrichment/")
setwd(dmb_dir)
enrich <- read.table(paste0(dmb_dir,"all_feature.txt"), header=TRUE)
enrich$bed <- gsub(".all_sample.*4cpg.|.bed","",enrich$bed)
enrich <- enrich %>% 
  pivot_wider(id_cols = c(bed,feature),names_from = c(test), values_from =c(num))
  
enrich$p <- apply(enrich, 1, function(x){x <- as.numeric(x[3:103]);sum(x[-1]>x[1])})/100
enrich$odds <- apply(enrich, 1, function(x){x <- as.numeric(x[3:103]);x[1]/mean(x[-1])})
enrich$p_level <- cut(enrich$p,c(0,0.01,1),include.lowest = TRUE,)
enrich$observe <- apply(enrich, 1, function(x){x <- as.numeric(x[3:103]);x[1]})
enrich$random <- apply(enrich, 1, function(x){x <- as.numeric(x[3:103]);mean(x[-1])})

enrich$feature <- gsub(".bed","",enrich$feature)
enrich$feature <- gsub("hg38lift_genome_100_segments/","chrhmm ",enrich$feature)
enrich$feature <- gsub("cpgIsland", "CGI", enrich$feature)
enrich$feature <- gsub("GapArtf", "assembly gaps and alignment artifacts", enrich$feature)
enrich$feature <- gsub("HET", "heterochromatin", enrich$feature)
enrich$feature <- gsub("Quies", "quiescent states", enrich$feature)
enrich$feature <- gsub("Acet", "acetylations marks", enrich$feature)
enrich$feature <- gsub("znf", "ZNF genes", enrich$feature)
enrich$feature <- gsub("ReprPC", "polycomb repressed", enrich$feature)
enrich$feature <- gsub("TSS", "TSS", enrich$feature)
enrich$feature <- gsub("BivProm", "bivalent promoters", enrich$feature)
enrich$feature <- gsub("PromF", "flanking promoters", enrich$feature)
enrich$feature <- gsub("DNase", "only DNase I hypersensitivity", enrich$feature)
enrich$feature <- gsub("EnhA", "strong enhancers", enrich$feature)
enrich$feature <- gsub("EnhWk", "weak enhancers", enrich$feature)
enrich$feature <- gsub("TxEnh", "transcribed candidate enhancers", enrich$feature)
enrich$feature <- gsub("TxEx", "exon", enrich$feature)
enrich$feature <- gsub("TxWk", "weak transcription", enrich$feature)
enrich$feature <- gsub("Tx", "transcription", enrich$feature)

dmb_num <- enrich %>% select(bed,feature,odds,observe,random) %>% melt(id.vars=c("bed","feature","odds"))


dmb_num$odds <- round(dmb_num$odds,2)
dmb_num$odds[dmb_num$variable=="observe"] <- ""
dmb_num$enriched <- ifelse(dmb_num$odds>1,"yes","not")

pdf(paste0(dmb_dir,"feature_enrichment.num_or.pdf"),width = 10, height = 7)
dmb_num[grepl("normal",dmb_num$bed)&!grepl("Satellite_centromeric|lincRNA|Retroposon_SVA",dmb_num$feature),] %>%
  ggplot(aes(x = feature, y = value, fill=variable)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("number of DMB in genomic features") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  geom_text(
    aes(label = odds,y=8000,colour = enriched),
    size = 3, 
    vjust = 0.5, position = position_dodge(0.05)) + 
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
dmb_num[grepl("tissue_specific_gene",dmb_num$bed)&!grepl("Satellite_centromeric|lincRNA|Retroposon_SVA",dmb_num$feature),] %>%
  ggplot(aes(x = feature, y = value, fill=variable)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("number of DMB in genomic features") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  geom_text(
    aes(label = odds,y=2000,colour = enriched),
    size = 3, 
    vjust = 0.5, position = position_dodge(0.05)) + 
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
dmb_num[grepl("tumour",dmb_num$bed)&!grepl("Satellite_centromeric|lincRNA|Retroposon_SVA",dmb_num$feature),] %>%
  ggplot(aes(x = feature, y = value,fill=variable)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("number of DMB in genomic features") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  geom_text(
    aes(label = odds,y=4000,colour = enriched),
    size = 3, 
    vjust = 0.5, position = position_dodge(0.05)) + 
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
dev.off()

pdf(paste0(dmb_dir,"feature_enrichment.num_or.selfeature.pdf"),width = 10, height = 5)
dmb_num[grepl("normal",dmb_num$bed)&grepl("2k|gene$|segments.Enh|repeatmasker$",dmb_num$feature),] %>%
  ggplot(aes(x = feature, y = value, fill=variable)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("number of DMB in genomic features") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  geom_text(
    aes(label = odds,y=8000,colour = enriched),
    size = 3, 
    vjust = 0.5, position = position_dodge(0.05)) + 
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
dmb_num[grepl("tissue_specific_gene",dmb_num$bed)&grepl("2k|gene$|segments.Enh|repeatmasker$",dmb_num$feature),] %>%
  ggplot(aes(x = feature, y = value, fill=variable)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("number of DMB in genomic features") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  geom_text(
    aes(label = odds,y=2000,colour = enriched),
    size = 3, 
    vjust = 0.5, position = position_dodge(0.05)) + 
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
dmb_num[grepl("tumour",dmb_num$bed)&grepl("2k|gene$|segments.Enh|repeatmasker$",dmb_num$feature),] %>%
  ggplot(aes(x = feature, y = value,fill=variable)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("number of DMB in genomic features") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  geom_text(
    aes(label = odds,y=4000,colour = enriched),
    size = 3, 
    vjust = 0.5, position = position_dodge(0.05)) + 
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
dev.off()
pdf(paste0(dmb_dir,"feature_enrichment.num_or.selfeature.pdf"),width = 12, height = 5)
p1 <- dmb_num[grepl("tissue_specific_gene",dmb_num$bed)&grepl("2k|gene$|segments.Enh|repeatmasker$",dmb_num$feature),] %>%
  filter(variable=="random") %>%
  ggplot(aes(x = feature, y = as.numeric(odds), fill=as.numeric(odds)<1)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=1) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom") +
  ylab("odds ratio") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
p2 <- dmb_num[grepl("tumour",dmb_num$bed)&grepl("2k|gene$|segments.Enh|repeatmasker$",dmb_num$feature),] %>%
  filter(variable=="random") %>%
  ggplot(aes(x = feature, y = as.numeric(odds), fill=as.numeric(odds)<1)) +
  geom_bar(stat = "identity",width = 0.6, position=position_dodge(width=0.7)) +
  facet_wrap(~bed,ncol=2) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom") +
  ylab("odds ratio") + xlab("chromatin states")+
  scale_fill_manual(values=c("darkred","grey")) +
  scale_color_manual(values = c("grey", "blue")) +
  coord_flip()
cowplot::plot_grid(p1,p2,ncol=2, rel_widths = c(1.5,2))
dev.off()
write.table(dmb_num,"feature_enrichment.num_or.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)



dmb_dir <- paste0(WORKDIR,"feature_enrichment/")
setwd(dmb_dir)
features <- fread("TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Hypo.tumour.shuffle.feature.bed")
features$V4 <- gsub(".bed","",features$V4)
features$V4 <- gsub("hg38lift_genome_100_segments/","chrhmm ",features$V4)
features$V4 <- gsub("cpgIsland", "CGI", features$V4)
features$V4 <- gsub("GapArtf", "assembly gaps and alignment artifacts", features$V4)
features$V4 <- gsub("HET", "heterochromatin", features$V4)
features$V4 <- gsub("Quies", "quiescent states", features$V4)
features$V4 <- gsub("Acet", "acetylations marks", features$V4)
features$V4 <- gsub("znf", "ZNF genes", features$V4)
features$V4 <- gsub("ReprPC", "polycomb repressed", features$V4)
features$V4 <- gsub("TSS", "TSS", features$V4)
features$V4 <- gsub("BivProm", "bivalent promoters", features$V4)
features$V4 <- gsub("PromF", "flanking promoters", features$V4)
features$V4 <- gsub("DNase", "only DNase I hypersensitivity", features$V4)
features$V4 <- gsub("EnhA", "strong enhancers", features$V4)
features$V4 <- gsub("EnhWk", "weak enhancers", features$V4)
features$V4 <- gsub("TxEnh", "transcribed candidate enhancers", features$V4)
features$V4 <- gsub("TxEx", "exon", features$V4)
features$V4 <- gsub("TxWk", "weak transcription", features$V4)
features$V4 <- gsub("Tx", "transcription", features$V4)

p <- features[!grepl("Satellite_centromeric|lincRNA|Retroposon_SVA",features$V4),] %>%
  filter(V10=="0") %>% 
  select(V4,V5,V6,V7,V9) %>% 
  distinct() %>% select(V4,V9) %>%
  table() %>% as.data.frame() %>%
  ggplot(aes(x = V9, y = V4)) +
  geom_tile(aes(fill = Freq),color="white") +
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)) +
  geom_text(aes(label = sprintf("%d", Freq)), color = "black", size = 3) +
  guides(col = guide_legend(ncol = 3)) +
  xlab("Tissue") + ylab("Feature") +
  theme(axis.text=element_text(size=10), legend.position = "right") +
  theme_minimal()
ggsave("TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Hypo.tumour.shuffle.feature.pdf",
       p,width = 15, height = 8)

p <- features[grepl("PMD",features$V4),] %>%
  filter(V10=="0") %>% 
  select(V4,V5,V6,V7,V9) %>% 
  distinct() %>% select(V4,V9) %>%
  table() %>% as.data.frame() %>%
  ggplot(aes(x = V9, y = Freq/300)) +
  geom_bar(stat="identity") +
  coord_flip()+
  theme_minimal() +
  ylab("% DMB overlap with PMD")
ggsave("TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.Hypo.tumour.shuffle.selfeature.pdf",
       p,width = 6, height = 4)
#### 2. DMB vs. gene enrich ####
```{bash}
##### dmb genelist #####
module load BEDTools
mkdir genelist
faidx=/gpfs3/well/ludwig/users/cfo155/tissueMap/others/mC/resource/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.fai
gene=/gpfs3/well/ludwig/users/cfo155/tissueMap/others/mC/resource/MANE.GRCh38.v1.0.refseq_genomic.gene.bed

tapsb_dmb=TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
caps_dmb=CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv

cat $caps_dmb|awk 'BEGIN{OFS="\t"}{print $1,($(NF-1)),($NF)}'|\
    tail -n +2 |sed 's/_/\t/g'|bedtools sort -i - -g $faidx |\
    bedtools closest -a - -b <(awk 'BEGIN{OFS="\t"}{if($6=="+")print $1,$2,$2+1,$4;else print $1,$3-1,$3,$4 }' $gene|\
    bedtools sort -i - -g $faidx) -g $faidx -wa -wb -d|\
    tee ${caps_dmb}.all_genelist|\
    awk '$10 < 200000'|cut -f4,5,9|sort -u|awk 'BEGIN{OFS="\t"}{print $3>"genelist/caps"$1"_"$2"_genelist"}'

cat $tapsb_dmb|awk 'BEGIN{OFS="\t"}{print $1,($(NF-1)),($NF)}'|\
    tail -n +2 |sed 's/_/\t/g'|bedtools sort -i - -g $faidx |\
    bedtools closest -a - -b <(awk 'BEGIN{OFS="\t"}{if($6=="+")print $1,$2,$2+1,$4;else print $1,$3-1,$3,$4 }' $gene|\
    bedtools sort -i - -g $faidx) -g $faidx -wa -wb -d|\
    tee ${tapsb_dmb}.all_genelist|\
    awk '$10 < 200000'|cut -f4,5,9|sort -u|awk 'BEGIN{OFS="\t"}{print $3>"genelist/tapsbeta"$1"_"$2"_genelist"}'


```

ts_gene_enrich <- function(target_gene, tissue_specific_gene, all_gene){
  # Read data from files
  
  x1 <- read.table(target_gene, header = FALSE,col.names = c("gene"))
  ts_gene <- read.table(tissue_specific_gene, header=FALSE, sep="\t", col.names = c("gene","tissue"))
  total_gene <- read.table(all_gene, header=FALSE, sep="\t", col.names = c("chr","start","end","gene","score","strand"))
  enrich_list <- NULL
  for(tissue in unique(ts_gene$tissue)){
    x2 <- ts_gene[ts_gene$tissue==tissue, ]
    # Calculate overlap and lengths
    overlap <- length(intersect(x1$gene, x2$gene))
    overlap_gene <- paste(intersect(x1$gene, x2$gene), collapse = ";")
    nx1 <- length(x1$gene)
    nx2 <- length(x2$gene)
    # Create a 2x2 contingency table
    contingency_table <- matrix(c(
      overlap, nx1 - overlap,
      nx2 - overlap, nrow(total_gene) + overlap - nx1 - nx2
    ), nrow = 2, dimnames = list(
      list1 = c("list1", "nlist1"),
      list2 = c("list2", "nlist2")
    ))
    
    # Perform Fisher's exact test
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    
    # Extract relevant statistics
    p_value <- fisher_result$p.value
    odds_ratio <- fisher_result$estimate %>% round(2)
    
    # Store results in a list
    enrich_list <- c(
      enrich_list,
      list(
        basename(target_gene),
        tissue,
        overlap_gene,
        overlap,
        nx1,
        nx2,
        odds_ratio,
        p_value
        
      )
    )
  }
  enrich <-  data.frame(matrix(do.call(rbind, enrich_list)%>%unlist,
                               nrow=length(enrich_list)/8,byrow=TRUE))
  colnames(enrich)  <- c("list1","list2","overlap_gene","overlap","nlist1","nlist2", "odds_ratio","pvalue")
  enrich[, -c(1:3)] <- apply(enrich[, -c(1:3)], 2, function(x) {
    as.numeric(as.character(x))
  })
  enrich <- enrich[order(enrich$odds_ratio, decreasing = TRUE),]
  return(enrich)
}

library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
##### tumour-specific gene #####
resource_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/"
genelist_path <- paste0(WORKDIR,"genelist/")
setwd(genelist_path)
genelists <- grep("Tumour",list.files(pattern="*genelist$", path = genelist_path),value=TRUE, invert = FALSE)

result_list <- lapply(genelists, function(genelist) {
  ts_gene_enrich(
    target_gene = paste0(genelist_path,genelist),
    tissue_specific_gene = paste0(resource_path,"toil_tumour_specific_gene.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed")
  )
})


result <- do.call(rbind, result_list)
result <- data.frame(result)
write.table(result,"toil_tumour_specific_gene.dmb_enrichment.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)


for(dmb in c("capsHyper_","tapsbetaHypo_")){
  dat <-  result[grep(dmb,result$list1),]
  dat$logp <- -log(dat$pvalue ,10)
  dat_w <- dat %>%
    select(list1, list2, logp) %>%
    pivot_wider(names_from = list2, values_from = logp) %>% as.data.frame()
  dat_w$list1 <- gsub(".*Hyper_|.*Hypo_|_genelist","",dat_w$list1) %>% as.factor()
  tissue_order <-paste0(c("Brain", "Breast", "Colon", "Kidney", "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Stomach"),"-Tumour")
  dat_w$list1 <- factor(dat_w$list1,levels=tissue_order)
  dat_w <-dat_w[order(dat_w$list1),]
  colnames(dat_w) <- gsub(" $","",colnames(dat_w))
  dat_w <- dat_w %>%
    select("list1", "Brain","Breast","Colon", "Kidney", "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Stomach")
  rownames(dat_w) <- dat_w$list1;dat_w$list1<-NULL
  pheatmap(
    dat_w %>% t(),
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
    breaks = seq(0,5,5/100), fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
    main = paste0(dmb,"logP"),  
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.logP.pdf"),  width = 4.5, height = 4.5,
  )
  pheatmap(
    dat_w %>% t(),
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
    fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
    main = paste0(dmb,"logP"),  
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.logP.scale.pdf"),  width = 4.5, height = 4.5,scale="row"
  )
  dat_w <- dat %>%
    select(list1, list2, odds_ratio) %>%
    pivot_wider(names_from = list2, values_from = odds_ratio) %>% as.data.frame()
  dat_w$list1 <- gsub(".*Hyper_|.*Hypo_|_genelist","",dat_w$list1) %>% as.factor()
  tissue_order <-paste0(c("Brain", "Breast", "Colon", "Kidney", "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Stomach"),"-Tumour")
  dat_w$list1 <- factor(dat_w$list1,levels=tissue_order)
  dat_w <-dat_w[order(dat_w$list1),]
  colnames(dat_w) <- gsub(" $","",colnames(dat_w))
  dat_w <- dat_w %>%
    select("list1", "Brain","Breast","Colon", "Kidney", "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Stomach")
  rownames(dat_w) <- dat_w$list1;dat_w$list1<-NULL
  pheatmap(
    dat_w %>% t(),
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
    breaks = seq(0,5,5/100), fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
    main = paste0(dmb,"OR"), 
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.OR.pdf"),  width = 4.5, height = 4.5,
  )
}

##### tissue-specific gene #####
genelists <- grep("Tumour",grep("capsHyper|tapsbetaHypo",list.files(pattern="*genelist$", path = genelist_path),value=TRUE),invert = TRUE,value = TRUE)
result_list <- lapply(genelists, function(genelist) {
  ts_gene_enrich(
    target_gene = paste0(genelist_path,genelist),
    tissue_specific_gene = paste0(resource_path,"PanglaoDB_immune_GTEx_Tissues_specific_genelist.MANE.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed")
  )
})

result <- do.call(rbind, result_list)
result <- data.frame(result)
result$list1 <- gsub("Erythroblasts","Erythroid-precursors",result$list1)
write.table(result,"PanglaoDB_immune_GTEx_Tissues_specific_gene.dmb_enrichment.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)


for(dmb in c("capsHyper_","tapsbetaHypo_")){
  dat <-  result[grep(dmb,result$list1),]
  dat$logp <- -log(dat$pvalue ,10)
  dat_w <- dat %>%
    select(list1, list2, logp) %>%
    pivot_wider(names_from = list2, values_from = logp) %>% as.data.frame()
  dat_w$list1 <- gsub("capsHyper_|tapsbetaHypo_|_genelist","",dat_w$list1) %>% as.factor()
  tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
  dat_w <- dat_w[dat_w$list1%in%tissue_order,]
  dat_w$list1 <- factor(dat_w$list1,levels=tissue_order)
  dat_w <-dat_w[order(dat_w$list1),]
  colnames(dat_w) <- gsub(" $","",colnames(dat_w))
  dat_w <- dat_w %>%
    select("list1", "Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", "T Helper Cells", "T Cytotoxic Cells", "NK Cells", "B Cells", "Neutrophils", "Eosinophils", "Monocytes",  "Erythroid Precursor Cells","Megakaryocytes")
  rownames(dat_w) <- dat_w$list1;dat_w$list1<-NULL
  pheatmap(
    dat_w %>% t(),
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
    breaks = seq(0,5,5/100), fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
    filename = paste0(dmb,"DMB.PanglaoDB_immune_GTEx_Tissues.logP.pdf"),  width = 6, height = 4.5,
  )
  pheatmap(
    dat_w %>% t(),
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
    fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
    filename = paste0(dmb,"DMB.PanglaoDB_immune_GTEx_Tissues.logP.scale.pdf"),  width = 6, height = 4.5,scale="row"
  )
}

##### feature enrichment on tissue specific genes #####
setwd(WORKDIR)

enrich_res <- fread("genelist/PanglaoDB_immune_GTEx_Tissues_specific_gene.dmb_enrichment.txt")
# enrich_res[!grepl("Cirrhosis|Pancreatitis", enrich_res$V1),]$V1 %>%
#   lapply(.,function(x) unlist(strsplit(x, "_"))[[2]]) %>%
#   unlist() %>%
#   unique() %>% sort()
# unique(enrich_res$V2) %>% sort()
caps_dmb_gene <- fread("CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv.all_genelist", col.names =c("chr","start","end","type","cell_type","g_chr","g_start","g_end","gene","dis") )


enrich_res$V2 <- gsub("B Cells","B-cells", enrich_res$V2)
enrich_res$V2 <- gsub("T Helper Cells","CD4-T-cells", enrich_res$V2)
enrich_res$V2 <- gsub("T Cytotoxic Cells","CD8-T-cells", enrich_res$V2)
enrich_res$V2 <- gsub("NK Cells","NK-cells", enrich_res$V2)
enrich_res$V2 <- gsub("Erythroid Precursor Cells","Erythroid-precursors", enrich_res$V2)
enrich_res$V2 <- gsub("Erythroid Precursor Cells","Erythroid-precursors", enrich_res$V2)

enrich_res %>%
  filter(V1 == paste0("capsHyper_", V2, "_genelist")) %>%
  select(V2,V3)%>%
  separate_rows(V3, sep = ";")  %>%
  rename(cell_type = V2,
         gene = V3) %>%
  merge(., caps_dmb_gene, by=c("cell_type","gene")) %>%
  filter(type=="Hyper")%>%
  select(chr,start,end,type,cell_type) %>%
  write.table(.,"CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.tissue_specific_gene.bed", quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE) 
  
tapsbeta_dmb_gene <- fread("TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv.all_genelist", col.names =c("chr","start","end","type","cell_type","g_chr","g_start","g_end","gene","dis") )

enrich_res %>%
  filter(V1 == paste0("tapsbetaHypo_", V2, "_genelist")) %>%
  select(V2,V3)%>%
  separate_rows(V3, sep = ";")  %>%
  rename(cell_type = V2,
         gene = V3) %>%
  merge(., tapsbeta_dmb_gene, by=c("cell_type","gene")) %>%
  filter(type=="Hypo")%>%
  select(chr,start,end,type,cell_type) %>%
  write.table(.,"TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.tissue_specific_gene.bed", quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE) 

#### 3. ATAC-seq ####
```
mkdir -p atac_seq/normal atac_seq/tumour
cd atac_seq/normal
sbatch /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/code_tissue_methylome_atlas/dmb_downstream_analysis/atac_seq_normal.sh
cd -;cd atac_seq/tumour
sbatch /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/code_tissue_methylome_atlas/dmb_downstream_analysis/atac_seq_tumour.sh
```


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


##### for tumour #####
setwd(WORKDIR)
tapsb_dmb <- fread("TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv")
caps_dmb <- fread("CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv")
setwd("atac_seq/tumour/")

pdf("atac_seq_tumour_dmb.pdf",width = 6,height = 5)
dmrs <- c("TAPSbeta.Hypo.tumour","TAPSbeta.Hyper.tumour","CAPS.Hypo.tumour","CAPS.Hyper.tumour")
for(dmr in dmrs){
  list.files(pattern=paste0(dmr,".*"))
  
  dat <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(dat) <- c("chr", "start", "end", "score","smp")
  
  for (i in list.files(pattern=paste0(dmr,".*inser.*"))){
    temp <- fread(i)
    temp$id <- gsub(".insertions.bed","",i)
    dat <- rbind(dat,temp)
  }
  dat <- as.data.frame(dat)
  dat$cancer <- gsub("_.*","",dat$id)
  colnames(dat)<- c("chr", "start", "end","dmb","selected_tissue", "score","smp","cancer")
  dat_sum <- dat %>%
    select(chr,start,end,score,cancer) %>%
    group_by(chr,start,end,cancer) %>%
    summarize(mean_value = mean(score, na.rm = TRUE)) %>%
    pivot_wider(names_from = "cancer", values_from = "mean_value") %>%
    as.data.frame()
  
  dat_sum$chr_start_end <- paste(dat_sum$chr,dat_sum$start,dat_sum$end,sep="_")
  if(dmr %in% c("TAPSbeta.Hypo.tumour","TAPSbeta.Hyper.tumour")){
    dat_sum <- tapsb_dmb %>% select(chr_start_end,Hyper_or_Hypo,selected_tissue) %>%
      filter(grepl("Tumour", selected_tissue)) %>%
      merge(., dat_sum[,-c(1:3)],by=c("chr_start_end"))
  }else{
    dat_sum <- caps_dmb %>% select(chr_start_end,Hyper_or_Hypo,selected_tissue) %>%
      filter(grepl("Tumour", selected_tissue)) %>%
      merge(., dat_sum[,-c(1:3)],by=c("chr_start_end"))
  }
  
  colnames(dat_sum)<- gsub(paste0(dmr,"."),"",colnames(dat_sum))
  dat_sum_sel <- dat_sum[order(dat_sum$selected_tissue),] %>%
    filter(grepl("Brain|Breast|Colon|Kidney|Liver|Lung|Prostate|Stomach", selected_tissue))%>%
    filter(grepl(strsplit(dmr,"\\.")[[1]][2], Hyper_or_Hypo))%>%
    select(chr_start_end, selected_tissue, GBMx, BRCA, COAD, KIRC, LIHC, LUAD, PRAD, STAD) %>%
    as.data.frame()
  rownames(dat_sum_sel) <- paste0(dat_sum_sel$chr_start_end,dat_sum_sel$selected_tissue)
  
  row_annotations <- data.frame(cancer=dat_sum_sel$selected_tissue)
  rownames(row_annotations) <- paste0(dat_sum_sel$chr_start_end,dat_sum_sel$selected_tissue)
  
  pheatmap(dat_sum_sel[,-c(1:2)],
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           show_rownames = FALSE,
           annotation_row = row_annotations,
           main = dmr)
  
  pheatmap(dat_sum_sel[,-c(1:2)],
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           show_rownames = FALSE,
           scale="row",
           annotation_row = row_annotations,
           main = dmr)
}
dev.off()


# "Ovary-Tumour" "Pancreas-Tumour" are not in https://www.ncbi.nlm.nih.gov/pubmed/30361341


##### for normal #####
setwd(WORKDIR)
tapsb_dmb <- fread("TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv")
caps_dmb <- fread("CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv")
setwd("atac_seq/normal/")


pdf("atac_seq_normal_dmb.pdf",width = 6,height = 5)
dmrs <- c("TAPSbeta.Hypo.normal","TAPSbeta.Hyper.normal","CAPS.Hypo.normal","CAPS.Hyper.normal")

for(dmr in dmrs){
  list.files(pattern=paste0(dmr,".*"))
  
  dat <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(dat) <- c("chr", "start", "end", "score","smp")
  
  for (i in list.files(pattern=paste0(dmr,".*EN.*"))){
    temp <- fread(i)
    temp$id <- gsub(".insertions.bed","",i)
    dat <- rbind(dat,temp)
  }
  dat <- as.data.frame(dat)
  dat$tissue <- lapply(dat$id,function(x)strsplit(x,"\\.") %>% unlist() %>% .[5]) %>% unlist()
  colnames(dat)<- c("chr", "start", "end","dmb","selected_tissue", "score","smp","tissue")
  
  
  dat_sum <- dat %>%
    select(chr,start,end,score,selected_tissue,tissue) %>%
    group_by(chr,start,end,selected_tissue,tissue) %>%
    summarize(mean_value = mean(score, na.rm = TRUE)) %>%
    pivot_wider(names_from = "tissue", values_from = "mean_value") %>%
    as.data.frame()
  
  dat_sum$chr_start_end <- paste(dat_sum$chr,dat_sum$start,dat_sum$end,sep="_")
  if(dmr %in% c("TAPSbeta.Hypo.normal","TAPSbeta.Hyper.normal")){
    dat_sum <- tapsb_dmb %>% filter(!grepl("Pancreatitis|Cirrhosis",selected_tissue)) %>% select(chr_start_end,Hyper_or_Hypo,selected_tissue) %>%
      merge(., dat_sum[,-c(1:3)],by=c("chr_start_end","selected_tissue"))
  }else{
    dat_sum <- caps_dmb  %>% filter(!grepl("Pancreatitis|Cirrhosis",selected_tissue)) %>%select(chr_start_end,Hyper_or_Hypo,selected_tissue) %>%
      merge(., dat_sum[,-c(1:3)],by=c("chr_start_end","selected_tissue"))
  }
  colnames(dat_sum) <- gsub(",","", colnames(dat_sum))
  tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells", "CD8-T-cells", "Neutrophils", "NK-cells", "B-cells", "Eosinophils", "Monocytes", "Megakaryocytes","Erythroblasts")
  
  dat_sum$selected_tissue <- factor(dat_sum$selected_tissue, levels = tissue_order)
  dat_sum <- dat_sum[order(dat_sum$selected_tissue),]
  # [1] "chr_start_end"                   "Hyper_or_Hypo"                  
  # [3] "selected_tissue"                 "body_of_pancreas"               
  # [9] "esophagus_mucosa"                "esophagus_muscularis_mucosa"    
  # [11] "esophagus_squamous_epithelium"  "heart_right_ventricle"                          
  
  dat_sum_sel <- dat_sum %>%
    filter(!grepl("Neutrophils|Megakaryocytes|Erythroblasts|Eosinophils|Monocytes", selected_tissue)) %>%
    select(chr_start_end, selected_tissue, cerebellum, breast_epithelium, heart_left_ventricle, kidney, liver, lung, ovary, pancreas, prostate_gland, transverse_colon, stomach, esophagus_mucosa, spleen, `CD4-positive_alpha-beta_T_cell`, `CD8-positive_alpha-beta_T_cell`, natural_killer_cell, naive_B_cell) %>%
    as.data.frame()
  
  rownames(dat_sum_sel) <- paste0(dat_sum_sel$chr_start_end,dat_sum_sel$selected_tissue)
  
  row_annotations <- data.frame(tissue=dat_sum_sel$selected_tissue)
  rownames(row_annotations) <- paste0(dat_sum_sel$chr_start_end,dat_sum_sel$selected_tissue)
  
  print(pheatmap(dat_sum_sel[,-c(1:2)],
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           show_rownames = FALSE,
           annotation_row = row_annotations,
           main = dmr)[[4]])
  
  print(pheatmap(dat_sum_sel[,-c(1:2)],
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           show_rownames = FALSE,
           scale="row",
           annotation_row = row_annotations,
           main=dmr)[[4]])
}
dev.off()

#### 4. DMB number ####
WORKDIR <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr15/"
setwd(WORKDIR)
dir.create("dmb_num")
for (infile  in c("TAPSbeta.all_samples.Tissue_Group_DMBs.all_delta_quants0.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
                  "CAPS.all_samples.Tissue_Group_DMBs.all_delta_quants0.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv")) {
  prefix <- gsub(".csv","",infile)
  dmb <- fread(paste0("raw/",infile))
  dmb$selected_tissue <- gsub("Erythroblasts","Erythroid-precursors",dmb$selected_tissue)
  tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                    "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
  seldmb <- dmb[dmb$selected_tissue%in%tissue_order,] 
  seldmb$selected_tissue <- factor(seldmb$selected_tissue,levels = tissue_order)
  seldmb <- seldmb[order(seldmb$selected_tissue),]
  p <- seldmb %>%
    select(delta_quants,significant_dmrs,Hyper_or_Hypo,selected_tissue) %>%
    ggplot(aes(x=selected_tissue,y=delta_quants,color=significant_dmrs)) +
    geom_jitter(width = 0.1,size=0.2)+
    facet_wrap(~Hyper_or_Hypo,nrow=2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom")
  
  ggsave(paste0("dmb_num/",prefix,".dot.png"),p,width = 10, height = 8)
  
  selN <- 300
  seldmb_sta <- seldmb %>%
    select(delta_quants,significant_dmrs,Hyper_or_Hypo,selected_tissue) %>%
    group_by(selected_tissue) %>%
    arrange(desc(delta_quants)) %>%
    slice_head(n = selN) %>%
    select(Hyper_or_Hypo,selected_tissue) %>% 
    table() %>% as.data.frame() 
  seldmb_sta$selected_tissue <- factor(seldmb_sta$selected_tissue,levels = rev(tissue_order))
  seldmb_sta <- seldmb_sta[order(seldmb_sta$selected_tissue),]
  p <- seldmb_sta %>%
    ggplot(aes(x=selected_tissue,y=Freq,fill=Hyper_or_Hypo)) +
    geom_bar(stat="identity")+
    theme_minimal() +
    scale_fill_manual(values=  brewer.pal(4,"RdBu")[c(1,4)])+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom") +
    coord_flip()
  ggsave(paste0("dmb_num/",prefix,".barplot.pdf"),p,width = 4, height = 5)
}

for (infile  in c("TAPSbeta.all_samples.Tissue_Group_DMBs.all_delta_quants0.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv",
                  "CAPS.all_samples.Tissue_Group_DMBs.all_delta_quants0.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv")) {
  prefix <- gsub(".csv","",infile)
  dmb <- fread(paste0("raw/",infile))
  dmb$selected_tissue <- gsub("Erythroblasts","Erythroid-precursors",dmb$selected_tissue)
  tissue_order <- c("Brain-Tumour","Breast-Tumour", "Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Colon-Tumour","Stomach-Tumour",
                    "Liver-Cirrhosis","Pancreas-Pancreatitis")
  seldmb <- dmb[dmb$selected_tissue%in%tissue_order,] 
  seldmb$selected_tissue <- factor(seldmb$selected_tissue,levels = tissue_order)
  seldmb <- seldmb[order(seldmb$selected_tissue),]
  p <- seldmb %>%
    select(delta_quants,significant_dmrs,Hyper_or_Hypo,selected_tissue) %>%
    ggplot(aes(x=selected_tissue,y=delta_quants,color=significant_dmrs)) +
    geom_jitter(width = 0.1,size=0.2)+
    facet_wrap(~Hyper_or_Hypo,nrow=2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom")
  
  ggsave(paste0("dmb_num/",prefix,".dot.tumour.png"),p,width = 10, height = 8)
  
  selN <- 300
  seldmb_sta <- seldmb %>%
    select(delta_quants,significant_dmrs,Hyper_or_Hypo,selected_tissue) %>%
    group_by(selected_tissue) %>%
    arrange(desc(delta_quants)) %>%
    slice_head(n = selN) %>%
    select(Hyper_or_Hypo,selected_tissue) %>% 
    table() %>% as.data.frame() 
  seldmb_sta$selected_tissue <- factor(seldmb_sta$selected_tissue,levels = rev(tissue_order))
  seldmb_sta <- seldmb_sta[order(seldmb_sta$selected_tissue),]
  p <- seldmb_sta %>%
    ggplot(aes(x=selected_tissue,y=Freq,fill=Hyper_or_Hypo)) +
    geom_bar(stat="identity")+
    theme_minimal() +
    scale_fill_manual(values=  brewer.pal(4,"RdBu")[c(1,4)])+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom") +
    coord_flip()
  ggsave(paste0("dmb_num/",prefix,".barplot.tumour.pdf"),p,width = 4, height = 5)
}
