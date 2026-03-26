.libPaths(c(.libPaths(),"/gpfs3/users/ludwig/cfo155/R/x86_64-pc-linux-gnu-library/4.2"))
library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)
data(XenaData)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/gene_expr")

#### download data from toil ####
tumourInfo <- read.table("tumour_cohort.txt", header=FALSE)
all_tcga <- XenaGenerate(subset = XenaHostNames == "tcgaHub") 

options(warn=-1)
download_tumour_data <- function(tissue, paraHistologicalType){
  all_tcga <- XenaGenerate(subset = XenaHostNames == "tcgaHub")
  cohort <- tumourInfo$V1[tumourInfo$V2==tissue]
  paraCohort = gsub(" \\(.*", "",grep(paste0("\\(",cohort,"\\)"),all_tcga@cohorts,value=TRUE)) 
  paraDatasets = grep(paste0(cohort,"_clini"),all_tcga@datasets,value=TRUE) 
  if(paraHistologicalType=="NA"){
    paraHistologicalType = paste0(tissue, " Adenocarcinoma"); 
  }
  
  GeneExpectedCnt_toil = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
    XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
    XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count");
  XenaQuery(GeneExpectedCnt_toil) %>%
    XenaDownload(destdir = "./")
  
  Clin_TCGA = XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
    XenaFilter(filterCohorts = paraCohort) %>%
    XenaFilter(filterDatasets = paraDatasets);
  XenaQuery(Clin_TCGA) %>%
    XenaDownload(destdir = "./")
  Pheno_GTEx = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
    XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
    XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype");
  XenaQuery(Pheno_GTEx) %>%
    XenaDownload(destdir = "./")
  
  filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz");
  names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01));
  paraStudy = "GTEX"; #Setting "GTEx" as the study of interest.
  paraPrimarySiteGTEx = tissue; 
  paraPrimaryTissueGTEx = paste0("^",tissue); 
  filterGTEx02 = subset(filterGTEx01,
                        study == paraStudy &
                          primarysite == paraPrimarySiteGTEx &
                          grepl(paraPrimaryTissueGTEx, filterGTEx01$`primary disease or tissue`))
  filterTCGA01 = fread(paraDatasets);
  names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01));
  paraSampleType = "Primary Tumor"; #Setting "Primary Tumor" as the sample type of interest.
  paraPrimarySiteTCGA = tissue; 
  print(table(filterTCGA01$histologicaltype))
  
  filterTCGA02 = subset(filterTCGA01,
                        sampletype == paraSampleType &
                          primarysite == paraPrimarySiteTCGA &
                          grepl(paraHistologicalType, filterTCGA01$histologicaltype))
  
  filterExpr = c(filterGTEx02$sample, filterTCGA02$sampleID, "sample");
  
  ExprSubsetBySamp = fread("TcgaTargetGtex_gene_expected_count.gz",
                           select = filterExpr)
  suffix1 <- paste0(tissue,"_gtex_n",length(grep("GTEX",colnames(ExprSubsetBySamp))),"_mean")
  suffix2 <- paste0(tissue,"_tcga_n",length(grep("TCGA",colnames(ExprSubsetBySamp))),"_mean")
  print(suffix1)
  print(suffix2)
  ExprSubsetBySamp[, suffix1] <- ExprSubsetBySamp %>%
    select(contains("GTEX")) %>%
    rowMeans()
  ExprSubsetBySamp[, suffix2] <- ExprSubsetBySamp %>%
    select(contains("TCGA")) %>%
    rowMeans()
  rownames(ExprSubsetBySamp) <- ExprSubsetBySamp$sample
  ExprSubsetBySamp$sample <- NULL
  write.csv(ExprSubsetBySamp%>%select(contains("_mean")), paste0(tissue,"_ExpectedCnt_mean.csv"))
  write.csv(ExprSubsetBySamp%>%select(!contains("_mean")), paste0(tissue,"_ExpectedCnt.csv"))
  
  ExprSubsetBySamp = fread("TcgaTargetGtex_rsem_gene_tpm.gz",
                           select = filterExpr)
  suffix1 <- paste0(tissue,"_gtex_n",length(grep("GTEX",colnames(ExprSubsetBySamp))),"_mean")
  suffix2 <- paste0(tissue,"_tcga_n",length(grep("TCGA",colnames(ExprSubsetBySamp))),"_mean")
  ExprSubsetBySamp[, suffix1] <- ExprSubsetBySamp %>%
    select(contains("GTEX")) %>%
    rowMeans()
  ExprSubsetBySamp[, suffix2] <- ExprSubsetBySamp %>%
    select(contains("TCGA")) %>%
    rowMeans()
  rownames(ExprSubsetBySamp) <- ExprSubsetBySamp$sample
  ExprSubsetBySamp$sample <- NULL
  write.csv(ExprSubsetBySamp%>%select(contains("_mean")), paste0(tissue,"_rsem_gene_tpm_mean.csv"))
  write.csv(ExprSubsetBySamp%>%select(!contains("_mean")), paste0(tissue,"_rsem_gene_tpm.csv"))
}

tumourInfo <- read.table("tumour_cohort.txt", header=FALSE)
for(tissue in tumourInfo$V2){
  dat <- fread(paste0(tissue,"_rsem_gene_tpm.csv"))
  tmp <- lapply(dat[, -1, with = FALSE], function(x) 2^x - 0.001) %>%
    as.data.frame()
  tmp <- lapply(tmp, function(x) ifelse(x < 0, 0, x)) %>%
    as.data.frame()
  tmp <- cbind(dat[,1],tmp)
  
  suffix1 <- paste0(tissue,"_gtex_n",length(grep("GTEX",colnames(tmp))),"_mean")
  suffix2 <- paste0(tissue,"_tcga_n",length(grep("TCGA",colnames(tmp))),"_mean")
  # suffix1 <- paste0(tissue,"_gtex_n",length(grep("GTEX",colnames(tmp))),"_median")
  # suffix2 <- paste0(tissue,"_tcga_n",length(grep("TCGA",colnames(tmp))),"_median")
  print(suffix1)
  print(suffix2)
  tmp[, suffix1] <- tmp %>%
    select(contains("GTEX")) %>%
    # rowMeans()
    apply(1,median)
  tmp[, suffix2] <- tmp %>%
    select(contains("TCGA")) %>%
    # rowMeans()
    apply(1,median)
  rownames(tmp) <- tmp$V1
  tmp$V1 <- NULL
  # write.csv(tmp%>%select(!contains("_mean")), paste0(tissue,"_rsem_gene_tpm.raw.csv"))
  write.csv(tmp%>%select(contains("_mean")), paste0(tissue,"_rsem_gene_tpm_mean.raw.csv"))
  # write.csv(tmp%>%select(contains("_median")), paste0(tissue,"_rsem_gene_tpm_median.raw.csv"))
}


tpm <- fread("all.rsem_gene_tpm_mean.raw.csv")
probmap <- fread("probeMap%2Fgencode.v23.annotation.gene.probemap")
tpm <- merge(probmap[,1:2],tpm,by.x=c("id"),by.y=c("V1"))
write.csv(tpm, "all.rsem_gene_tpm_mean.raw.genename.csv", quote=FALSE, row.names=FALSE)

tpm <- fread("all.rsem_gene_tpm_median.raw.csv")
probmap <- fread("probeMap%2Fgencode.v23.annotation.gene.probemap")
tpm <- merge(probmap[,1:2],tpm,by.x=c("id"),by.y=c("V1"))
write.csv(tpm, "all.rsem_gene_tpm_median.raw.genename.csv", quote=FALSE, row.names=FALSE)


#### compare toil with HPA ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/gene_expr")
toil_tpm <- fread("all.rsem_gene_tpm_mean.raw.genename.csv")
# toil_tpm <- fread("all.rsem_gene_tpm_median.raw.genename.csv")
# toil_tpm_log <- lapply(toil_tpm[, -c(1:2), with = FALSE], function(x)log2(x+0.001)) %>% as.data.frame()
# toil_tpm_log <- cbind(toil_tpm[,c(1:2)], toil_tpm_log)

gtex_tpm <- fread("rna_tissue_gtex.tsv")
# gtex_tpm$TPM <- log2(gtex_tpm$TPM+0.001)
gtex_tpm <- gtex_tpm %>%
  pivot_wider(
    id_cols = c("Gene", "Gene name"), 
    names_from = "Tissue",
    values_from = c("TPM")
  )
immune_tpm <- fread("rna_immune_cell.tsv")
# immune_tpm$TPM <- log2(immune_tpm$TPM+0.001)
immune_tpm <- immune_tpm %>%
  pivot_wider(
    id_cols = c("Gene", "Gene name"), 
    names_from = "Immune cell",
    values_from = c("TPM")
  )

gtex_tissue <- c("breast","colon","kidney","liver","lung","ovary","pancreas","prostate","stomach")
# all_tpm <- merge(toil_tpm_log %>%select(contains("gene") | contains("gtex")),
#                  gtex_tpm[,-1] %>% select(c("Gene name",gtex_tissue)),by.x="gene",by.y="Gene name")
all_tpm <- merge(toil_tpm %>%select(contains("gene") | contains("gtex")),
                 gtex_tpm[,-1] %>% select(c("Gene name",gtex_tissue)),by.x="gene",by.y="Gene name")
all_tpm <- all_tpm[complete.cases(all_tpm),]
all_tpm <- all_tpm[!grep("^MT-",all_tpm$gene),]
all_tpm_cor <- cor(all_tpm[,-c(1:2)]) %>% round(2)
all_tpm_cor <- all_tpm_cor[order(rownames(all_tpm_cor)),order(rownames(all_tpm_cor))]
pheatmap(all_tpm_cor,display_numbers = TRUE,cluster_rows = FALSE, cluster_cols = FALSE,main ="correlation between toil_gtex_median(raw) and hpa_gtex")
pheatmap(all_tpm_cor,display_numbers = TRUE,cluster_rows = FALSE, cluster_cols = FALSE,main ="correlation between toil_gtex_mean(raw) and hpa_gtex")
pheatmap(all_tpm_cor,display_numbers = TRUE,cluster_rows = FALSE, cluster_cols = FALSE,main ="correlation between toil_gtex_mean(log2) and hpa_gtex")


ggplot(all_tpm,aes(x=Kidney_gtex_n28_mean,y=kidney)) + geom_point() + xlim(0,1000) + ylim(0,1000) +
ggplot(all_tpm,aes(x=Colon_gtex_n308_mean,y=colon)) + geom_point() + xlim(0,1000) + ylim(0,1000) +
ggplot(all_tpm,aes(x=Liver_gtex_n110_mean,y=liver)) + geom_point() + xlim(0,1000) + ylim(0,1000)
order_col <- c(order(colnames(all_tpm)[-1]) + 1)
png("gtex_toil_hpa.png", width = 20, height = 12, units = "in", res = 150)
par(mar = c(4, 4, 1, 1))  # Adjust the margins as needed
print(as.data.frame(all_tpm)[, order_col] %>%
  pairs(method="lm"),xlim=c(0,1000), ylim = c(0,1000), lower.panel=NULL)
dev.off()


#### tumour specific gene from toil ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/gene_expr")
toil_tpm <- fread("all.rsem_gene_tpm_mean.raw.genename.csv")
# gtex_tpm <- fread("rna_tissue_gtex.tsv")
all_gene <- read.table("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/resource/MANE.GRCh38.v1.0.refseq_genomic.gene.bed",
                       col.names = c("chr","start","end","gene","score","strand"))
toil_tpm <- toil_tpm[toil_tpm$gene%in%all_gene$gene,]
tumourInfo <- read.table("tumour_cohort.txt", header=FALSE)
selgene_up <- data.frame()
selgene_down <- data.frame()
for(tissue in tumourInfo$V2){
  fc_cutoff <- 1.5
  diff_cutoff <- 1 # 0
  tmp <- cbind(toil_tpm %>% select("gene"), 
               toil_tpm %>% select(contains(paste0(tissue,"_tcga"))),
               toil_tpm %>% select(contains("mean") & !contains(paste0(tissue, "_tcga")))%>% apply(1,max),
               toil_tpm %>% select(contains("mean") & !contains(paste0(tissue, "_tcga")))%>% apply(1,min))
  tmp$fc1 <- tmp[,2]/tmp[,3]
  tmp$fc2 <- tmp[,4]/tmp[,2]
  # write.csv(tmp$gene[which(tmp$fc>=fc_cutoff)], paste0(tissue,"_cancer_gene.csv"),row.names = FALSE, quote = FALSE)
  selgene_up <- rbind(selgene_up,data.frame("tissue"=tissue,"gene"=tmp$gene[which(tmp$fc1>=fc_cutoff & tmp[,2] - tmp[,3]>diff_cutoff)]))
  if(length(which(tmp$fc2>=fc_cutoff & tmp[,4] - tmp[,2]>diff_cutoff))>0){
    selgene_down <- rbind(selgene_down,data.frame("tissue"=tissue,"gene"=tmp$gene[which(tmp$fc2>=fc_cutoff & tmp[,4] - tmp[,2]>diff_cutoff)]))
  }
}
selgene_up_tpm <- merge(selgene_up, toil_tpm %>% select(-id), by="gene")
selgene_up_tpm <- selgene_up_tpm[order(selgene_up_tpm$tissue),]
selgene_down_tpm <- merge(selgene_down, toil_tpm %>% select(-id), by="gene")
selgene_down_tpm <- selgene_down_tpm[order(selgene_down_tpm$tissue),]
pdf("heatmap_for_tumour_specific_genes.raw.diff1.pdf",width = 10, height = 5)
pheatmap(selgene_up_tpm %>% select(contains("mean")),
         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,breaks = seq(0,100,1), color=colorRampPalette(c("navy", "white", "red"))(50))
pheatmap(selgene_down_tpm %>% select(contains("mean")),
         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,breaks = seq(0,100,1), color=colorRampPalette(c("navy", "white", "red"))(50))
merge(as.data.frame(table(selgene_up_tpm$tissue)),
      as.data.frame(table(selgene_down_tpm$tissue)),by="Var1",all.x=TRUE, all.y=TRUE) %>%
  melt(id.vars=("Var1")) %>%
  ggplot(aes(x=Var1,y=value,fill=variable)) +
  geom_bar(stat="identity",position="dodge") +
  geom_text(aes(label=value,position="dodge"), vjust=0) +
  scale_fill_discrete(breaks=c("Freq.x","Freq.y"),labels=c("Up","Down"))+
  theme_bw()
dev.off()



# write.table(selgene[,c(2,1)],"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/toil_tumour_specific_gene.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)
write.table(selgene_up[,c(2,1)],"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/toil_tumour_specific_gene.diff1.up.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)
write.table(selgene_down[,c(2,1)],"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/toil_tumour_specific_gene.diff1.down.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)
write.table(selgene_up[,c(2,1)],"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/toil_tumour_specific_gene.diff0.up.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)
write.table(selgene_down[,c(2,1)],"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/toil_tumour_specific_gene.diff0.down.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)

#### tumour dmb enrich ####
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

resource_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/"
###### list from tissue_atlas_v3_dmr10 ######
genelist_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr10/genelist/"
setwd(genelist_path)
genelists <- grep("Tumour",grep("capsHyper|tapsbetaHypo",list.files(pattern="*genelist$", path = genelist_path),value=TRUE),value = TRUE)
result_list <- lapply(genelists, function(genelist) {
  ts_gene_enrich(
    target_gene = paste0(genelist_path,genelist),
    tissue_specific_gene = paste0(resource_path,"toil_tumour_specific_gene.up.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed")
  )
})


result <- do.call(rbind, result_list)
result <- data.frame(result)

for(dmb in c("capsHyper","tapsbetaHypo")){
  dat <-  result[grep(dmb,result$list1),]
  dat$logp <- -log(dat$pvalue ,10)
  dat_w <- dat %>%
    select(list1, list2, logp) %>%
    pivot_wider(names_from = list2, values_from = logp) %>% as.data.frame()
  dat_w$list1 <- gsub("capsHyper_|tapsbetaHypo_|_genelist","",dat_w$list1) %>% as.factor()
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
    main = paste0(dmb,"logP"),  
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.logP.png"),  width = 4.5, height = 4.5,
  )
  dat_w <- dat %>%
    select(list1, list2, odds_ratio) %>%
    pivot_wider(names_from = list2, values_from = odds_ratio) %>% as.data.frame()
  dat_w$list1 <- gsub("capsHyper_|tapsbetaHypo_|_genelist","",dat_w$list1) %>% as.factor()
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
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.OR.png"),  width = 4.5, height = 4.5,
  )
}

write.table(result,"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr10/genelist/toil_tumour_specific_gene.dmb_enrichment.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)
###### list from Masato (toil)######
resource_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/"
genelist_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/dmb_tumour_mi/genelist/"
setwd(genelist_path)
genelists <- list.files(pattern="*genelist$", path = genelist_path)
                  
result_list <- lapply(genelists, function(genelist) {
  ts_gene_enrich(
    target_gene = paste0(genelist_path,genelist),
    tissue_specific_gene = paste0(resource_path,"toil_tumour_specific_gene.up.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed")
  )
})


result <- do.call(rbind, result_list)
result <- data.frame(result)

for(dmb in unique(gsub("Hyper.*","Hyper",result$list1) %>% gsub("Hypo.*","Hypo",.))){
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
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
    breaks = seq(0,5,5/100), fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
    main = paste0(dmb,"logP"),  
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.logP.png"),  width = 4.5, height = 4.5,
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
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.OR.png"),  width = 4.5, height = 4.5,
  )
}

###### list from tissue_atlas_v3_dmr13  (toil)######
resource_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/"
genelist_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr13/genelist/"
setwd(genelist_path)
genelists <- list.files(pattern="*Tumour.*genelist$", path = genelist_path)

result_list <- lapply(genelists, function(genelist) {
  ts_gene_enrich(
    target_gene = paste0(genelist_path,genelist),
    tissue_specific_gene = paste0(resource_path,"toil_tumour_specific_gene.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed")
  )
})


result <- do.call(rbind, result_list)
result <- data.frame(result)

for(dmb in unique(gsub("Hyper.*","Hyper",result$list1) %>% gsub("Hypo.*","Hypo",.))){
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
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
    breaks = seq(0,5,5/100), fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
    main = paste0(dmb,"logP"),  
    filename = paste0(dmb,"DMB.toil_tumour_specific_gene.logP.pdf"),  width = 4.5, height = 4.5,
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

###### list from Masato(enrichr)######
setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/dmb_tumour_mi/genelist")
dat <- fread("Cancer_Cell_Line_Encyclopedia.out")
dat$logp <- -log10(dat$P.value)
for(i in gsub("Hyper.*","Hyper",unique(dat$genelist)) %>% gsub("Hypo.*","Hypo",.) %>% unique()){
  seldat <- dat[grep(i,dat$genelist),]
  selterm <- seldat %>% group_by(genelist) %>%
    top_n(1, logp) %>% 
    select(Term)
  
  seldat_selterm <- seldat[seldat$Term%in%selterm$Term,] %>%
    select(genelist,Term, logp) %>% 
    pivot_wider(names_from = Term, values_from = logp) %>% as.data.frame()
  seldat_selterm <- seldat_selterm[,order(factor(colnames(seldat_selterm),levels=selterm$Term))]
  seldat_selterm[is.na(seldat_selterm)]<-0
  rownames(seldat_selterm)<- seldat_selterm$genelist
  seldat_selterm$genelist<- NULL
  pheatmap(seldat_selterm,cluster_cols = FALSE, cluster_rows = FALSE, filename = paste0(i,".top1_enriched.pdf"),width = 12, height = 8)
  
}

#### from HPA ####
# download from this science paper Table S2: https://www.science.org/doi/10.1126/science.aan2507
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/gene_expr")
dat <- fread("Science_2017_TableS2.txt")
tumourInfo <- read.table("tumour_cohort.txt", header=FALSE)
selgene <- data.frame()
for(tissue in tumourInfo$V2){
  selgene <- rbind(selgene,data.frame("tissue"=tissue,"gene"=dat$Symbols[grep(tumourInfo$V1[tumourInfo$V2==tissue],dat$`Category in cancer`)]))
}

write.table(selgene[,c(2,1)],"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/Science_2017_TableS2_tumour_specific_gene.txt", col.names =FALSE, quote = FALSE,sep="\t",row.names = FALSE)

resource_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/"
genelist_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr10/genelist/"
setwd(genelist_path)
genelists <- grep("Tumour",grep("capsHyper|tapsbetaHypo",list.files(pattern="*genelist$", path = genelist_path),value=TRUE),value = TRUE)
result_list <- lapply(genelists, function(genelist) {
  ts_gene_enrich(
    target_gene = paste0(genelist_path,genelist),
    tissue_specific_gene = paste0(resource_path,"Science_2017_TableS2_tumour_specific_gene.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed")
  )
})


result <- do.call(rbind, result_list)
result <- data.frame(result)

for(dmb in c("capsHyper","tapsbetaHypo")){
  dat <-  result[grep(dmb,result$list1),]
  dat$logp <- -log(dat$pvalue ,10)
  dat_w <- dat %>%
    select(list1, list2, logp) %>%
    pivot_wider(names_from = list2, values_from = logp) %>% as.data.frame()
  dat_w$list1 <- gsub("capsHyper_|tapsbetaHypo_|_genelist","",dat_w$list1) %>% as.factor()
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
    main = dmb, 
    filename = paste0(dmb,"DMB.Science_2017_TableS2_tumour_specific_gene.logP.png"),  width = 4.5, height = 4.5,
  )
}








seldat <- dat %>% filter(str_detect(`Category in cancer`, "Enriched in")) 
colnames(seldat)[2] <- c("gene")
seldat <- merge(seldat, toil_tpm, by=c("gene") )
seldat <- seldat[order(seldat$`Category in cancer`),]
tumour_types <- c("GBM","BRCA","COAD","KIRC","LIHC","LUAD","OV","PAAD","PRAD","STAD")
seldat <- seldat[grep("GBM|BRCA|COAD|KIRC|LIHC|LUAD|OV|PAAD|PRAD|STAD",seldat$`Category in cancer`),]
seldat$`Category in cancer` <- factor(seldat$`Category in cancer`,levels=paste0("Enriched in ",tumour_types))
seldat <- seldat[order(seldat$`Category in cancer`),]
pheatmap(seldat %>% select(contains("mean")),
         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, scale = "row")
pheatmap(seldat %>% select(contains("mean")),
         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,breaks = seq(0,100,1))




