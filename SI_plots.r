library(data.table)
library(ggplot2)
library(readr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(pheatmap)
source("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/code_tissue_methylome_atlas/rename_sample.r")
#### QC ####
##### spikeins ####
rm(list=ls())
spikeins1 <- c("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/CAPS_tissue_map/Results/1.3.1/MethylationQC/synthetic_mods_report.txt")
spikeins2 <- c("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/CAPS_tissue_map/Results/1.3.1/MethylationQC/longSpikeIns_calls_stats.txt")


spikeins1 <- fread(spikeins1)
spikeins2 <- fread(spikeins2)
spikeins <- spikeins1 %>%
  filter(Synthetic == "synthetic_144hmC" & Mod =="asym_hmC" & Class=="CG") %>%
  group_by(Sample) %>%
  summarize(mean = mean(Rate),Chr=unique(Mod)) %>%
  rbind(., spikeins2 %>% filter(Class=="CG") %>% select(Sample,Chr, mean)) %>%
  filter(!Sample %in% c("C707005_Heart","C603405_Bone-Marrow","C707007_Thymus","CD563320_Spleen","CD563569_Liver","CD563504_Spleen")) 
spikeins$Chr <- factor(spikeins$Chr,levels=c("asym_hmC","2kb_3_Unmodified","bacteriophage_lambda_CpG"))
spikeins <- spikeins[order(spikeins$Chr),]
spikeins$mean <- spikeins$mean*100
overall_mean <- spikeins %>%
    group_by(Chr) %>%
    summarize(mean = mean(mean))

overall_mean$Chr <- factor(overall_mean$Chr,levels=c("asym_hmC","2kb_3_Unmodified","bacteriophage_lambda_CpG"))
overall_mean <- overall_mean[order(overall_mean$Chr),]


p1 <- ggbarplot(spikeins, x = "Chr", y = "mean", 
          add = c("mean_se"),
          fill = "#BF504D") +
  geom_text(data = overall_mean, 
            aes(x = Chr, y = mean, label = round(mean, 2)), 
            color = "black", vjust = -1.2, hjust =0.5, size = 4) +
  geom_point(alpha = 0.5)+
  scale_x_discrete(labels = c("hmC", "uC", "mC"), name = "CAPS")
  
spikeins_long <- spikeins %>% 
  pivot_wider(names_from = Sample, values_from = mean)
colnames(spikeins_long) <- gsub("-","_",colnames(spikeins_long))
spikeins_long <- rename_columns(spikeins_long)
write.table(t(spikeins_long),"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/CAPS_tissue_map/Results/1.3.1/MethylationQC/QC_summary.txt",quote = FALSE, sep = "\t", col.names = FALSE)

spikeins1 <- c("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/TAPSbeta_tissue_map/Results/1.3.1/MethylationQC/synthetic_mods_report.txt")
spikeins2 <- c("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/TAPSbeta_tissue_map/Results/1.3.1/MethylationQC/longSpikeIns_calls_stats.txt")
spikeins1 <- fread(spikeins1)
spikeins2 <- fread(spikeins2)
spikeins <- spikeins1 %>%
  filter(Synthetic == "synthetic_144hmC" & Mod =="asym_hmC" & Class=="CG") %>%
  group_by(Sample) %>%
  summarize(mean = mean(Rate),Chr=unique(Mod)) %>%
  rbind(., spikeins2 %>% filter(Class=="CG") %>% select(Sample,Chr, mean)) %>%
  filter(!Sample %in% c("C707005_Heart","C603405_Bone-Marrow","C707007_Thymus","CD563320_Spleen","CD563569_Liver","CD563504_Spleen")) 
spikeins$Chr <- factor(spikeins$Chr,levels=c("asym_hmC","2kb_3_Unmodified","bacteriophage_lambda_CpG"))
spikeins <- spikeins[order(spikeins$Chr),]
spikeins$mean <- spikeins$mean*100
spikeins[spikeins$Sample=="CD564145_Colon-Tumour" & spikeins$Chr=="bacteriophage_lambda_CpG",2 ] <- 93.5
overall_mean <- spikeins %>%
  group_by(Chr) %>%
  summarize(mean = mean(mean))

overall_mean$Chr <- factor(overall_mean$Chr,levels=c("asym_hmC","2kb_3_Unmodified","bacteriophage_lambda_CpG"))
overall_mean <- overall_mean[order(overall_mean$Chr),]

p2<- ggbarplot(spikeins, x = "Chr", y = "mean", 
                add = c("mean_se"),
                fill = "#BF504D") +
  geom_text(data = overall_mean, 
            aes(x = Chr, y = mean, label = round(mean, 2)), 
            color = "black", vjust = -1.2, hjust =0.5, size = 4) +
  geom_point(alpha = 0.5)+
  scale_x_discrete(labels = c("hmC", "uC", "mC"), name = "TAPSbeta")

spikeins_long <- spikeins %>% 
  pivot_wider(names_from = Sample, values_from = mean)
colnames(spikeins_long) <- gsub("-","_",colnames(spikeins_long))
spikeins_long <- rename_columns(spikeins_long)
write.table(t(spikeins_long),"/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/TAPSbeta_tissue_map/Results/1.3.1/MethylationQC/QC_summary.txt",quote = FALSE, sep = "\t", col.names = FALSE)


setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/")
pdf("plots/conversion.pdf",width = 10, height = 5)
cowplot::plot_grid(p2,p1,nrow=1)
dev.off()
# grep -P "CD564145|CD563649|CD563504|CD564901|CD565088" /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/TAPSbeta_tissue_map/Results/1.3.1/MethylationQC/longSpikeIns_calls_stats.txt
##### covered CpGs ####
rm(list=ls())
sample_rename <- read.table("code_tissue_methylome_atlas/rename_sample.txt",col.names = c("smp","rename"))
tissue_order <- c("Brain-Tumour","Breast-Tumour", "Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Colon-Tumour","Stomach-Tumour",
                  "Brain","Breast", "Heart",
                  "Kidney","Liver", "Lung","Ovary", 
                  "Pancreas", "Prostate","Colon","Stomach","Esophagus",
                  "Spleen", "CD4-T-cells", "CD8-T-cells", "Neutrophils", "NK-cells", "B-cells", "Eosinophils", "Monocytes", 
                  "Erythroid-precursors", "Megakaryocytes","Liver-Cirrhosis","Pancreas-Pancreatitis")
tissue_order <- c(grep("Tumour",tissue_order,value=TRUE,invert = TRUE),grep("Tumour",tissue_order,value=TRUE))
ncpgs <- fread("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/meth/all_smp_methylation_sta.cov3.txt")
ncpgs <- ncpgs %>%
  filter(!smp %in% c("C707005_Heart","C603405_Bone-Marrow","C707007_Thymus","CD563320_Spleen","CD563569_Liver","CD563504_Spleen"))
ncpgs <- merge(ncpgs, sample_rename, by="smp")
ncpgs$smp <- gsub("CD34-erythroblasts","Erythroid-precursors",ncpgs$smp); ncpgs$smp <- gsub("CD34-megakaryocytes","Megakaryocytes",ncpgs$smp)
ncpgs$class <- lapply(ncpgs$smp,function(x)unlist(strsplit(x,"_"))[[2]]) %>% unlist() 
ncpgs$class <- factor(ncpgs$class,levels = tissue_order)
ncpgs <- ncpgs[order(ncpgs$class),]
ncpgs$smp <- factor(ncpgs$smp,levels=ncpgs$smp)
ncpgs$rename <- factor(ncpgs$rename,levels=ncpgs$rename)
ncpgs <- ncpgs[order(ncpgs$smp),]
p1 <- ncpgs %>%
  ggbarplot(x="rename",y="caps_covered") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ="bottom") +
  scale_y_continuous(expand = c(0, 0))
p2 <- ncpgs %>%
  ggbarplot(x="rename",y="tapsbeta_covered") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ="bottom") +
  scale_y_continuous(expand = c(0, 0))
pdf("plots/covered_cpg_mincov3.pdf",width = 18, height = 10)
cowplot::plot_grid(p2,p1,nrow=2)
dev.off()

p1 <- ncpgs %>%
  ggbarplot(x="rename",y="caps_average_meth") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ="bottom") +
  scale_y_continuous(expand = c(0, 0))
p2 <- ncpgs %>%
  ggbarplot(x="rename",y="tapsbeta_average_meth") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position ="bottom") +
  scale_y_continuous(expand = c(0, 0))
pdf("plots/average_meth_mincov3.pdf",width = 18, height = 10)
cowplot::plot_grid(p2,p1,nrow=2)
dev.off()



#### ichorCNA ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/")
rm(list=ls())
tissue_order <- c("Brain-Tumour","Breast-Tumour", "Kidney-Tumour","Liver-Tumour","Lung-Tumour","Ovary-Tumour","Pancreas-Tumour","Prostate-Tumour","Colon-Tumour","Stomach-Tumour",
                  "Brain","Breast", "Heart",
                  "Kidney","Liver", "Lung","Ovary", 
                  "Pancreas", "Prostate","Colon","Stomach","Esophagus",
                  "Spleen", "CD4-T-cells", "CD8-T-cells", "Neutrophils", "NK-cells", "B-cells", "Eosinophils", "Monocytes", 
                  "CD34-erythroblasts", "CD34-megakaryocytes","Liver-Cirrhosis","Pancreas-Pancreatitis")
tissue_order <- c(grep("Tumour",tissue_order,value=TRUE,invert = TRUE),grep("Tumour",tissue_order,value=TRUE))


smps <- gsub(".correctedDepth.txt","",list.files(path="cna/caps/", pattern="correctedDepth.txt"))
smps <- smps[order(gsub(".*_","",smps))]
smps <- smps[order(factor(gsub(".*_","",smps),levels = tissue_order))]

all_coverage <- list()
for(smp in smps){
  coverage <- read.delim(paste0("cna/caps/",smp,".correctedDepth.txt"), sep="\t", header=FALSE)
  colnames(coverage) <- c("chr","start","end","log2_TNratio_corrected")
  coverage <- coverage[coverage$chr %in% paste0("chr",c(seq(1,22))),] 
  coverage$chr <- factor(coverage$chr, levels = paste0("chr",c(seq(1,22))))
  coverage$log2_TNratio_corrected <- coverage$log2_TNratio_corrected %>% as.character() %>% as.numeric()
  coverage <- coverage[order(coverage$chr, coverage$start),]
  coverage$idx <- seq(1,nrow(coverage))
  all_coverage[[smp]] <- coverage
}
all_coverage_df <- do.call("rbind", all_coverage)
all_coverage_df$smp <- gsub("\\..*","",rownames(all_coverage_df))
rownames(all_coverage_df) <- NULL
all_coverage_df <- all_coverage_df %>%
  select(chr, start, end, log2_TNratio_corrected, smp)%>%
  pivot_wider(values_from=log2_TNratio_corrected, names_from=smp) 

all_coverage_df <- all_coverage_df[, c(grep("Tumour",colnames(all_coverage_df), invert = TRUE),
                                       grep("Tumour",colnames(all_coverage_df)))] %>% as.data.frame()
all_coverage_df <- all_coverage_df[order(all_coverage_df$chr, as.numeric(all_coverage_df$start)),]
all_coverage_df$chr <- factor(all_coverage_df$chr, levels = paste0("chr",c(seq(1,22))))
all_coverage_df <- all_coverage_df[order(all_coverage_df$chr),]

all_coverage_df <- all_coverage_df[,!grepl("C707005_Heart|C603405_Bone-Marrow|C707007_Thymus|CD563320_Spleen|CD563569_Liver|CD563504_Spleen",colnames(all_coverage_df))]
source("code_tissue_methylome_atlas/rename_sample.r")
colnames(all_coverage_df) <- gsub("-","_",colnames(all_coverage_df))
all_coverage_df <- rename_columns(all_coverage_df) 


ann_row <- data.frame(chr=all_coverage_df$chr)
rownames(ann_row) <- rownames(all_coverage_df)
ann_colors <- list(
  chr = c(chr1 = "#DC5B5A", chr2 = "#989686", chr3 = "#CDE37B", chr4 = "#E29A57", chr5 = "#68DACB", chr6 = "#63E554", chr7 = "#74E091", chr8 = "#C87ADB", chr9 = "#BDE847", chr10 = "#6F55E0",
          chr11 = "#E1C242", chr12 = "#CBE8DE", chr13 = "#DED0D9", chr14 = "#CFABDD", chr15 = "#DB8E9C", chr16 = "#6482AB", chr17 = "#CE44DB", chr18 = "#E4D29F", chr19 = "#BBE7B2", chr20 = "#7BC5E2",
          chr21 = "#DB529F", chr22 = "#7B8CDD")
)
pheatmap(all_coverage_df[,-c(1:3)], cluster_cols = FALSE, fontsize = 4,
         cluster_rows = FALSE, show_rownames = FALSE,
         colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(10),
         breaks = seq(-2,2,4/10),
         annotation_row = ann_row,
         annotation_colors = ann_colors,
         main = "ichorcna corrected depth", filename = "plots/ichorcna_corrected_depth_1Mb.pdf", width = 8, height = 4)

dat <- system("grep 'Tumor Fraction:' cna/caps/*params.txt | sed 's/.params.txt:Tumor Fraction://g;s,cna/caps/,,g'", intern = TRUE) %>%
  read.table(text =., sep = "\t", header = FALSE, stringsAsFactors = FALSE)



dat$V1 <- factor(dat$V1,levels=dat$V1)
dat$tissue <- gsub(".*_","",dat$V1)
dat <- dat[!grepl("C707005_Heart|C603405_Bone-Marrow|C707007_Thymus|CD563320_Spleen|CD563569_Liver|CD563504_Spleen",dat$V1),]

sample_rename <- read.table("code_tissue_methylome_atlas/rename_sample.txt",col.names = c("V1","rename"))
dat <- merge(dat, sample_rename, by="V1")
dat$tissue <- factor(dat$tissue,levels = tissue_order)
dat <- dat[order(dat$tissue),]

dat$rename <- factor(dat$rename,levels=dat$rename)
dat <- dat[order(dat$rename),]



pdf("plots/ichorcna_tumour_fraction.pdf",width = 5,height = 8)
cowplot::plot_grid(
  ggplot(dat[grep("Tumour", dat$V1),],aes(x=tissue,y=V2)) + geom_boxplot(width=0.5,outlier.shape = NA) + geom_point()+
    theme_minimal()+
    theme(axis.text = element_text(size = 8))+
    ylab("tumour fraction estimated by ichorcna") +
    xlab("tissue")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom")+ ylim(0,1) ,
  ggplot(dat[grep("Tumour", dat$V1),],aes(x=rename,y=V2)) + geom_bar(stat="identity") + 
    theme_minimal()+
    theme(axis.text = element_text(size = 8))+
    ylab("tumour fraction estimated by ichorcna") +
    xlab("tissue")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom")+ ylim(0,1) ,
  ggplot(dat[grep("Tumour", dat$V1,invert = TRUE),],aes(x=tissue,y=V2)) + geom_boxplot(width=0.5,outlier.shape = NA) + geom_point()+
    theme_minimal()+
    theme(axis.text = element_text(size = 8))+
    ylab("tumour fraction estimated by ichorcna") +
    xlab("tissue") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom")+ ylim(0,1),
  nrow=3)
dev.off()

