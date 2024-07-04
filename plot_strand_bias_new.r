setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/")
###### by hmC ######
samples <- c("C505016_Brain","C603405_Bone-Marrow","C707004_Brain","C707005_Heart","C707006_Heart","C707007_Thymus","C707008_Brain-Tumour","C707009_Brain-Tumour","CD563106_Liver-Cirrhosis","CD563137_Breast","CD563162_Stomach","CD563176_Liver-Tumour","CD563257_Liver","CD563267_Prostate","CD563278_Liver-Cirrhosis","CD563304_Breast","CD563320_Spleen","CD563419_Colon","CD563430_Stomach","CD563503_Breast-Tumour","CD563504_Spleen","CD563544_Ovary","CD563553_Pancreas-Pancreatitis","CD563569_Liver","CD563610_Prostate","CD563649_Colon-Tumour","CD563663_Colon","CD563678_Esophagus","CD563685_Prostate","CD563761_Prostate-Tumour","CD563776_Kidney","CD563778_Kidney","CD563808_Ovary-Tumour","CD563854_Lung","CD563857_Kidney-Tumour","CD563880_Spleen","CD563882_Prostate-Tumour","CD563930_Lung","CD563984_Kidney","CD564011_Pancreas","CD564068_Breast","CD564076_Breast-Tumour","CD564082_Liver","CD564127_Pancreas-Tumour","CD564145_Colon-Tumour","CD564146_Liver-Tumour","CD564159_Colon","CD564191_Ovary","CD564208_Liver-Tumour","CD564241_Kidney","CD564242_Prostate","CD564285_Lung-Tumour","CD564295_Ovary","CD564368_Breast","CD564404_Pancreas","CD564414_Prostate-Tumour","CD564422_Pancreas-Tumour","CD564511_Lung","CD564572_Lung-Tumour","CD564575_Breast-Tumour","CD564588_Stomach-Tumour","CD564596_Stomach","CD564659_Stomach-Tumour","CD564797_Liver","CD564801_Liver-Cirrhosis","CD564844_Pancreas","CD564847_Liver-Cirrhosis","CD564876_Kidney-Tumour","CD564888_Ovary-Tumour","CD564901_Lung","CD564934_Liver-Tumour","CD564953_Kidney-Tumour","CD564971_Colon-Tumour","CD564980_Prostate-Tumour","CD564986_Esophagus","CD565017_Spleen","CD565018_Ovary-Tumour","CD565042_Stomach","CD565056_Breast-Tumour","CD565088_Ovary-Tumour","CD565099_Lung-Tumour","CD565136_Esophagus","CD565153_Kidney-Tumour","CD565189_Colon","CD565203_Colon-Tumour","CD565252_Esophagus","CD565285_Lung-Tumour","CD565341_Pancreas","D30_CD34-erythroblasts","D30_CD34-megakaryocytes","D37_CD34-erythroblasts","D66_CD34-erythroblasts","D66_CD34-megakaryocytes","D70_CD34-erythroblasts","UKVAC-001-5_B-cells","UKVAC-001-5_CD4-T-cells","UKVAC-001-5_CD8-T-cells","UKVAC-001-5_Eosinophils","UKVAC-001-5_Monocytes","UKVAC-001-5_NK-cells","UKVAC-001-5_Neutrophils","UKVAC-003-6_B-cells","UKVAC-003-6_CD4-T-cells","UKVAC-003-6_CD8-T-cells","UKVAC-003-6_Eosinophils","UKVAC-003-6_Monocytes","UKVAC-003-6_NK-cells","UKVAC-003-6_Neutrophils","UKVAC-049-3_B-cells","UKVAC-049-3_CD4-T-cells","UKVAC-049-3_CD8-T-cells","UKVAC-049-3_Eosinophils","UKVAC-049-3_Monocytes","UKVAC-049-3_NK-cells","UKVAC-049-3_Neutrophils","UKVAC-140_B-cells","UKVAC-140_CD4-T-cells","UKVAC-140_CD8-T-cells","UKVAC-140_Eosinophils","UKVAC-140_Monocytes","UKVAC-140_NK-cells","UKVAC-140_Neutrophils")
samples <- samples[!samples%in%c('CD563569_Liver','CD563504_Spleen','CD563320_Spleen','C707005_Heart','C707007_Thymus','C603405_Bone-Marrow')]
samples <- samples[!grepl("Tumour|Cirrhosis|Pancreatitis|Esophagus|Spleen|Heart",samples)]
sel_samples <- samples[grep("^CD|^C|^UK",samples)]
dat1 <- fread("meth/Genebody_CpG_meth_nonstrand.txt") %>% melt(id.vars=c("gene","gene_strand"))
dat2 <- fread("meth/Genebody_CpG_meth_strand.txt") %>% 
  pivot_longer(cols = contains(c("_CAPS","_TAPS")), names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = c(cg_strand),values_from = "value") %>%
  as.data.frame()



dat <- merge(dat1,dat2,by=c("gene","gene_strand","variable"))
colnames(dat) <- c("gene","strand","smp","non_strand_specific","pos","neg")
dat <- dat[dat$smp%in%c(paste0(sel_samples,"_CAPS"),paste0(sel_samples,"_TAPSbeta")),]

dat$mark <- "nonBrain"
dat$mark[grep("Brain",dat$smp)] <- "Brain"
seldat <- dat[grep("CAPS",dat$smp),]
seldat <- seldat[complete.cases(seldat),]
seldat$diff <- seldat$pos-seldat$neg
seldat$level <- cut(seldat$non_strand_specific, breaks=c(0,10,20,Inf),include.lowest = TRUE,labels=c("0≤hmC≤10","10<hmC≤20","hmC>20"))
seldat$strand <- factor(seldat$strand,levels=c("+","-"))
p1 <- ggplot(seldat[seldat$mark=="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("Brain\nhmC") + xlab("hmC(pos-neg)") +
  scale_fill_aaas()
p2 <- ggplot(seldat[seldat$mark!="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("nonBrain\nhmC") + xlab("hmC(pos-neg)")+
  scale_fill_aaas()
level_table <- table(seldat%>%select(mark,level)) %>% t() %>% as.data.frame() %>%
  pivot_wider(names_from = c(mark),values_from = "Freq")
level_table$Brain_p <- level_table$Brain/sum(level_table$Brain)
level_table$nonBrain_p <- level_table$nonBrain/sum(level_table$nonBrain)
p3 <- level_table %>% select(level, nonBrain, Brain) %>%
  melt(by=("level")) %>%
  ggplot(aes(x = level, y = variable)) +
  geom_tile(fill="white") +
  geom_text(aes(label = sprintf("%d", value)), color = "black", size = 3) +
  guides(col = guide_legend(ncol = 3)) +
  theme(legend.position = "right") +
  xlab("hmC level") + ylab("") +
  theme_classic()+
  theme(axis.text=element_text(size=10))
p123 <- cowplot::plot_grid(p1,p2,p3,nrow=3,rel_heights = c(1,1,0.4))

seldat <- dat[grep("TAPSbeta",dat$smp),]
seldat <- seldat[complete.cases(seldat),]
seldat$diff <- seldat$pos-seldat$neg
seldat$level <- cut(seldat$non_strand_specific, breaks=c(0,30,60,Inf),include.lowest = TRUE,labels=c("0≤mC≤30","30<mC≤60","mC>60"))
seldat$strand <- factor(seldat$strand,levels=c("+","-"))
p4 <- ggplot(seldat[seldat$mark=="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("Brain\nmC")+ xlab("mC(pos-neg)") +
  scale_fill_aaas()
p5 <- ggplot(seldat[seldat$mark!="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("nonBrain\nmC")+ xlab("mC(pos-neg)")+
  scale_fill_aaas()
p45 <- cowplot::plot_grid(p3,p4,nrow=2)
level_table <- table(seldat%>%select(mark,level)) %>% t() %>% as.data.frame() %>%
  pivot_wider(names_from = c(mark),values_from = "Freq")
level_table$Brain_p <- level_table$Brain/sum(level_table$Brain)
level_table$nonBrain_p <- level_table$nonBrain/sum(level_table$nonBrain)

p6 <- level_table %>% select(level, nonBrain, Brain) %>%
  melt(by=("level")) %>%
  ggplot(aes(x = level, y = variable)) +
  geom_tile(fill="white") +
  geom_text(aes(label = sprintf("%d", value)), color = "black", size = 3) +
  guides(col = guide_legend(ncol = 3)) +
  theme(legend.position = "right") +
  xlab("mC level") + ylab("") +
  theme_classic()+
  theme(axis.text=element_text(size=10))
p456 <- cowplot::plot_grid(p4,p5,p6,nrow=3,rel_heights = c(1,1,0.4))

pdf("strand_bias/strand_bias_on_allgene_sort_byhmC.pdf", width = 20, height = 5)
cowplot::plot_grid(p123,p456,nrow=1)
dev.off()

###### by tpm ######
meth <- fread("meth/Genebody_CpG_meth_strand.txt") %>% 
  pivot_longer(cols = contains(c("_CAPS","_TAPS")), names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = c(cg_strand),values_from = "value") %>%
  as.data.frame()

expr <- fread("gene_expr/select_tissue_tpm.txt") %>% 
  pivot_longer(cols = contains("tpm"), names_to = "variable", values_to = "value")
meth$tissue <- lapply(meth$variable,function(x){strsplit(x,"_") %>% unlist() %>% .[2]}) %>% unlist()
expr$tissue <- lapply(expr$variable,function(x){strsplit(x,"_") %>% unlist() %>% .[1]}) %>% unlist()

dat <- merge(meth,expr%>%select(gene,value,tissue),by=c("gene","tissue"))
colnames(dat) <- c("gene","tissue","strand","smp","pos","neg","tpm")
dat <- dat[dat$smp%in%c(paste0(sel_samples,"_CAPS"),paste0(sel_samples,"_TAPSbeta")),]

dat$mark <- "nonBrain"
dat$mark[grep("Brain",dat$smp)] <- "Brain"

seldat <- dat[grep("CAPS",dat$smp),]
seldat <- seldat[complete.cases(seldat),]
seldat$diff <- seldat$pos-seldat$neg
seldat$level <- cut(seldat$tpm, breaks=c(0,1,10,Inf),include.lowest = TRUE)
seldat$strand <- factor(seldat$strand,levels=c("+","-"))
p1 <- ggplot(seldat[seldat$mark=="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("Brain\nhmC") + xlab("hmC(pos-neg)") +
  scale_fill_aaas()
p2 <- ggplot(seldat[seldat$mark!="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("nonBrain\nhmC") + xlab("hmC(pos-neg)") +
  scale_fill_aaas() 
level_table <- table(seldat%>%select(mark,level)) %>% t() %>% as.data.frame() %>%
  pivot_wider(names_from = c(mark),values_from = "Freq")
level_table$Brain_p <- level_table$Brain/sum(level_table$Brain)
level_table$nonBrain_p <- level_table$nonBrain/sum(level_table$nonBrain)
p3 <- level_table %>% select(level, nonBrain, Brain) %>%
  melt(by=("level")) %>%
  ggplot(aes(x = level, y = variable)) +
  geom_tile(fill="white") +
  geom_text(aes(label = sprintf("%d", value)), color = "black", size = 3) +
  guides(col = guide_legend(ncol = 3)) +
  theme(legend.position = "right") +
  xlab("TPM") + ylab("") +
  theme_classic()+
  theme(axis.text=element_text(size=10))
p123 <- cowplot::plot_grid(p1,p2,p3,nrow=3,rel_heights = c(1,1,0.4))

seldat <- dat[grep("TAPSbeta",dat$smp),]
seldat <- seldat[complete.cases(seldat),]
seldat$diff <- seldat$pos-seldat$neg
seldat$level <- cut(seldat$tpm, breaks=c(0,1,10,Inf),include.lowest = TRUE)
seldat$strand <- factor(seldat$strand,levels=c("+","-"))
p4 <- ggplot(seldat[seldat$mark=="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("Brain\nmC")+ xlab("mC(pos-neg)") +
  scale_fill_aaas()
p5 <- ggplot(seldat[seldat$mark!="Brain",], aes(x=diff, fill=strand)) + 
  geom_density(alpha=0.5) +
  xlim(-5, 5) +
  facet_grid(~level)+
  theme(legend.position = "bottom") +
  theme_classic()+
  ylab("nonBrain\nmC")+ xlab("mC(pos-neg)") +
  scale_fill_aaas()
p45 <- cowplot::plot_grid(p3,p4,nrow=2)
level_table <- table(seldat%>%select(mark,level)) %>% t() %>% as.data.frame() %>%
  pivot_wider(names_from = c(mark),values_from = "Freq")
level_table$Brain_p <- level_table$Brain/sum(level_table$Brain)
level_table$nonBrain_p <- level_table$nonBrain/sum(level_table$nonBrain)

p6 <- level_table %>% select(level, nonBrain, Brain) %>%
  melt(by=("level")) %>%
  ggplot(aes(x = level, y = variable)) +
  geom_tile(fill="white") +
  geom_text(aes(label = sprintf("%d", value)), color = "black", size = 3) +
  guides(col = guide_legend(ncol = 3)) +
  theme(legend.position = "right") +
  xlab("TPM") + ylab("") +
  theme_classic()+
  theme(axis.text=element_text(size=10))
p456 <- cowplot::plot_grid(p4,p5,p6,nrow=3,rel_heights = c(1,1,0.4))

pdf("strand_bias/strand_bias_on_allgene_sort_bytpm.pdf", width = 25, height = 5)
cowplot::plot_grid(p123,p456,nrow=1)
dev.off()
