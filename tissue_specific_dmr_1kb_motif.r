library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(cowplot)
library(reshape2)
library(stringr)
library(purrr)
library(tibble)
options(bitmapType='cairo-png')
##### Enrich motif (Normal) #####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/motif/")
tissue_order1 <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                  "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
# tissue_order <- c(tissue_order1) #, tissue_order2)
# tissue_order <- tissue_order[order(tissue_order)]
tissue_order <- rev(tissue_order1)
motif_cluster <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/homer/motifs/all_motif_tomtom.all.cluster.name.txt",col.names=c("id","motif","cluster"))
allmotif <- fread("all_motif.txt")

allmotif <- merge(allmotif,motif_cluster,by.x=c("Motif_Name"),by.y = "motif")
allmotif <- allmotif %>%
  mutate(
    bed=gsub(".bgQ0.05-0.05.tgQ0.25.|_dmrs.bg_quant_modegroups.|alltop500.txt|.bed|_Hyper|_Hypo","",bed)
    
  ) %>%
  separate(bed, into = c("prefix","tissue"),
           sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
  filter(tissue %in% tissue_order) %>%
  mutate(tissue=factor(tissue,levels=rev(tissue_order))) %>%
  dplyr::rename(., c("logP"="Log_P-value", "q_value"="q-value(Benjamini)"))


allmotif$odds <- gsub("%","",allmotif$Pct_of_targetSequenceswithMotif) %>% as.numeric()/
  gsub("%","",allmotif$Pct_of_BackgroundSequenceswithMotif) %>% as.numeric()

motif_select <- c(
  "K562-GATA1-ChIP-Seq(GSE18829)",
  "CEBP(bZIP)/ThioMac-CEBPb-ChIP-Seq(GSE21512)/Homer",
  "PU.1-IRF(ETS:IRF)/Bcell-PU.1-ChIP-Seq(GSE21512)/Homer",
  "CD8-Tbet-ChIP-Seq(GSE33802)",
  "Jurkat-RUNX1-ChIP-Seq(GSE29180)",
  "GM12878-TCF7-ChIP-Seq(Encode)",
  "VCaP-ERG-ChIP-Seq(GSE14097)",
  "Keratinocyte-p63-ChIP-Seq(GSE17611)",
  "Heart-Gata4-ChIP-Seq(GSE35151)",
  "mES-Cdx2-ChIP-Seq(GSE14586)",
  "ProstateTumor-HOXB13-ChIP-Seq(GSE56288)",
  "PDAC-ZEB1-ChIP-Seq(GSE64557)",
  "H295R-Nr5a1-ChIP-Seq(GSE44220)",
  "LungAC-Nkx2.1-ChIP-Seq(GSE43252)",
  "DR1/HepG2-HNF4a-ChIP-Seq(GSE25021)",
  "PDAC-HNF1B-ChIP-Seq(GSE64557)",
  "HEK293-Mef2b.V5-ChIP-Seq(GSE67450)",
  "MCF7-TFAP2C-ChIP-Seq(GSE21234)",
  "ESC-SOX21-ChIP-Seq(GSE110505)")

enriched_motif <- unlist(
  lapply(motif_select, function(x) {
    grep(x, unique(allmotif$Motif_Name), value = TRUE, fixed = TRUE)
  })
) %>%
  as.data.frame() %>%
  dplyr::rename(Motif_Name = 1)

logp_cap <- 30
odds_cap <- 4
dat <- merge(enriched_motif%>%select("Motif_Name"), 
             allmotif, 
             by=c("Motif_Name")) %>%
  mutate(logP=-logP)%>%
  mutate(logp_c = pmin(logP, logp_cap), odds_c = pmin(odds, odds_cap), sig = logP >= 3) 

p <- dat %>%
  mutate(Motif_Name=factor(Motif_Name, levels=rev(enriched_motif$Motif_Name)),
         prefix=factor(prefix, c("mChypo","umChyper","hmChyper")),
         tissue=factor(tissue, levels=tissue_order)) %>%
  ggplot(aes(x = Motif_Name, y = tissue)) +
  facet_wrap(prefix~., nrow=1)+
  geom_point(aes(size = odds_c, fill = logp_c, alpha = sig),
             shape = 21, stroke = 0.2) +
  scale_size_area(max_size = 7, name = "odds_ratio") +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), name = "-log10(p)") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.1), guide = "none") +
  coord_fixed() +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        panel.grid = element_blank())
p
pdf("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig4/motif_enrich.pdf",width = 12, height = 12)
print(p)
dev.off()
odds_merged <- merge(dat %>%filter(prefix=="mChypo") %>% select(Motif_Name, prefix,odds,tissue),
      dat %>%filter(prefix=="hmChyper") %>% select(Motif_Name, prefix,odds,tissue), by=c("Motif_Name","tissue"))
  
ggplot(odds_merged, aes(x=odds.x, y=odds.y)) +geom_point()
cor(odds_merged$odds.x, odds_merged$odds.y, method ="pearson")
# 0.4181961

##### Enrich motif (Tumour) #####
tissue_order2 <-paste0(c("Brain", "Breast", "Colon", "Kidney",  "Lung", "Ovary", "Pancreas", "Prostate"),"-Tumour") # "Liver",, "Stomach"
# tissue_order2 <-paste0(c("Brain", "Breast", "Colon", "Kidney",  "Lung", "Ovary", "Pancreas", "Prostate")) # "Liver",, "Stomach"

allmotif <- fread("all_motif.txt")

allmotif <- merge(allmotif,motif_cluster,by.x=c("Motif_Name"),by.y = "motif")
allmotif <- allmotif %>%
  mutate(
    bed=gsub(".bgQ0.05-0.05.tgQ0.25.|_dmrs.bg_quant_modegroups.|alltop500.txt|.bed|_Hyper|_Hypo","",bed)
    
  ) %>%
  separate(bed, into = c("prefix","tissue"),
           sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
  filter(tissue %in% tissue_order2) %>%
  mutate(tissue=factor(tissue,levels=rev(tissue_order2))) %>%
  dplyr::rename(., c("logP"="Log_P-value", "q_value"="q-value(Benjamini)"))
allmotif$odds <- gsub("%","",allmotif$Pct_of_targetSequenceswithMotif) %>% as.numeric()/
  gsub("%","",allmotif$Pct_of_BackgroundSequenceswithMotif) %>% as.numeric()

motif_select <- c(
  # "HFSC-Lhx2-ChIP-Seq(GSE48068)",
  "GBM-ATF3-ChIP-Seq(GSE33912)",
  "GATA3(Zf),DR8/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer",
  "GM12878-TCF7-ChIP-Seq(Encode)",
  "PDAC-HNF1B-ChIP-Seq(GSE64557)",
  "MCF7-HIF1a-ChIP-Seq(GSE28352)",
  "LungAC-Nkx2.1-ChIP-Seq(GSE43252)",
  "Kidney-WT1-ChIP-Seq(GSE90016)",
  "HepG2-TEAD3-ChIP-Seq(Encode)",
  "Striatum-Fra2-ChIP-Seq(GSE43429)",
  "ProstateTumor-HOXB13-ChIP-Seq(GSE56288)",
  "FOXA1:AR(Forkhead,NR)/LNCAP-AR-ChIP-Seq(GSE27824)/Homer",
  "MCF7-FOXM1-ChIP-Seq(GSE72977)"
) %>% rev()

enriched_motif <- unlist(
  lapply(motif_select, function(x) {
    grep(x, unique(allmotif$Motif_Name), value = TRUE, fixed = TRUE)
  })
) %>%
  as.data.frame() %>%
  dplyr::rename(Motif_Name = 1)

logp_cap <- 30
odds_cap <- 4
dat <- merge(enriched_motif%>%select("Motif_Name"), 
             allmotif, 
             by=c("Motif_Name")) %>%
  mutate(logP=-logP)%>%
  mutate(logp_c = pmin(logP, logp_cap), odds_c = pmin(odds, odds_cap), sig = logP >= 3) 

p <- dat %>%
  mutate(Motif_Name=factor(Motif_Name, levels=rev(enriched_motif$Motif_Name)),
         prefix=factor(prefix, c("mChypo","umChyper","hmChyper")),
         tissue=factor(tissue, levels=rev(tissue_order2))) %>%
  ggplot(aes(x = Motif_Name, y = tissue)) +
  facet_wrap(prefix~., nrow=1)+
  geom_point(aes(size = odds_c, fill = logp_c, alpha = sig),
             shape = 21, stroke = 0.2) +
  scale_size_area(max_size = 7, name = "odds_ratio") +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), name = "-log10(p)") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.1), guide = "none") +
  coord_fixed() +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        panel.grid = element_blank())

pdf("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig6/motif_enrich.pdf",width = 9, height = 8)
print(p)
dev.off()

