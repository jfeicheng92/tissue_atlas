library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(cowplot)
library(reshape2)
library(eulerr)
library(stringr)
options(bitmapType='cairo-png')
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1")



#### Enrichment of DMR linked gene with tissue specific gene ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1")
dmr_files <- grep("s500.hmC.*hyper|s500.mC.*hypo|s500.umC.*hyper", list.files(pattern = "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.*groups.all.txt"),value = TRUE)
for(dmr_file in dmr_files){
  dmr <- fread(dmr_file) %>%
    group_by(selected_tissue) %>%
    slice_max(order_by = delta_quants, n = 500, with_ties = FALSE) 
  write.table(dmr, gsub(".txt","top500.bed",dmr_file),sep="\t",quote = FALSE,row.names = FALSE)
}
for(dmr_file in dmr_files){
  dmr <- fread(dmr_file) %>%
    group_by(selected_tissue) %>%
    slice_max(order_by = delta_quants, n = 500, with_ties = FALSE)
  write.table(dmr[,1:3], gsub(".txt","top500.bed",dmr_file),sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
  dmr <- fread(dmr_file) %>%
    group_by(selected_tissue) %>%
    slice_max(order_by = delta_quants, n = 500, with_ties = FALSE) 
  write.table(dmr%>%select(chr,start,end,selected_tissue,delta_quants,Hyper_or_Hypo), gsub(".txt","top500.txt",dmr_file),sep="\t",quote = FALSE,row.names = FALSE)
  
}





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

setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/dmr_genes/")
genelist_path <- "genelist/"
resource_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/"

genelists <- grep("Tumour|-Pancreatitis|-Cirrhosis",list.files(pattern="*genelist$", path = genelist_path),invert=TRUE,value=TRUE)
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

# for(dmb in c("hmCHyper","umCHyper","mCHypo")){
#   dat <-  result[grep(dmb,result$list1),]
#   dat$logp <- -log(dat$pvalue ,10)
#   dat_w <- dat %>%
#     select(list1, list2, logp) %>%
#     pivot_wider(names_from = list2, values_from = logp) %>% as.data.frame()
#   dat_w$list1 <- gsub("hmCHyper_|umCHyper_|mCHypo_|_genelist","",dat_w$list1) %>% as.factor()
#   tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
#   dat_w <- dat_w[dat_w$list1%in%tissue_order,]
#   dat_w$list1 <- factor(dat_w$list1,levels=tissue_order)
#   dat_w <-dat_w[order(dat_w$list1),]
#   colnames(dat_w) <- gsub(" $","",colnames(dat_w))
#   dat_w <- dat_w %>%
#     select("list1", "Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", "T Helper Cells", "T Cytotoxic Cells", "NK Cells", "B Cells", "Neutrophils", "Eosinophils", "Monocytes",  "Erythroid Precursor Cells","Megakaryocytes")
#   rownames(dat_w) <- dat_w$list1;dat_w$list1<-NULL
#   pheatmap(
#     dat_w %>% t(),
#     color = colorRampPalette(rev(brewer.pal(n = 7, name ="Spectral")))(100),
#     breaks = seq(0,4,4/100), fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
#     display_numbers = round(dat_w,2) %>% t(),
#     filename = paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig5/",dmb,"DMB.PanglaoDB_immune_GTEx_Tissues.logP.pdf"),  width = 6, height = 4.5,
#   )
#   pheatmap(
#     dat_w %>% t(),
#     color = colorRampPalette(rev(brewer.pal(n = 7, name ="Spectral")))(100),
#     scale="column",
#     fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE, 
#     display_numbers = round(dat_w,2) %>% t(),
#     filename = paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig5/",dmb,"DMB.PanglaoDB_immune_GTEx_Tissues.logP.scale.pdf"),  width = 6, height = 4.5
#   )
# }





tissue_order1 <- c("Brain","Breast","Heart","Kidney","Liver","Lung","Ovary",
                  "Pancreas","Prostate","Colon","Stomach","Esophagus","Spleen",
                  "CD4-T-cells","CD8-T-cells","NK-cells","B-cells","Neutrophils",
                  "Eosinophils","Monocytes","Erythroid-precursors","Megakaryocytes")
tissue_order2 <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", 
                   "Prostate","Colon","Stomach","Esophagus","Spleen", 
                   "T Helper Cells", "T Cytotoxic Cells", "NK Cells", "B Cells", "Neutrophils", 
                   "Eosinophils", "Monocytes",  "Erythroid Precursor Cells","Megakaryocytes")
pdf(paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig5/all.DMB.PanglaoDB_immune_GTEx_Tissues.logP.scale.pdf"), width = 8, height = 18)
dat <- result %>%
  # filter(str_detect(list1, fixed(dmb))) %>%
  mutate(
    logp   = -log10(pvalue),
    dmr=str_remove(list1, "_.*"),
    dmr=factor(dmr, levels=c("umCHyper","hmCHyper","mCHypo")),
    tissue = str_remove(str_remove(list1, "^(hmCHyper|umCHyper|mCHypo)_"), "_genelist"),
    tissue = factor(tissue, levels = rev(tissue_order1)),
    gene_set = list2,
    gene_set = factor(gene_set, levels = tissue_order2),
  ) %>%
  filter(!is.na(tissue))

logp_cap <- 30  # cap extreme values
odds_cap <- 25
dat <- mutate(dat, logp_c = pmin(logp, logp_cap), odds_c = pmin(odds_ratio, odds_cap), sig = logp >= 1.3) # 1.3 ~ p=0.05

ggplot(dat, aes(x = gene_set, y = tissue)) +
  facet_wrap(dmr~., nrow=3)+
  geom_point(aes(size = odds_c, fill = logp_c, alpha = sig),
             shape = 21, stroke = 0.2) +
  scale_size_area(max_size = 7, name = "odds_ratio") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(11, "Spectral"), name = "-log10(p)") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
  coord_fixed() +
  # labs(x = NULL, y = NULL, title = paste(dmb, "enrichment")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())
dev.off()















#### Plot raw expression ####
##### Solid tissue #####
# Load genomic annotations linking windows/DMRs to genes.
RNA_anno <- read.table("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/hg38_ws1000.s500.protein_coding_gene.5k.bed")
colnames(RNA_anno)  <- c("chr", "start","end","gene_chr","gene_start","gene_end","strand","gene","info","dis","pos" )
# Load gene expression matrix (TPM medians), then keep only ID, gene name, and GTEx columns.
RNA_exp <- data.table::fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/gene_expr/all.rsem_gene_tpm_median.raw.genename.csv") |>
  as_tibble() |>
  select(id, gene, contains("gtex"))

mat <- as.matrix(RNA_exp |> select(contains("gtex")))
row_z <- t(scale(t(mat)))   # z per row
RNA_exp_nor <- bind_cols(RNA_exp |> select(id, gene), as.data.frame(row_z))

setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1")
bg_tg_q <- "bgQ0.05-0.05.tgQ0.25"
dmr_files <- grep("s500.hmC.*hyper|s500.mC.*hypo|s500.umC.*hyper", list.files(pattern = paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.*",bg_tg_q,".*groups.all.txt")),value = TRUE)
sel_N <- 400
z_list <- list()
for( dmr_file in dmr_files){
  print(dmr_file)
  dmr <- fread(dmr_file) %>%
    group_by(selected_tissue) %>%
    slice_max(order_by = delta_quants, n = sel_N, with_ties = FALSE) %>%
    select(chr, start, end, selected_tissue) %>%
    as.data.frame()
  # link dmr to gene
  dmr_gene <- merge(dmr,RNA_anno, by= c("chr", "start","end")) %>%
    select(selected_tissue, gene) %>%
    unique()
  # get gene expression
  dmr_gene <- dmr_gene[order(dmr_gene$selected_tissue),]
  dmr_gene_expr <- merge(dmr_gene, RNA_exp_nor, by=c("gene"))
  
  select_tissue_order <- colnames(dmr_gene_expr) %>% gsub("_gtex.*","",.) %>% .[-c(1:3)]
  
  dmr_gene_expr <- dmr_gene_expr %>%
    filter(selected_tissue %in% select_tissue_order) %>%
    mutate(selected_tissue = factor(selected_tissue, levels = select_tissue_order)) %>%
    arrange(selected_tissue) %>%
    select(selected_tissue,contains("median"))
  df <- dmr_gene_expr
  numeric_cols <- grep("median", names(df), value = TRUE)
  
  df_zscore <- df %>%
    group_by(selected_tissue) %>%
    summarise(across(all_of(numeric_cols), mean, na.rm = TRUE)) %>% as.data.frame()
  rownames(df_zscore) <- df_zscore$selected_tissue
  pheatmap(df_zscore[,-1], cluster_rows = FALSE, cluster_cols = FALSE,scale="row", 
           main=gsub(".*s500.|bg_quant_mode.*","",dmr_file),
           filename = paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig5/",gsub(".txt","linked_gene.solid_tissue.pdf", dmr_file)))
  z_list[[dmr_file]] <- df_zscore
}

df_zscore_all <- do.call(rbind,z_list) %>%
  as.data.frame() %>%
  rownames_to_column("dmr") %>%
  mutate(
    dmr=str_extract(dmr, "(?<=s500\\.).*?(?=\\.bg)")
  ) %>%
  pivot_longer(
    cols = ends_with("median"),
    values_to = "value"
  ) %>%
  mutate(
    selected_tissue = factor(selected_tissue, level=rev(c("Brain","Breast","Colon","Kidney","Liver","Lung","Ovary","Pancreas","Prostate","Stomach"))),
    dmr=factor(dmr, level=c("umC","hmC","mC"))
    )
  
pdf(paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig5/",bg_tg_q,".topN.",sel_N,".solid_tissue.pdf"),width=15, height=5)
ggplot(df_zscore_all, aes(x = name, y = selected_tissue , fill=value)) +
  facet_wrap(dmr~., nrow=1)+
  geom_tile() +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdBu")), name = "z-score expression") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())
dev.off()





















##### Blood cell #####
# Load genomic annotations linking windows/DMRs to genes.
RNA_anno <- read.table("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/hg38_ws1000.s500.protein_coding_gene.5k.bed")
colnames(RNA_anno)  <- c("chr", "start","end","gene_chr","gene_start","gene_end","strand","gene","info","dis","pos" )
# Load gene expression matrix (TPM medians), then keep only ID, gene name, and GTEx columns.
RNA_exp <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/immune_cells/rna_immune_cell.tsv")
select_cell_types <- c("memory CD4 T-cell","memory CD8 T-cell",
                       "NK-cell","memory B-cell", "neutrophil","classical monocyte","eosinophil"
                       )
sel_N <- 400
RNA_exp <- RNA_exp %>%
  select(Gene, `Gene name`, `Immune cell`, nTPM) %>%
  filter(`Immune cell` %in% select_cell_types) %>%   # fix typo here
  pivot_wider(
    names_from = `Immune cell`,
    values_from = nTPM,
    values_fill = 0
  ) %>%
  relocate(any_of(select_cell_types), .after = `Gene name`)

mat <- as.matrix(RNA_exp[,-c(1,2)])
row_z <- t(scale(t(mat)))   # z per row
RNA_exp_nor <- bind_cols(RNA_exp |> select(Gene, `Gene name`), as.data.frame(row_z))


dup_names <- RNA_exp_nor %>%
  count(`Gene name`) %>%
  filter(n > 1) %>%
  arrange(desc(n))
RNA_exp_nor <- RNA_exp_nor[!RNA_exp_nor$`Gene name` %in% dup_names$`Gene name`, ]
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1")
bg_tg_q <- "bgQ0.05-0.05.tgQ0.25"
dmr_files <- grep("s500.hmC.*hyper|s500.mC.*hypo|s500.umC.*hyper", list.files(pattern = paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.*",bg_tg_q,".*groups.all.txt")),value = TRUE)
select_tissue_order <- c("CD4-T-cells","CD8-T-cells",
                           "NK-cells","B-cells","Neutrophils","Monocytes","Eosinophils")
z_list <- list() 
for( dmr_file in dmr_files){
  print(dmr_file)
  # choose top 400 dmr
  dmr <- fread(dmr_file) %>%
    group_by(selected_tissue) %>%
    slice_max(order_by = delta_quants, n = sel_N, with_ties = FALSE) %>%
    select(chr, start, end, selected_tissue) %>%
    as.data.frame()
  # link dmr to gene
  dmr_gene <- merge(dmr,RNA_anno, by= c("chr", "start","end")) %>%
    select(selected_tissue, gene) %>%
    unique()
  # get gene expression
  dmr_gene <- dmr_gene[order(dmr_gene$selected_tissue),]
  dmr_gene_expr <- merge(dmr_gene, RNA_exp_nor, by.x=c("gene"), by.y=c("Gene name"),all.x = TRUE)
  dmr_gene_expr <- dmr_gene_expr %>%
    filter(selected_tissue %in% select_tissue_order) %>%
    mutate(selected_tissue = factor(selected_tissue, levels = select_tissue_order)) %>%
    arrange(selected_tissue) 
  df <- dmr_gene_expr
  numeric_cols <- names(df)[-c(1:3)]
  
  df_zscore <- df %>%
    group_by(selected_tissue) %>%
    summarise(across(all_of(numeric_cols), mean, na.rm = TRUE)) %>% as.data.frame()
  rownames(df_zscore) <- df_zscore$selected_tissue
  pheatmap(df_zscore[,-1], cluster_rows = FALSE, cluster_cols = FALSE,scale="row", 
           main=gsub(".*s500.|bg_quant_mode.*","",dmr_file),
           filename = paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig5/",gsub(".txt","linked_gene.immune.pdf", dmr_file)))
  z_list[[dmr_file]] <- df_zscore
}
df_zscore_all <- do.call(rbind,z_list) %>%
  as.data.frame() %>%
  rownames_to_column("dmr") %>%
  mutate(
    dmr=str_extract(dmr, "(?<=s500\\.).*?(?=\\.bg)")
  ) %>%
  pivot_longer(
    cols = select_cell_types,
    values_to = "value"
  ) %>%
  mutate(
    selected_tissue = factor(selected_tissue, level=rev(select_tissue_order)),
    dmr=factor(dmr, level=c("umC","hmC","mC")),
    name=factor(name,level=select_cell_types)
  )

p <- ggplot(df_zscore_all, aes(x = name, y = selected_tissue , fill=value)) +
  facet_wrap(dmr~., nrow=1)+
  geom_tile() +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdBu")), name = "z-score expression") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())
p
pdf(paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig5/",bg_tg_q,".topN.",sel_N,".immune_cells.pdf"),width=10, height=4)
print(p)
dev.off()


