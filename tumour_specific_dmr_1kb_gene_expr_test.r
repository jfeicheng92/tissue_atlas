.libPaths(c(.libPaths(),"/gpfs3/users/ludwig/cfo155/R/x86_64-pc-linux-gnu-library/4.2"))

#### tumour dmb enrich ####
ts_gene_enrich <- function(dmr_name, tissue_specific_gene, all_gene, sel_tissues, out_name, gene_set_order, logp_cap, odds_cap, file_size){
  dmr <- merge(fread(dmr_name),
               fread(paste0(dmr_name,".all_genelist"), 
                     col.names = c("chr","start","end", "selected_tissue", "Hyper_or_Hypo","g_chr","g_start","g_end","gene","dis")),
                     by=c("chr","start","end"))
  print(quantile(dmr$dis))
  enrich_list <- NULL
  max_dis <- 20000
  selN <- 400
  dmr %>%
    filter(dis < max_dis) %>%
    group_by(selected_tissue.x) %>%
    slice_max(order_by = delta_quants, n = selN, with_ties = FALSE)%>%
    select(chr,start,end,gene) %>%
    unique() %>%
    write.table(., gsub(".txt$",paste0(".distancefilter.top",selN,".txt"),dmr_name),sep="\t",quote = FALSE, row.names = FALSE)
  
  
  for(sel_tissue in sel_tissues){
    x1 <- dmr %>%
      filter(selected_tissue.x==sel_tissue) %>%
      filter(dis < max_dis) %>%
      slice_max(order_by = delta_quants, n = selN, with_ties = FALSE)%>%
      select(gene) %>%
      unique()
    
    
    ts_gene <- read.table(tissue_specific_gene, header=FALSE, sep="\t", col.names = c("gene","tissue"))
    total_gene <- read.table(all_gene, header=FALSE, sep="\t", col.names = c("chr","start","end","gene","score","strand"))
    
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
          sel_tissue,
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
  }
  
  result <- data.frame(enrich)
  result <- result %>%
    mutate(tissue=factor(list1,levels=sel_tissues))
  write.table(result, gsub(".txt$",paste0(".distancefilter.top",selN,".enrich.txt"),dmr_name),sep="\t",quote = FALSE, row.names = FALSE)
  dat <- result %>%
    mutate(
      logp   = -log10(pvalue),
      dmr=str_remove(list1, "_.*"),
      tissue = factor(tissue, levels = sel_tissues),
      gene_set = list2,
      gene_set = factor(gene_set, levels = gene_set_order),
    ) %>%
    filter(!is.na(tissue))
  
  dat <- mutate(dat, logp_c = pmin(logp, logp_cap), odds_c = pmin(odds_ratio, odds_cap), sig = logp >= 1.3) # 1.3 ~ p=0.05
  
  p <- ggplot(dat, aes(x = gene_set, y = tissue)) +
    geom_point(aes(size = odds_c, fill = logp_c, alpha = sig),
               shape = 21, stroke = 0.2) +
    scale_size_area(max_size = 7, name = "odds_ratio") +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), name = "-log10(p)") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.1), guide = "none") +
    ggtitle(dmr_name) +
    coord_fixed() +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid = element_blank())
  ggsave(out_name,p, width = file_size, height = file_size)
}


###### list from dmr_call_new1 ######
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/dmr_genes_test/")
dmr_files <- grep("s500.hmC.*hyper|s500.mC.*hypo|s500.umC.*hyper", list.files(pattern = "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.*.bgQ0.05-0.05.tgQ0.25.*groups.all.txt"),value = TRUE)
dmr_files <- grep("s500.hmC.*hypo|s500.mC.*hyper|s500.umC.*hypo", list.files(pattern = "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.*.bgQ0.05-0.05.tgQ0.25.*groups.all.txt"),value = TRUE)

selN <- 2000
for(dmr_file in dmr_files){
  dmr <- fread(dmr_file) %>%
    group_by(selected_tissue) %>%
    slice_max(order_by = delta_quants, n = selN, with_ties = FALSE)
  write.table(dmr[,1:3], gsub(".txt",paste0("top", selN, ".bed"),dmr_file),sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
  dmr <- fread(dmr_file) %>%
    group_by(selected_tissue) %>%
    slice_max(order_by = delta_quants, n = selN, with_ties = FALSE) 
  write.table(dmr%>%select(chr,start,end,tg_quant,delta_quants,selected_tissue,Hyper_or_Hypo), gsub(".txt",paste0("top", selN, ".txt"),dmr_file),sep="\t",quote = FALSE,row.names = FALSE)
}



setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmr_call_new1/dmr_genes_test")
resource_path <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/"


dmr_names <- c(
  "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop2000.txt",
  "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop2000.txt",
  "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.bgQ0.05-0.05.tgQ0.25.hypo_dmrs.bg_quant_modegroups.alltop2000.txt"
  # "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop2000.txt"
  
)
for(dmr_name in dmr_names){
  tissue_order <-c("Brain", "Breast", "Colon", "Kidney", "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Stomach")
  ts_gene_enrich(
    dmr_name = dmr_name,
    tissue_specific_gene = paste0(resource_path,"toil_tumour_specific_gene.diff1.up.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed"),
    sel_tissues = rev(paste0(tissue_order,"-Tumour")),
    gene_set_order = paste0(tissue_order),
    out_name = paste0("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig6/",dmr_name,"tumour.pdf"),
    logp_cap = 20,
    odds_cap = 20,
    file_size = 5
    )
    
}

for(dmr_name in dmr_names){
  tissue_order1 <- c("Brain","Breast","Heart","Kidney","Liver","Lung","Ovary",
                     "Pancreas","Prostate","Colon","Stomach","Esophagus","Spleen",
                     "CD4-T-cells","CD8-T-cells","NK-cells","B-cells","Neutrophils",
                     "Eosinophils","Monocytes","Erythroid-precursors","Megakaryocytes")
  tissue_order2 <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", 
                     "Prostate","Colon","Stomach","Esophagus","Spleen", 
                     "T Helper Cells", "T Cytotoxic Cells", "NK Cells", "B Cells", "Neutrophils", 
                     "Eosinophils", "Monocytes",  "Erythroid Precursor Cells","Megakaryocytes")
  
  ts_gene_enrich(
    dmr_name = dmr_name,
    tissue_specific_gene = paste0(resource_path,"PanglaoDB_immune_GTEx_Tissues_specific_genelist.MANE.txt"),
    all_gene = paste0(resource_path,"MANE.GRCh38.v1.0.refseq_genomic.gene.bed"),
    sel_tissues = rev(tissue_order1),
    gene_set_order=tissue_order2,
    out_name = paste0("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig5/",dmr_name,"normal.pdf"),
    logp_cap = 30,
    odds_cap = 30,
    file_size = 7
  )
  
}


#### DMR distribution around gene ####
library(GenomicRanges)
enrich_res_filenames <- c(
  "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop2000.distancefilter.top400.enrich.txt",
  "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.bgQ0.05-0.05.tgQ0.25.hypo_dmrs.bg_quant_modegroups.alltop2000.distancefilter.top400.enrich.txt",
  "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop2000.distancefilter.top400.enrich.txt"
)
gene_windows <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/MANE.GRCh38.v1.0.refseq_genomic.gene.flank20000.nwin20.bed",col.names = c("chr","start","end","gene","pos","index"))
for(enrich_res_filename in enrich_res_filenames){
  enrich_res <- fread(enrich_res_filename)
  tissue_order1 <- c("Brain","Breast","Heart","Kidney","Liver","Lung","Ovary",
                     "Pancreas","Prostate","Colon","Stomach","Esophagus","Spleen",
                     "CD4-T-cells","CD8-T-cells","NK-cells","B-cells","Neutrophils",
                     "Eosinophils","Monocytes","Erythroid-precursors","Megakaryocytes")
  tissue_order2 <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", 
                     "Prostate","Colon","Stomach","Esophagus","Spleen", 
                     "T Helper Cells", "T Cytotoxic Cells", "NK Cells", "B Cells", "Neutrophils", 
                     "Eosinophils", "Monocytes",  "Erythroid Precursor Cells","Megakaryocytes")
  enrich_res_filter <- data.frame()
  for(i in 1:length(tissue_order1)){
    enrich_res_filter <- rbind(enrich_res_filter,
                               enrich_res %>%
                                 filter(list1==tissue_order1[i] & list2==tissue_order2[i]))
  }
  
  
  gr1 <- enrich_res_filter %>%
    select(list1,overlap_gene) %>%
    separate_rows(overlap_gene, sep = ";") %>%
    mutate(overlap_gene = trimws(overlap_gene)) %>%
    filter(overlap_gene != "") %>%
    dplyr::rename(gene=overlap_gene,
                  tissue=list1)%>%
    merge(.,gene_windows,by=c("gene")) %>%
    select(chr,start,end,tissue,gene,pos,index) %>%
    makeGRangesFromDataFrame(
    ., keep.extra.columns = TRUE,
    seqnames.field = "chr", start.field = "start", end.field = "end",
    starts.in.df.are.0based = FALSE)
  
  gr2 <- fread(gsub(".dis.*",".txt",enrich_res_filename)) %>%
    makeGRangesFromDataFrame(.,
                             keep.extra.columns = TRUE,
                             seqnames.field = "chr", start.field = "start", end.field = "end",
                             starts.in.df.are.0based = FALSE
    )
  # 2) Find overlaps (like bedtools intersect -wa -wb)
  hits <- findOverlaps(gr1, gr2, type = "any", ignore.strand = TRUE)
  
  # 3) Build a joined data.frame with both records (+ exact intersect)
  ov_rng <- pintersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])  # the actual overlap intervals
  
  ov_df <- cbind(
    as.data.frame(gr1[queryHits(hits)])[, c("seqnames","start","end", setdiff(names(mcols(gr1)), character(0)))],
    setNames(
      as.data.frame(gr2[subjectHits(hits)])[, c("seqnames","start","end", setdiff(names(mcols(gr2)), character(0)))],
      c("seqnames_b","start_b","end_b", paste0(names(mcols(gr2)), "_b"))
    )
  )
  ov_df$ov_start <- start(ov_rng)
  ov_df$ov_end   <- end(ov_rng)
  ov_df$ov_width <- width(ov_rng)
  
  # 4) (Optional) Require tissue match
  ov_df_tissue_match <- subset(ov_df, tissue == selected_tissue_b)
  
  
  win <- 3  # window size (odd number recommended)
  
  plot_df <- ov_df_tissue_match %>%
    count(index, name = "Freq") %>%
    mutate(index = as.integer(as.character(index))) %>%
    arrange(index) %>%
    complete(index = seq(min(index), max(index))) %>%
    mutate(Freq = tidyr::replace_na(Freq, 0L)) %>%
    mutate(
      # centered moving average; partial windows allowed at the edges
      ma = slide_dbl(Freq, mean,
                     .before = floor((win-1)/2),
                     .after  = floor((win-1)/2),
                     .complete = FALSE)
    )
  
  p <- ggplot(plot_df, aes(x = index, y = ma, group = 1)) +
    geom_line() +
    ggtitle(enrich_res_filename) +
    coord_cartesian(ylim = c(0, 350)) +  # avoids clipping data
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid = element_blank())
  ggsave(paste0("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig5/",gsub(".txt",".genebody_distribution.pdf",enrich_res_filename)),p,width = 5, height =5 )
}




