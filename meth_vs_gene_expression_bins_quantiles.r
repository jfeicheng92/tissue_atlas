library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(cowplot)
library(reshape2)  # for melt()
options(bitmapType='cairo-png')
setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/")

#### README ####
# This script generates plots of average DNA methylation levels
# stratified by gene expression quantiles.
#################

#### Functions ####
#----Methylation----#
deeptools_load <- function(fname){
  mat <- fread(fname,skip = 1)
  tmp <-  readLines(fname, n = 1) %>% strsplit('],"') %>% unlist()
  params <- data.frame(par=lapply(gsub('\\{|\\}|\\"|\\[|\\]',"", tmp), function(x)unlist(strsplit(x,":"))[[1]]) %>% unlist,
                       value=lapply(gsub('\\{|\\}|\\"|\\[|\\]',"", tmp), function(x)unlist(strsplit(x,":"))[[2]]) %>% unlist)
  # Extract samples from comma-separated string
  samples <- unlist(strsplit(params[10, 2], ","))
  
  # Extract nbins from comma-separated value (2nd element)
  nbins <- strsplit(params[11, 2], ",") %>% unlist() %>% .[2] %>% as.numeric()
  
  # Generate repeated names with bin index
  sample_bins <- unlist(lapply(samples, function(s) {
    paste0(s, "_bin", seq_len(nbins))
  }))
  
  colnames(mat) <- c("chr", "start", "end", "gene", "info", "strand", sample_bins)
  return(mat)
}


fname_lists <- c("meth_around_gene.umC.mlml.mincov10_common.mat.gz", "meth_around_gene.hmC.mlml.mincov10_common.mat.gz", "meth_around_gene.mC.mlml.mincov10_common.mat.gz")
for(fname in fname_lists){
  meth_marker <- gsub(".*gene.|.mlml.*","",fname)
  meth <- deeptools_load(paste0("meth/combine_mC_hmC/",fname)) %>% as.data.frame()
  
  ##### calculate mean methylation per tissue type #####
  tissues_bins <- sub("^.*?_([^_].*?)$", "\\1", colnames(meth)[-c(1:6)]) %>% unique()
  
  for (pattern in tissues_bins) {
    col_name <- paste0(pattern, "_avg")
    matched_cols <- grep(paste0("_", pattern, "$"), colnames(meth), value = TRUE)
    
    #cat("Matched columns for", pattern, ":", matched_cols, "\n")
    
    if (length(matched_cols) > 1) {
      meth[[col_name]] <- rowMeans(meth[, matched_cols], na.rm = TRUE)
    } else if (length(matched_cols) == 1) {
      meth[[col_name]] <- meth[[matched_cols]]
    } else {
      meth[[col_name]] <- NA
    }
    
  }

  #### compare methylation and expression ####
  ## solid tissue ##
  select_tissues <- grep("Tumour.*_avg", colnames(meth), value = TRUE) %>%
    gsub("-.*", "", .) %>% unique()
  
  gene_exp <- read.csv("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/gene_expr/all.rsem_gene_tpm_median.raw.genename.csv") # 

  
  summary_df <- data.frame()
  for (select_tissue in select_tissues) {
    ### ---- Non-tumour ---- ###
    select_gene_exp <- gene_exp[, c(1, 2, grep(paste0("^", select_tissue, "_gtex"), names(gene_exp)))]
    select_gene_meth <- meth[, c(4, grep(paste0("^", select_tissue, "_.*avg$"), names(meth)))]
    
    select_tissue_mat_n <- merge(select_gene_exp, select_gene_meth, by = "gene")
    
    # Safely compute logTPM (adding +1 to avoid log(0))
    median_col_n <- grep("_median", colnames(select_tissue_mat_n), value = TRUE)
    print(median_col_n)
    
    select_tissue_mat_n$logtpm <- log2(select_tissue_mat_n[[median_col_n]] + 1)
    
    select_tissue_mat_n$expr_level <- cut(
      select_tissue_mat_n[[median_col_n]],
      breaks = c(0, quantile(select_tissue_mat_n[[median_col_n]][select_tissue_mat_n[[median_col_n]] > 0], probs = seq(0,1,1/5))),
      labels = c("no", "1", "2", "3","4","5"),
      include.lowest = TRUE
    )
    
    summary_df1 <- select_tissue_mat_n %>%
      group_by(expr_level) %>%
      summarise(across(
        all_of(grep("mC_bin\\d+_avg$", colnames(select_tissue_mat_n), value = TRUE)),
        mean,
        na.rm = TRUE
      ), .groups = "drop") %>%
      melt(id.vars = "expr_level") %>%
      mutate(tissue = select_tissue, condition = "Normal")
    
    ### ---- Tumour ---- ###
    select_gene_exp <- gene_exp[, c(1, 2, grep(paste0("^", select_tissue, "_tcga"), names(gene_exp)))]
    select_gene_meth <- meth[, c(4, grep(paste0("^", select_tissue, "-Tumour_.*avg$"), names(meth)))]
    
    select_tissue_mat_t <- merge(select_gene_exp, select_gene_meth, by = "gene")
    
    median_col_t <- grep("_median", colnames(select_tissue_mat_t), value = TRUE)
    print(median_col_t)
    select_tissue_mat_t$logtpm <- log2(select_tissue_mat_t[[median_col_t]] + 1)
    
    select_tissue_mat_t$expr_level <- cut(
      select_tissue_mat_t[[median_col_t]],
      breaks = c(0, quantile(select_tissue_mat_t[[median_col_t]][select_tissue_mat_t[[median_col_t]] > 0], probs = seq(0,1,1/5))),
      labels = c("no", "1", "2", "3","4","5"),
      include.lowest = TRUE
    )
    
    summary_df2 <- select_tissue_mat_t %>%
      group_by(expr_level) %>%
      summarise(across(
        all_of(grep("mC_bin\\d+_avg$", colnames(select_tissue_mat_t), value = TRUE)),
        mean,
        na.rm = TRUE
      ), .groups = "drop") %>%
      melt(id.vars = "expr_level") %>%
      mutate(tissue = select_tissue, condition = "Tumour")
    
    # Combine both into main summary
    summary_df <- rbind(summary_df, summary_df1, summary_df2)
  }
  ## blood cell ##
  gene_exp <- read.table("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/gene_expr/rna_immune_cell_selected.tpm.tsv",sep="\t",header=TRUE,check.names = FALSE)
  select_tissues <- colnames(gene_exp)[-c(1:2)]
  for (select_tissue in select_tissues) {
    select_gene_exp <- gene_exp[, c(1, 2, grep(select_tissue, names(gene_exp)))]
    select_gene_meth <- meth[, c(4, grep(paste0("^", select_tissue, "_.*avg$"), names(meth)))]
    
    select_tissue_mat_n <- merge(select_gene_exp[,-1], select_gene_meth, by.x="Gene.name",by.y = "gene")
    
    # Safely compute logTPM (adding +1 to avoid log(0))
    median_col <- colnames(select_tissue_mat_n)[2]
    print(median_col)
    select_tissue_mat_n$logtpm <- log2(select_tissue_mat_n[[median_col]] + 1)
    
    select_tissue_mat_n$expr_level <- cut(
      select_tissue_mat_n[[median_col]],
      breaks = c(0, quantile(select_tissue_mat_n[[median_col]][select_tissue_mat_n[[median_col]] > 0], probs = seq(0,1,1/5))),
      labels = c("no", "1", "2", "3","4","5"),
      include.lowest = TRUE
    )
    
    summary_df1 <- select_tissue_mat_n %>%
      group_by(expr_level) %>%
      summarise(across(
        all_of(grep("mC_bin\\d+_avg$", colnames(select_tissue_mat_n), value = TRUE)),
        mean,
        na.rm = TRUE
      ), .groups = "drop") %>%
      melt(id.vars = "expr_level") %>%
      mutate(tissue = select_tissue, condition = "blood")
    summary_df <- rbind(summary_df, summary_df1)
  }
  summary_df$bin <- gsub(".*bin|_avg","",summary_df$variable) %>% as.numeric()
  
  p1 <- ggplot(summary_df %>% filter(condition!="blood"), aes(x = bin, y = value, color = expr_level)) +
    geom_line() +
    facet_grid(condition ~ tissue, scales = "free_y") +
    scale_color_manual(values=rev(brewer.pal(6,"Spectral")))+
    scale_x_continuous(breaks = c(0, 20, 40, 60),
                       labels = c("-10K","TSS","TES","10K")) +
    theme_bw() +
    # ylim(0,ifelse(max(summary_df$value)<90, max(summary_df$value)+5, max(summary_df$value))) +
    ylab(paste0(meth_marker, " group by exp")) 
  summary_df_blood <- summary_df %>% filter(condition=="blood") %>% as.data.frame()
  p2 <- ggplot(summary_df_blood, aes(x = bin, y = value, color = expr_level)) +
    geom_line() +
    facet_grid(. ~ tissue) +  # or facet_wrap(~ tissue)
    scale_color_manual(values = rev(brewer.pal(6,"Spectral"))) +
    scale_x_continuous(breaks = c(0, 20, 40, 60),
                       labels = c("-10K","TSS","TES","10K")) +
    theme_bw() +
    ylim(
      0,
      ifelse(max(summary_df_blood$value) < 90,
             max(summary_df_blood$value) + 5,
             max(summary_df_blood$value))
    ) +
    ylab(paste0(meth_marker, " group by exp"))
  p <- plot_grid(p1,p2,ncol = 1, rel_heights = c(2,1))
  ggsave(paste0("figs/fig5/meth_vs_expression_bin_quantile_freey", meth_marker, ".pdf"), p, width=15,height = 6)
  
}
