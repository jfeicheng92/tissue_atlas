### deeptools_plot
options(bitmapType='cairo-png')
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(forecast)
library(zoo)  
options(ggrepel.max.overlaps = Inf)


tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes")
tissue_pairs <- fread("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/histone/tissue_pairs.txt", header=FALSE, col.names = c("tissue","cell"))
tissue_pairs <- tissue_pairs[complete.cases(tissue_pairs) & tissue_pairs$tissue %in% tissue_order,]
tissue_pairs$tissue <- factor(tissue_pairs$tissue, levels=tissue_order)

setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmr_call_new1")
setwd("histone/")


plot_dmb_histone <- function(rname, fname, oname, selN = 200, tissue_pairs, seltissue = "None") {
  # read matrix and robustly parse header
  tmp <-  readLines(fname, n = 1) %>% strsplit('],"') %>% unlist()
  params <- data.frame(par=lapply(gsub('\\{|\\}|\\"|\\[|\\]',"", tmp), function(x)unlist(strsplit(x,":"))[[1]]) %>% unlist,
                       value=lapply(gsub('\\{|\\}|\\"|\\[|\\]',"", tmp), function(x)unlist(strsplit(x,":"))[[2]]) %>% unlist)
  n <- params[params$par=='sample_boundaries',]$value %>% strsplit(",",3) %>% unlist %>% as.numeric() %>% `[`(2)
  samples <- gsub(",_","_",params[params$par=='sample_labels',]$value) %>% strsplit(",",3) %>% unlist %>% as.character()
  
  mat <- data.table::fread(fname, skip = 1)
  
  # top regions per tissue
  region <- data.table::fread(rname) |>
    dplyr::group_by(selected_tissue) |>
    dplyr::slice_max(order_by = delta_quants, n = selN, with_ties = FALSE) |>
    dplyr::ungroup()
  
  region$pos <- paste0(region$chr, ":", region$start, "-", region$end)
  
  # keep V4 (pos label) and value columns
  mat <- merge(
    region |> dplyr::select(selected_tissue, pos),
    mat[, -c(1:3, 5, 6)],                 # keep V4 + value columns
    by.x = "pos", by.y = "V4"
  )
  
  # PDF for per-tissue raw/smoothed/normalized
  pdf(oname, width = 20, height = 4)
  
  mat_agg <- data.frame(matrix(ncol = 6, nrow = 0))
  
  for (i in tissue_pairs$tissue) {
    # summarise per bin across selected regions

    mat_mean <- mat |>
      dplyr::filter(selected_tissue == i) |>
      dplyr::select(dplyr::starts_with("V")) |>
      dplyr::summarise(dplyr::across(dplyr::everything(), median, na.rm = TRUE))
    
    # name columns as sample:bin
    colnames(mat_mean) <- paste0(rep(samples, each = n), ":", rep(seq_len(n), length(samples)))
    
    mat_mean <- tidyr::pivot_longer(mat_mean, dplyr::everything(),
                                    names_to = "idx", values_to = "level")
    
    mat_mean$pos   <- as.integer(sub(".*:", "", mat_mean$idx))
    mat_mean$label <- sub(":.*", "", mat_mean$idx)
    mat_mean$cell  <- sub(".*\\.", "", mat_mean$label)
    
    # optional filter by tissue-to-cell mapping
    if (seltissue != "None") {
      keep <- tissue_pairs$cell[tissue_pairs$tissue %in% unique(c(seltissue, i))]
      mat_mean <- mat_mean[Reduce(`|`, lapply(keep, function(k) grepl(k, mat_mean$cell, fixed = TRUE))), ]
    }
    
    mat_mean_sel <- mat_mean |>
      dplyr::select(level, pos, cell) |>
      dplyr::group_by(cell, pos) |>
      dplyr::summarise(mean_value = level, .groups = "drop")
    
    mat_mean_sel$target <- "2"
    mat_mean_sel$target[grep(tissue_pairs$cell[tissue_pairs$tissue == i], mat_mean_sel$cell)] <- "1"
    
    anno_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique(mat_mean_sel$cell)))
    
    p1 <- ggplot2::ggplot(mat_mean_sel, ggplot2::aes(x = pos, y = mean_value, color = cell)) +
      ggplot2::geom_line(ggplot2::aes(linetype = target)) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste0(i, " ", sub("\\.gz$", "", fname))) +
      ggplot2::scale_x_continuous(breaks = c(1, 150, 300), labels = c("-15K", "0", "15K")) +
      ggplot2::scale_color_manual(values = anno_colors) +
      ggplot2::ylab(paste0(sub("\\.bg.*", "", fname), "\n raw")) +
      ggplot2::theme(legend.position ="None" )
    
    # -------- smoothing (do NOT touch pos) --------
    ma_order <- 15L
    
    wide <- mat_mean_sel |>
      dplyr::select(pos, cell, mean_value) |>
      tidyr::pivot_wider(names_from = cell, values_from = mean_value) |>
      dplyr::arrange(pos)
    
    num_cols <- setdiff(names(wide), "pos")
    smoothed <- wide
    smoothed[num_cols] <- lapply(wide[num_cols], function(x) {
      v <- zoo::rollmean(x, k = ma_order, align = "center", fill = NA)
      edge <- ma_order %/% 2
      # pad edges with original values
      v[1:edge] <- x[1:edge]
      v[(length(x) - edge + 1):length(x)] <- x[(length(x) - edge + 1):length(x)]
      as.numeric(v)
    })
    
    mat_mean_sel_avg <- smoothed |>
      tidyr::pivot_longer(-pos, names_to = "cell", values_to = "mean_value") |>
      dplyr::mutate(
        target = ifelse(grepl(tissue_pairs$cell[tissue_pairs$tissue == i], cell), "1", "2"),
        seltissue = i,
        value = "average"
      )
    
    p2 <- ggplot2::ggplot(mat_mean_sel_avg, ggplot2::aes(x = pos, y = mean_value, color = cell)) +
      ggplot2::geom_line(ggplot2::aes(linetype = target)) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste0(i, " ", sub("\\.gz$", "", fname))) +
      ggplot2::scale_x_continuous(breaks = c(1, 150, 300), labels = c("-15K", "0", "15K")) +
      ggplot2::scale_color_manual(values = anno_colors) +
      ggplot2::ylab(paste0(sub("\\.bg.*", "", fname), "\n smooth average")) +
      ggplot2::theme(legend.position ="None" )
    
    # -------- normalize by median (per cell), then smooth --------
    norm <- smoothed
    med <- vapply(norm[num_cols], stats::median, numeric(1), na.rm = TRUE)
    for (j in num_cols) norm[[j]] <- norm[[j]] - med[[j]]
    
    norm_sm <- norm

    mat_mean_sel_avg_nor <- norm_sm |>
      tidyr::pivot_longer(-pos, names_to = "cell", values_to = "mean_value") |>
      dplyr::mutate(
        target = ifelse(grepl(tissue_pairs$cell[tissue_pairs$tissue == i], cell), "1", "2"),
        seltissue = i,
        value = "nor"
      )
    
    p3 <- ggplot2::ggplot(mat_mean_sel_avg_nor, ggplot2::aes(x = pos, y = mean_value, color = cell)) +
      ggplot2::geom_line(ggplot2::aes(linetype = target)) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste0(i, " ", sub("\\.gz$", "", fname))) +
      ggplot2::scale_x_continuous(breaks = c(1, 150, 300), labels = c("-15K", "0", "15K")) +
      ggplot2::scale_color_manual(values = anno_colors) +
      ggplot2::ylab(paste0(sub("\\.bg.*", "", fname), "\n normalized by median"))
    
    print(cowplot::plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1, 1, 2)))
    
    mat_agg <- rbind(mat_agg, mat_mean_sel_avg, mat_mean_sel_avg_nor)
  }
  
  dev.off()
  
  utils::write.table(mat_agg, sub("\\.pdf$", ".agg.txt", oname), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Aggregate line + heat
  p4 <- mat_agg |>
    dplyr::filter(target == 1, value == "average") |>
    ggplot2::ggplot(ggplot2::aes(x = pos, y = mean_value, color = cell)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(sub("\\.gz$", "", fname)) +
    ggplot2::scale_x_continuous(breaks = c(1, 150, 300), labels = c("-15K", "0", "15K")) +
    ggplot2::ylab(fname)
  
  mat_agg_wide <- mat_agg |>
    dplyr::filter(target == 1, value == "average") |>
    dplyr::select(pos, cell, mean_value) |>
    tidyr::pivot_wider(names_from = pos, values_from = mean_value)
  
  long_heat <- cbind(mat_agg_wide["cell"],
                     t(apply(as.matrix(mat_agg_wide[,-1]), 1, scale))) |>
    as.data.frame() |>
    tidyr::pivot_longer(-cell, names_to = "variable", values_to = "value")
  
  # ensure y breaks match discrete labels
  long_heat$variable <- as.character(long_heat$variable)
  
  p5 <- cbind(mat_agg_wide[,1],mat_agg_wide[,-1] %>% apply(., 1, scale) %>% t()) %>%
    melt(id.vars=c("cell")) %>%
    ggplot(aes(x = cell, y = variable)) +
    geom_tile(aes(fill = value)) +
    coord_flip() +
    theme(legend.position = "right") +
    xlab("Tissue") + ylab("") +
    theme(axis.text=element_text(size=10)) +
    scale_fill_gradientn(
      colors = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100),
      limits = c(-4, 4),
      oob = scales::squish,               # clamp values outside [-4,4]
      breaks = seq(-4, 4, by = 2)         # optional: nicer legend ticks
    ) +
    scale_y_discrete(breaks=c("1","150","300"),labels = c("-15K","0","15K")) +  # if needed
    theme_minimal()
  
  
  pdf(sub("\\.pdf$", ".agg.pdf", oname), width = 20, height = 4)
  print(cowplot::plot_grid(p4, p5, nrow = 1, rel_widths = c(3, 1)))
  dev.off()
}



for(dmb_histone in list.files(pattern = "histone.*.top500.mat.gz")){
  plot_dmb_histone(rname=paste0(gsub("histone_[^_]*_|.mat.gz","",dmb_histone),".txt"),
                   fname=dmb_histone,
                   selN=500,
                   oname=gsub(".gz","median.pdf",dmb_histone),
                   tissue_pairs=tissue_pairs)
}






for(fname in list.files(pattern=".*agg.txt")){
  mat_agg <- fread(fname)
  anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(mat_agg$cell))))
  cell_order <- sapply(tissue_pairs$cell , function(p) {grep(patt=p, unique(mat_agg$cell))} ) %>% unlist() %>% unique(mat_agg$cell)[.]
  mat_agg$cell <- factor(mat_agg$cell, levels=cell_order)
  mat_agg$seltissue <- factor(mat_agg$seltissue,tissue_order)
  p <- mat_agg%>%
    filter(!seltissue %in% c("Erythroblasts","Megakaryocytes")) %>%
    filter(!cell%in%c("heart_left_ventricle","sigmoid_colon","lung")) %>%
    filter(value=="nor") %>%
    ggplot(aes(x=pos,y=mean_value,color=cell)) + 
    geom_line(aes(linetype=as.factor(target))) +
    theme_classic() + 
    facet_wrap( ~ seltissue, nrow = 5, scales="free")+
    scale_x_continuous(breaks = c(1,150,300),labels = c("-15K","0","15K"))+
    scale_color_manual(values=anno_colors)  +
    theme(legend.position = "bottom") +
    ylab(gsub("bg.*.txt","",fname)) 
  ggsave(gsub(".txt",".nor.all_tissue.pdf",fname),p,width = 8, height = 10)
}



for(fname in list.files(pattern=".*agg.txt")){
  mat_agg <- fread(fname)
  anno_colors<- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(mat_agg$cell))))
  cell_order <- sapply(tissue_pairs$cell , function(p) {grep(patt=p, unique(mat_agg$cell))} ) %>% unlist() %>% unique(mat_agg$cell)[.] 
  mat_agg$cell <- factor(mat_agg$cell, levels=cell_order)
  mat_agg$seltissue <- factor(mat_agg$seltissue,tissue_order)
  p4 <- mat_agg%>%
    filter(!seltissue %in% c("Erythroblasts","Megakaryocytes")) %>%
    filter(!cell%in%c("heart_left_ventricle","sigmoid_colon","lung")) %>%
    filter(target==1) %>%
    filter(value=="nor") %>%
    ggplot(aes(x=pos,y=mean_value,color=cell)) + 
    geom_line() +
    theme_bw() + 
    ggtitle(gsub(".txt","",fname)) +
    scale_x_continuous(breaks = c(1,150,300),labels = c("-15K","0","15K")) +
    ylab(gsub("bg.*.txt","",fname)) +
    scale_color_manual(values=anno_colors)
  
  mat_agg_wide <- mat_agg%>%
    filter(!seltissue %in% c("Erythroblasts","Megakaryocytes")) %>%
    filter(!cell%in%c("heart_left_ventricle","sigmoid_colon","lung")) %>%
    filter(target==1) %>%
    filter(value=="average") %>%
    select(pos,cell,mean_value) %>%
    pivot_wider(names_from = pos,values_from = mean_value)
  mat_agg_wide$cell <- factor(mat_agg_wide$cell,rev(cell_order))
  mat_agg_wide <- mat_agg_wide[order(mat_agg_wide$cell),]
  
  
  p5 <- cbind(mat_agg_wide[,1],mat_agg_wide[,-1] %>% apply(., 1, scale) %>% t()) %>%
    melt(id.vars=c("cell")) %>%
    ggplot(aes(x = cell, y = variable)) +
    geom_tile(aes(fill = value)) +
    coord_flip() +
    theme(legend.position = "right") +
    xlab("Tissue") + ylab("") +
    theme(axis.text=element_text(size=10)) +
    scale_fill_gradientn(
      colors = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
      limits = c(-3, 3),
      oob = scales::squish,               # clamp values outside [-4,4]
      breaks = seq(-3, 3, by = 2)         # optional: nicer legend ticks
    ) +
    scale_y_discrete(breaks=c(1,150,300),labels = c("-15K","0","15K")) +
    theme_minimal()
  
  pdf(gsub(".txt",".reorder.pdf",fname), width=12, height = 4)
  print(cowplot::plot_grid(p4,p5,nrow=1,rel_widths = c(1,0.6)))
  dev.off()
  
}


new_folder <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig4/"
list_of_files <- list.files(pattern=".pdf$") 
file.copy(file.path(list_of_files), new_folder)


