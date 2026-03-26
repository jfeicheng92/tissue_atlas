setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/atac_seq/")
#### Normal ####
for (prefix in c("mC.bgQ0.05-0.05.tgQ0.25.hypo_dmrs.bg_quant_modegroups.healthytop500","hmC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.healthytop500","umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.healthytop500")){
  file_list <- list.files(pattern=paste0(prefix,".[A-Z,a-z].*.bed"))
  tissue_order <- c("Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen","CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells")
  sel_N <- 500
  library(data.table)
  library(stringr)
  
  read_one <- function(f) {
    dt <- fread(f, col.names = c("chr","start","end","value"))
    # pull the sample/tissue name after "healthytop500." and before ".bed"
    tissue <- str_match(basename(f), "healthytop500\\.(.+)\\.bed$")[,2]
    # fallback: use the filename (sans .bed) if the pattern isn't present
    if (is.na(tissue)) tissue <- tools::file_path_sans_ext(basename(f))
    dt[, sample := tissue]
    dt
  }
  long_u <- rbindlist(lapply(file_list, read_one), use.names = TRUE, fill = TRUE) %>% unique(., by = c("chr","start","end","sample"))
  setDT(long_u)
  
  lib_size <- c(
    "CD8-positive,_alpha-beta_T_cell" = 78128,
    "CD4-positive,_alpha-beta_T_cell" = 79991,
    "naive_B_cell"                    = 94347,
    "natural_killer_cell"             = 165857,
    "pancreas"                        = 185268,
    "heart_right_ventricle"           = 189229,
    "heart_left_ventricle"            = 196990,
    "esophagus_mucosa"                = 203402,
    "cerebellum"                      = 204666,
    "body_of_pancreas"                = 208089,
    "ovary"                           = 208215,
    "transverse_colon"                = 212418,
    "sigmoid_colon"                   = 215019,
    "lung"                            = 215321,
    "stomach"                         = 220726,
    "kidney"                          = 220767,
    "spleen"                          = 222314,
    "prostate_gland"                  = 224455,
    "liver"                           = 227394,
    "breast_epithelium"               = 229249,
    "esophagus_muscularis_mucosa"     = 237434,
    "esophagus_squamous_epithelium"   = 238284
  )
  long_u[, lib_size := lib_size[sample]]
  
  stopifnot(!any(is.na(long_u$lib_size)))
  
  long_u[, value_cpm := value * 1e6 / lib_size]
  
  gmean <- exp(mean(log(unname(lib_size))))
  long_u[, value_rescaled := value * (gmean / lib_size)]
  
  wide <- pivot_wider(long_u %>% select(chr,start,end,sample,value_rescaled), names_from = sample, values_from = value_rescaled)
  
  
  meta <- fread(paste0(prefix,".bed")) %>% 
    select(chr,start,end,delta_quants,selected_tissue) %>%
    filter(selected_tissue %in% tissue_order) %>%
    mutate(selected_tissue=factor(selected_tissue,tissue_order))
  
  
  
  wide_annot <- inner_join(meta, wide, by = c("chr","start","end"))
  
  top_per_tissue <- wide_annot %>%
    group_by(selected_tissue) %>%
    slice_max(delta_quants, n = sel_N, with_ties = FALSE) %>%
    ungroup()
  
  
  
  col_map <- c(
    cerebellum                   = "Brain - Cerebellum",
    breast_epithelium            = "Breast",
    heart_left_ventricle         = "Heart - left ventricle",
    kidney                       = "Kidney",
    liver                        = "Liver",
    lung                         = "Lung",
    ovary                        = "Ovary",
    pancreas                     = "Pancreas",
    prostate_gland               = "Prostate",
    transverse_colon             = "Colon - transverse",
    stomach                      = "Stomach",
    esophagus_mucosa             = "Esophagus - mucosa",
    spleen                       = "Spleen",
    `CD4-positive,_alpha-beta_T_cell` = "CD4-T-cell",
    `CD8-positive,_alpha-beta_T_cell` = "CD8-T-cell",
    natural_killer_cell          = "NK-cells",
    naive_B_cell                 = "B-cells (naive)"
    
  )
  cols <- intersect(names(col_map), names(top_per_tissue))
  top_per_tissue <- top_per_tissue %>%
    filter(selected_tissue!="Brain")  #  cerebellum is a small part of brain
  mat <- top_per_tissue %>%
    select(names(col_map)) %>%
    rename(!!! setNames(names(col_map[cols]), col_map[cols])) %>%
    select(-c("Brain - Cerebellum"))%>% 
    as.matrix()
  
  
  
  rownames(mat) <- paste(top_per_tissue$selected_tissue,
                         top_per_tissue$chr,
                         top_per_tissue$start,
                         top_per_tissue$end, sep = ":")
  
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           scale ="row",angle_col =90, filename = paste0(prefix,".top",sel_N,".atac-seq.scaled.pdf"))
  
  top_per_tissue %>%
    group_by(selected_tissue) %>%
    summarise(
      across(-c(chr,start,end,delta_quants), ~ median(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    select(-selected_tissue) %>%
    select(names(col_map)) %>%
    rename(!!! setNames(names(col_map[cols]), col_map[cols])) %>%
    select(-c("Brain - Cerebellum"))%>%
    as.matrix() %>%
    pheatmap(.,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
             angle_col =90, cluster_rows = FALSE, cluster_cols = FALSE,scale="row", filename = paste0(prefix,".top",sel_N,".atac-seq.scaled_mean.pdf"))
}

new_folder <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig4/"
list_of_files <- list.files(pattern=".*bgQ0.05-0.05.*healthytop500.top500.atac-seq.scale.*pdf$") 
file.copy(file.path(list_of_files), new_folder)

#### Tumour ####
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/atac_seq/")
# multiBigwigSummary bins -b *.bw -o bw_summary.npz --outRawCounts bw_summary.tab
# for i in `seq 4 1 11`;do echo `cat  bw_summary.tab |head -1 |cut -f$i` `cat  bw_summary.tab |cut -f$i|tail -n +2|awk '{sum+=$0}END{print sum}'`;done
# 'BRCA.bw' 787264
# 'COAD.bw' 703116
# 'GBM.bw' 782975
# 'KIRC.bw' 696152
# 'LIHC.bw' 659520
# 'LUAD.bw' 800505
# 'PRAD.bw' 902279
# 'STAD.bw' 860841


for (prefix in c("mC.bgQ0.05-0.05.tgQ0.25.hypo_dmrs.bg_quant_modegroups.alltop500","hmC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop500","umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop500")){
  file_list <- list.files(pattern=paste0(prefix,".[A-Z,a-z].*.bed"))
  tissue_order <- c("Brain-Tumour","Breast-Tumour", "Colon-Tumour", "Kidney-Tumour","Liver-Tumour", "Lung-Tumour", "Prostate-Tumour","Stomach-Tumour")
  sel_N <- 500
  library(data.table)
  library(stringr)
  
  read_one <- function(f) {
    dt <- fread(f, col.names = c("chr","start","end","value"))
    # pull the sample/tissue name after "healthytop500." and before ".bed"
    tissue <- str_match(basename(f), "alltop500\\.(.+)\\.bed$")[,2]
    # fallback: use the filename (sans .bed) if the pattern isn't present
    if (is.na(tissue)) tissue <- tools::file_path_sans_ext(basename(f))
    dt[, sample := tissue]
    dt
  }
  
  long_u <- rbindlist(lapply(file_list, read_one), use.names = TRUE, fill = TRUE) %>% unique(., by = c("chr","start","end","sample"))
  setDT(long_u)
  
  lib_size <- c(
    BRCA = 787264, COAD = 703116, GBM  = 782975, KIRC = 696152,
    LIHC = 659520, LUAD = 800505, PRAD = 902279, STAD = 860841
  )
  
  long_u[, lib_size := lib_size[sample]]
  
  stopifnot(!any(is.na(long_u$lib_size)))
  
  long_u[, value_cpm := value * 1e6 / lib_size]
  
  gmean <- exp(mean(log(unname(lib_size))))
  long_u[, value_rescaled := value * (gmean / lib_size)]
  
  wide <- pivot_wider(long_u %>% select(chr,start,end,sample,value_rescaled), names_from = sample, values_from = value_rescaled)
  
  
  meta <- fread(paste0(prefix,".txt")) %>% 
    select(chr,start,end,delta_quants,selected_tissue) %>%
    filter(selected_tissue %in% tissue_order) %>%
    mutate(selected_tissue=factor(selected_tissue,tissue_order))
  
  
  
  wide_annot <- inner_join(meta, wide, by = c("chr","start","end"))
  
  top_per_tissue <- wide_annot %>%
    group_by(selected_tissue) %>%
    slice_max(delta_quants, n = sel_N, with_ties = FALSE) %>%
    ungroup()
  
  
  col_map <- c(
    "GBM"  = "Brain-Tumour(GMB)",
    "BRCA" = "Breast-Tumour(BRCA)",
    "COAD" = "Colon-Tumour(COAD)",
    "KIRC" = "Kidney-Tumour(KIRC)",
    "LIHC" = "Liver-Tumour(LIHC)",
    "LUAD" = "Lung-Tumour(LUAD)",
    "PRAD" = "Prostate-Tumour(PRAD)",
    "STAD" = "Stomach-Tumour(STAD)"
  )
  cols <- intersect(names(col_map), names(top_per_tissue))
  mat <- top_per_tissue %>%
    select(names(col_map)) %>%
    rename(!!! setNames(names(col_map[cols]), col_map[cols])) %>%
    as.matrix()
  
  
  
  rownames(mat) <- paste(top_per_tissue$selected_tissue,
                         top_per_tissue$chr,
                         top_per_tissue$start,
                         top_per_tissue$end, sep = ":")
  
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           scale ="row",angle_col =90, filename = paste0(prefix,".top",sel_N,".atac-seq.scaled.pdf"))
  
  top_per_tissue %>%
    group_by(selected_tissue) %>%
    summarise(
      across(-c(chr,start,end,delta_quants), ~ median(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    select(-selected_tissue) %>%
    select(names(col_map)) %>%
    rename(!!! setNames(names(col_map[cols]), col_map[cols])) %>%
    as.matrix() %>%
    pheatmap(.,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
             angle_col =90, cluster_rows = FALSE, cluster_cols = FALSE,scale="row", filename = paste0(prefix,".top",sel_N,".atac-seq.scaled_mean.pdf"))
  
}
new_folder <- "/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/figs/fig6/"
list_of_files <- list.files(pattern=".*bgQ0.05-0.05.*alltop500.top500.atac-seq.scaled.*pdf$") 
file.copy(file.path(list_of_files), new_folder)
