library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(cowplot)
library(reshape2)
library(stringr)
options(bitmapType='cairo-png')
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls")

prefix <- "dmr_call_grail/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500"
meth_mean <- fread(paste0(prefix,".mean.txt")) %>% as.data.frame()

meth_mean <- meth_mean %>%
  separate(V1, into = c("chr", "start", "end"), sep = "_", remove = FALSE) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end)
  )

#### tissue-specific-gene ####
#PLG (chr6:160,702,193-160,754,097), HS3ST1 (chr4:11,393,150-11,434,327) and CD300A (chr17:74,466,373-74,484,798).
library(forcats)
for(pos in c("chr6:160,702,193-160,754,097","chr4:11,393,150-11,434,327","chr17:74,466,373-74,484,798")){
  pos_chr <- gsub(":.*","",pos)
  pos_start <- gsub(".*:|,","",pos) %>% gsub("-.*","",.) %>% as.numeric()
  pos_end <- gsub(".*:|,","",pos) %>% gsub(".*-","",.) %>% as.numeric()
  sel_meth <- meth_mean %>%
    filter(chr==pos_chr & start > pos_start & end < pos_end)
  
  df_long <- sel_meth %>%
    select(starts_with("mean_")) %>%                 # e.g. mean_Liver_mC, ...
    mutate(sample = row_number()) %>%
    pivot_longer(cols = starts_with("mean_"),
                 names_to = "mod_full",
                 values_to = "value") %>%
    separate(mod_full, into = c("calc", "tissue", "mod"),
             sep = "_", remove = FALSE) %>%
    mutate(
      tissue = factor(tissue),                        # optional: control order
      mod    = factor(mod, levels = c("mC","hmC","umC"))
    )
  
  tissue_order <- c( "Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                     "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "CD34-erythroblasts", "CD34-megakaryocytes")
  df_long_sel <- df_long[df_long$tissue %in% tissue_order,]
  df_long_sel$tissue <- factor(df_long_sel$tissue, levels=tissue_order)

  
  p1 <- ggplot(df_long_sel, aes(x = factor(sample), y = value, fill = mod)) +
    geom_bar(stat = "identity", alpha = 0.6) +
    facet_grid(tissue ~ .) +
    scale_fill_manual(
      values = c("mC" = "#373C94", "hmC" = "#C03A2E", "umC" = "#76B7B2"),
      breaks = c("mC","hmC","umC"),
      labels = c("5mC","5hmC","umC")
    ) +
    labs(x = "Sample", y = "Proportion (%)", fill = "Modification") +
    theme_minimal() +
    theme(
      panel.spacing.y = unit(1, "mm"),   # ↓ space between facet rows
      panel.spacing.x = unit(1, "mm"),   # ↓ space between facet columns (if any)
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      strip.text.y = element_text(size = 9, margin = margin(0,0,0,0)), # slimmer strips
      strip.background = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
  df_area <- df_long_sel %>%
  mutate(sample_f = fct_inorder(factor(sample)),
         x = as.integer(sample_f))  
  p2 <- ggplot(df_area, aes(x = x, y = value, fill = mod)) +
  geom_area(position = "fill", color = "grey25", size = 0.2, alpha = 0.85) +
  facet_grid(tissue ~ .) +
  scale_fill_manual(values = c(mC="#373C94", hmC="#C03A2E", umC="#76B7B2"),
                    breaks = c("mC","hmC","umC"),
                    labels = c("5mC","5hmC","umC")) +
  scale_x_continuous(breaks = unique(df_area$x), labels = levels(df_area$sample_f)) +
  labs(x = "Sample", y = "Proportion", fill = "mod") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p3 <- ggplot(df_area, aes(x = x, y = value, fill = mod)) +
  geom_area(position = "fill", color = "grey25", size = 0.2, alpha = 0.85) +
  facet_grid(tissue ~ .) +
  scale_fill_manual(values = c(mC="#373C94", hmC="#C03A2E", umC="#76B7B2"),
                    breaks = c("mC","hmC","umC"),
                    labels = c("5mC","5hmC","umC")) +
  scale_x_continuous(breaks = unique(df_area$x), labels = levels(df_area$sample_f)) +
  labs(x = "Sample", y = "Proportion", fill = "mod") +
  theme_void() +                                 # no axes/grid/background
  theme(
    panel.spacing = grid::unit(0, "pt"),   # remove space between facet panels
    # or just the vertical spacing for facet_grid(rows ~ .):
    # panel.spacing.y = grid::unit(0, "pt"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)       # optional: no outer margin
  )
 pdf(paste0("figs/fig2/",pos,".area.pdf"), width = 10, height = 15)
 print(p3) 
 print(p2)
 print(p1)
 dev.off()
}

#### Overlap between select gene and DMR ####
# required packages
library(data.table)      # fread
library(dplyr)           # optional piping / select
library(GenomicRanges)   # GRanges, findOverlaps

# your inputs (example)
pos <- c(
  "chr6:160,702,193-160,754,097",
  "chr4:11,393,150-11,434,327",
  "chr17:74,466,373-74,484,798"
)
bg_tg_q <- "bgQ0.05-0.05.tgQ0.25"
dmr_files <- c(
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hypo_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.mC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt"),
  paste0("all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.",bg_tg_q,".hyper_dmrs.bg_quant_modegroups.all.txt")
)

# read one dmr file (you already did this; keep for completeness)
dmr <- fread(paste0("dmr_call_new1/", dmr_files[6])) %>%
  select(chr, start, end, selected_tissue) %>%
  as.data.frame()

# --- 1) parse pos into a data.frame of chr / start / end (remove commas) ---
pos_df <- data.frame(pos = pos, stringsAsFactors = FALSE) %>%
  # remove commas, then split
  mutate(pos_clean = gsub(",", "", pos)) %>%
  tidyr::separate(pos_clean, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  tidyr::separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  select(chr, start, end, pos)

# --- 2) make GRanges for pos regions ---
pos_gr <- GRanges(
  seqnames = pos_df$chr,
  ranges = IRanges(start = pos_df$start, end = pos_df$end),
  pos_label = pos_df$pos   # keep original label as metadata
)

# --- 3) make GRanges for DMRs (ensure numeric start/end) ---
# If your dmr columns are not numeric, convert them:
dmr$start <- as.integer(dmr$start)
dmr$end   <- as.integer(dmr$end)

dmr_gr <- GRanges(
  seqnames = dmr$chr,
  ranges = IRanges(start = dmr$start, end = dmr$end),
  selected_tissue = dmr$selected_tissue
)

# Optional: harmonize seqname style (e.g., "chr1" vs "1")
# seqlevelsStyle(pos_gr) <- seqlevelsStyle(dmr_gr)  # uncomment if needed

# --- 4) find overlaps ---
hits <- findOverlaps(pos_gr, dmr_gr, ignore.strand = TRUE)

# hits is a Hits object: queryHits -> index of pos_gr, subjectHits -> index of dmr_gr
if (length(hits) == 0) {
  message("No overlaps found.")
} else {
  # get a data.frame of overlaps with useful columns
  ov_df <- data.frame(
    pos_index = queryHits(hits),
    pos_label = mcols(pos_gr)$pos_label[queryHits(hits)],
    pos_chr   = as.character(seqnames(pos_gr))[queryHits(hits)],
    pos_start = start(pos_gr)[queryHits(hits)],
    pos_end   = end(pos_gr)[queryHits(hits)],
    dmr_index = subjectHits(hits),
    dmr_chr   = as.character(seqnames(dmr_gr))[subjectHits(hits)],
    dmr_start = start(dmr_gr)[subjectHits(hits)],
    dmr_end   = end(dmr_gr)[subjectHits(hits)],
    dmr_tissue = mcols(dmr_gr)$selected_tissue[subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  
  # if you want the original dmr rows joined:
  ov_joined <- cbind(
    ov_df,
    dmr[ov_df$dmr_index, , drop = FALSE]
  )
  
  # Results
  print(ov_df)
  # or view full joined DMR rows:
  # head(ov_joined)
}



library(forcats)

# chr1:1303567-1316677 PUSL1
# chr2:60933141-61023259 PUS10
# chr3:9832802-9849602 RPUSD3
# chr7:105434661-105527272 PUS7
# chr9:128300159-128327741 TRUB2
# chr10:114933195-114982676 TRUB1
# chr11:125888485-125908224 PUS3
# chr11:126197095-126216692 RPUSD4
# chr12:43713992-43764193 PUS7L
# chr12:131919277-131934646 PUS1-AS1
# chr12:131924200-131950896 PUS1
# chr15:40564299-40580171 RPUSD2
# chr16:779974-793406 RPUSD1
# chrX:154757742-154782697 DKC1
# "chrX:154757742-154782697;DKC1;+"
pos_list <- c("chr1:1303567-1316677;PUSL1;+","chr2:60933141-61023259;PUS10;-","chr3:9832802-9849602;RPUSD3;-","chr7:105434661-105527272;PUS7;-","chr9:128300159-128327741;TRUB2;-","chr10:114933195-114982676;TRUB1;+","chr11:125888485-125908224;PUS3;-","chr11:126197095-126216692;RPUSD4;-","chr12:43713992-43764193;PUS7L;-","chr12:131919277-131934646;PUS1-AS1;-","chr12:131924200-131950896;PUS1;+","chr15:40564299-40580171;RPUSD2;+","chr16:779974-793406;RPUSD1;-")
tissue_order <- c( "Brain","Brain-Tumour","Breast","Breast-Tumour", "Kidney","Kidney-Tumour","Liver","Liver-Tumour", "Lung","Lung-Tumour","Ovary","Ovary-Tumour", "Pancreas", "Pancreas-Tumour","Prostate","Prostate-Tumour","Colon","Colon-Tumour","Stomach","Stomach-Tumour")
pdf(paste0("pus_gene/pus_gene.meth.pdf"), width = 10, height = 10)
for(pos in pos_list){
  cat(pos)
  pos_chr <- gsub(":.*","",pos)
  pos_start <- gsub(".*:","",pos) %>% gsub(";.*","",.) %>% gsub("-[0-9].*","",.) %>% as.numeric()
  pos_end <- gsub(".*:","",pos) %>% gsub(";.*","",.) %>% gsub(".*[0-9]*-","",.)  %>% as.numeric()
  
  sel_meth <- meth_mean %>%
    filter(chr==pos_chr & start > pos_start & end < pos_end)
  
  df_long <- sel_meth %>%
    select(starts_with("mean_")) %>%                 # e.g. mean_Liver_mC, ...
    mutate(sample = row_number()) %>%
    pivot_longer(cols = starts_with("mean_"),
                 names_to = "mod_full",
                 values_to = "value") %>%
    separate(mod_full, into = c("calc", "tissue", "mod"),
             sep = "_", remove = FALSE) %>%
    mutate(
      tissue = factor(tissue),                        # optional: control order
      mod    = factor(mod, levels = c("mC","hmC","umC"))
    )
  
                      # "Esophagus","Spleen","Heart","CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "CD34-erythroblasts", "CD34-megakaryocytes")
  df_long_sel <- df_long[df_long$tissue %in% tissue_order,]
  df_long_sel$tissue <- factor(df_long_sel$tissue, levels=tissue_order)
  df_area <- df_long_sel %>%
    mutate(sample_f = fct_inorder(factor(sample)),
           x = as.integer(sample_f))  
  p2 <- ggplot(df_area, aes(x = x, y = value, fill = mod)) +
    geom_area(position = "fill", color = "grey25", size = 0.2, alpha = 0.85) +
    facet_grid(tissue ~ .) +
    scale_fill_manual(values = c(mC="#373C94", hmC="#C03A2E", umC="#76B7B2"),
                      breaks = c("mC","hmC","umC"),
                      labels = c("5mC","5hmC","umC")) +
    scale_x_continuous(breaks = unique(df_area$x), labels = levels(df_area$sample_f)) +
    labs(x = "Sample", y = "Proportion", fill = "mod") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(pos)
  print(p2)
}
dev.off()
