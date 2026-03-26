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
setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmr_call_new1/")

#### compare overall level in specific features ####
library(data.table)
library(tidyr)

# --- Load & prepare annotation (explode comma-separated features) ---
anno <- fread("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/resource/hg38_ws1000.s500.feature.bed",
              col.names = c("chr","start","end","feature"))
anno <- as_tibble(anno) |> separate_rows(feature, sep = "\\s*,\\s*")
setDT(anno); setkey(anno, chr, start, end)

# --- Read meth with only the columns you actually need ---
prefix <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500"
# prefix <- "all_sample.merged.mlml.mincov10_common.solo_WCGW.groupby.hg38_ws1000.s500"

# Grab header to pick *_mC columns (and optionally just selected tissues)
sel_tissues <- c("Brain","Breast","Colon","Kidney","Liver","Lung","Ovary","Pancreas","Prostate","Stomach")
meth_file <- paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/methregion/",prefix, ".bed")
hdr <- names(fread(meth_file, nrows = 0))
keep_mC <- grep(paste(c(paste(sel_tissues,"_",sep=""),
                        paste(sel_tissues,"-Tumour",sep="")),collapse="|"), hdr, value = TRUE)
keep_cols <- c("chr","start","end", keep_mC)
meth <- fread(meth_file, select = keep_cols)
setDT(meth); setkey(meth, chr, start, end)
meth[, rid := .I]  # row index for fast lookup

# --- Select features of interest ---
select_features <- c(
  "CGIshelves","CGIshore","cpgIsland",
  "GRCh38.Regulatory_Build.CTCF_binding_site","GRCh38.Regulatory_Build.enhancer",
  "GRCh38.Regulatory_Build.open_chromatin_region","GRCh38.Regulatory_Build.promoter",
  "GRCh38.Regulatory_Build.promoter_flanking_region","GRCh38.Regulatory_Build.TF_binding_site",
  "lincRNA","MANE.GRCh38.v1.0.refseq_genomic.gene","PMD_coordinates_hg38.commonPMD",
  "hg38.repeatmasker.DNA_hAT","hg38.repeatmasker.DNA_TcMar","hg38.repeatmasker.LINE_CR1",
  "hg38.repeatmasker.LINE_L1","hg38.repeatmasker.LINE_L2","hg38.repeatmasker.LINE_RTE",
  "hg38.repeatmasker.LTR_ERV","hg38.repeatmasker.SINE_Alu","hg38.repeatmasker.SINE_MIR"
)

anno_sel <- anno[feature %in% select_features, .(chr, start, end, feature)]

# --- Tiny join: map each (chr,start,end) to a meth row index (no big data copy) ---
# which=TRUE returns row indices in 'meth' for each row of 'anno_sel'
anno_sel[, rid := meth[.SD, on = .(chr, start, end), which = TRUE]]
idx <- anno_sel[!is.na(rid), .(rid), by = feature]  # just feature + row indices

# --- Compute means using indices (low peak RAM) ---
M <- as.matrix(meth[, ..keep_mC])  # one big matrix in memory, once

feature_meth_list <- lapply(split(idx$rid, idx$feature),
                            function(r) colMeans(M[r, , drop = FALSE], na.rm = TRUE))

# bind to a data.frame: one row per feature, columns = *_mC
feature_meth <- do.call(rbind, feature_meth_list)
feature_meth <- data.frame(feature = rownames(feature_meth), feature_meth, row.names = NULL)
colnames(feature_meth) <- gsub("\\.", "-", colnames(feature_meth))

feature_meth_summary <- feature_meth %>%
  pivot_longer(
    cols = keep_mC,
    names_to = c("tissue_raw","mark"),
    names_pattern = "^[^_]+_([^_]+)_(mC|hmC|umC)$",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  mutate(
    tumour = if_else(stringr::str_ends(tissue_raw, "-Tumour"), "yes", "no"),
    tissue = stringr::str_remove(tissue_raw, "-Tumour$"),
    mark   = dplyr::recode(mark, mC="5mC", hmC="5hmC", umC="5umC")
  ) %>%
  filter(tissue %in% sel_tissues) %>%
  group_by(feature, tissue, mark) %>%
  summarise(
    n_normal = sum(tumour=="no"),
    n_tumour = sum(tumour=="yes"),
    mean_normal = if (n_normal>0) mean(value[tumour=="no"]) else NA_real_,
    mean_tumour = if (n_tumour>0) mean(value[tumour=="yes"]) else NA_real_,
    delta = mean(value[tumour=="yes"]) - mean(value[tumour=="no"]),
    p_wilcox = { n_vals <- value[tumour=="no"]; t_vals <- value[tumour=="yes"]
    if (length(n_vals)>0 && length(t_vals)>0)
      wilcox.test(t_vals, n_vals, exact=FALSE)$p.value else NA_real_ },
    p_ttest = { n_vals <- value[tumour=="no"]; t_vals <- value[tumour=="yes"]
    .safe_ttest_p(t_vals, n_vals) },
    .groups = "drop"
  )

feature_meth_summary <- feature_meth_summary %>%
  mutate(sig = case_when(
    p_ttest < 0.001 ~ "***",
    p_ttest < 0.01  ~ "**",
    p_ttest < 0.05  ~ "*",
    TRUE ~ ""
  ))

feature_meth_summary <- feature_meth_summary %>%
  mutate(
    delta_lab = round(delta, 2),                # "0.12"
    sig_lab   = ifelse(is.na(sig) | sig == "", "", sig),       # "" or "*"
    label     = ifelse(sig_lab == "", delta_lab, paste0(delta_lab, "\n", sig_lab))
  )
feature_meth_summary$feature <- factor(feature_meth_summary$feature, levels = select_features)
feature_meth_summary <- feature_meth_summary[order(feature_meth_summary$feature),]
write.csv(feature_meth_summary, paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig6/",prefix,"tumour_normal_meth_delta_feature.txt"),quote=FALSE)


prefix <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500"
prefix <- "all_sample.merged.mlml.mincov10_common.solo_WCGW.groupby.hg38_ws1000.s500"
feature_meth_summary <- read.csv(paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig6/",prefix,"tumour_normal_meth_delta_feature.txt"))
feature_meth_summary <- feature_meth_summary[which(!feature_meth_summary$X %in% c("*","**","***")), ]
feature_meth_summary <- feature_meth_summary %>% mutate(label     = ifelse(sig_lab == "", delta_lab, paste0(delta_lab, "\n", sig_lab)))

L <- 10  
sel_feature <- "cpgIsland"
p1 <- ggplot(feature_meth_summary%>%filter(feature==sel_feature), aes(mark, tissue, fill = delta)) +
  geom_tile() +
  geom_text(aes(label = sig), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "RdBu", type = "div", direction = -1,            # <-- RdBu (blue low, red high)
    limits = c(-L, L),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p2 <- ggplot(feature_meth_summary%>%filter(feature==sel_feature), aes(mark, tissue, fill = delta)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "RdBu", type = "div", direction = -1,            # <-- RdBu (blue low, red high)
    limits = c(-L, L),
    oob = scales::squish,
    na.value = "grey90"
  )+
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  ggtitle(sel_feature)

L <- 30  
sel_feature <- "PMD_coordinates_hg38.commonPMD"
p3 <- ggplot(feature_meth_summary%>%filter(feature==sel_feature), aes(mark, tissue, fill = delta)) +
  geom_tile() +
  geom_text(aes(label = sig), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "RdBu", type = "div", direction = -1,            # <-- RdBu (blue low, red high)
    limits = c(-L, L),
    oob = scales::squish,
    na.value = "grey90"
  )+
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p4 <- ggplot(feature_meth_summary%>%filter(feature==sel_feature), aes(mark, tissue, fill = delta)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "RdBu", type = "div", direction = -1,            # <-- RdBu (blue low, red high)
    limits = c(-L, L),
    oob = scales::squish,
    na.value = "grey90"
  )+
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  ggtitle(sel_feature)

pdf(paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig6/",prefix,"tumour_normal_meth_delta_feature.sel.pdf"), width = 8, height = 5)
cowplot::plot_grid(p1,p2)
cowplot::plot_grid(p3,p4)
dev.off()




feature_meth_summary %>%
  select(feature, tissue, mark, mean_normal, mean_tumour) %>%
  pivot_longer(
    c(mean_normal, mean_tumour),
    names_to   = "group",
    names_prefix = "mean_",
    values_to  = "mean"
  ) %>%
  mutate(
    group  = recode(group, normal = "Normal", tumour = "Tumour"),
    # Optional: control orders if you have preferred ones
    tissue = factor(tissue, levels = sel_tissues),
    mark   = factor(mark, levels = c("5mC","5hmC","5umC"))
  ) %>%
  filter(feature %in% c("cpgIsland","PMD_coordinates_hg38.commonPMD")) %>%
  ggplot(., aes(x = tissue, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(feature ~mark, ncol=6) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +  # meth fractions
  labs(x = NULL, y = "Mean methylation", fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

##### bivalent promoter #####
both_feats <- c("cpgIsland", "hg38lift_genome_100_segments/BivProm")

# 1) Find (chr,start,end) present in BOTH features
both_bins <- anno %>%
  filter(feature %in% both_feats) %>%
  distinct(chr, start, end, feature) %>%
  group_by(chr, start, end) %>%
  summarise(n = n_distinct(feature), .groups = "drop") %>%
  filter(n == length(both_feats)) %>%
  select(chr, start, end)

# 2) Restrict meth to those bins, then summarise per tissue/mark
feature_meth_summary <- meth %>%
  # keep only bins annotated with BOTH features
  semi_join(both_bins, by = c("chr","start","end")) %>%
  # LONG once: parse sample/tissue and mark from column names
  pivot_longer(
    cols = matches("_(mC|hmC|umC)$"),
    names_to = c("tissue_raw","mark"),
    names_pattern = "^[^_]+_([^_]+)_(mC|hmC|umC)$",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  # split tumour flag and base tissue
  mutate(
    tumour = if_else(str_ends(tissue_raw, "-Tumour"), "yes", "no"),
    tissue = str_remove(tissue_raw, "-Tumour$"),
    mark   = recode(mark, mC = "5mC", hmC = "5hmC", umC = "5umC")
  ) %>%
  filter(tissue %in% sel_tissues) %>%
  group_by(tissue, mark) %>%
  summarise(
    n_normal   = sum(tumour == "no"),
    n_tumour   = sum(tumour == "yes"),
    mean_normal = if (n_normal > 0) mean(value[tumour == "no"], na.rm = TRUE) else NA_real_,
    mean_tumour = if (n_tumour > 0) mean(value[tumour == "yes"], na.rm = TRUE) else NA_real_,
    delta       = mean(value[tumour == "yes"], na.rm = TRUE) - mean(value[tumour == "no"], na.rm = TRUE),
    p_wilcox = {
      n_vals <- value[tumour == "no"]; t_vals <- value[tumour == "yes"]
      if (length(n_vals) > 0 && length(t_vals) > 0)
        suppressWarnings(stats::wilcox.test(t_vals, n_vals, exact = FALSE)$p.value)
      else NA_real_
    },
    p_ttest = {
      n_vals <- value[tumour == "no"]; t_vals <- value[tumour == "yes"]
      .safe_ttest_p(t_vals, n_vals)
    },
    .groups = "drop"
  ) %>%
  mutate(feature = "cpgIsland&BivProm") %>%
  relocate(feature, tissue, mark, n_normal, n_tumour, mean_normal, mean_tumour, delta)
feature_meth_summary <- feature_meth_summary %>%
  mutate(sig = case_when(
    p_ttest < 0.001 ~ "***",
    p_ttest < 0.01  ~ "**",
    p_ttest < 0.05  ~ "*",
    TRUE ~ ""
  ))

feature_meth_summary <- feature_meth_summary %>%
  mutate(
    delta_lab = round(delta, 2),                # "0.12"
    sig_lab   = ifelse(is.na(sig) | sig == "", "", sig),       # "" or "*"
    label     = ifelse(sig_lab == "", delta_lab, paste0(delta_lab, "\n", sig_lab))
  )
write.csv(feature_meth_summary,paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig6/",prefix,"tumour_normal_meth_delta_feature.cpgisland_bivprom.txt"),quote=FALSE)

feature_meth_summary <- read.csv(paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig6/",prefix,"tumour_normal_meth_delta_feature.cpgisland_bivprom.txt"))
feature_meth_summary <- feature_meth_summary[which(!feature_meth_summary$X %in% c("*","**","***")), ]
feature_meth_summary <- feature_meth_summary %>% mutate(label     = ifelse(sig_lab == "", delta_lab, paste0(delta_lab, "\n", sig_lab)))

L <- 10  # example: clamp to [-2, 2]
p1 <- ggplot(feature_meth_summary, aes(mark, tissue, fill = delta)) +
  geom_tile() +
  geom_text(aes(label = sig), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "RdBu", type = "div", direction = -1,            # <-- RdBu (blue low, red high)
    limits = c(-L, L),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p2 <- ggplot(feature_meth_summary, aes(mark, tissue, fill = delta)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "RdBu", type = "div", direction = -1,            # <-- RdBu (blue low, red high)
    limits = c(-L, L),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("CpGisland & BivProm")

pdf(paste0("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/figs/fig6/",prefix,"tumour_normal_meth_delta_feature.cpgisland_bivprom.pdf"), width = 8, height = 5)
cowplot::plot_grid(p1,p2)
dev.off()


# 
# feature_meth_summary %>%
#   select(feature, tissue, mark, mean_normal, mean_tumour) %>%
#   pivot_longer(
#     c(mean_normal, mean_tumour),
#     names_to   = "group",
#     names_prefix = "mean_",
#     values_to  = "mean"
#   ) %>%
#   mutate(
#     group  = recode(group, normal = "Normal", tumour = "Tumour"),
#     # Optional: control orders if you have preferred ones
#     tissue = factor(tissue, levels = sel_tissues),
#     mark   = factor(mark, levels = c("5mC","5hmC","5umC"))
#   ) %>%
#   filter(feature %in% c("cpgIsland&BivProm")) %>%
#   ggplot(., aes(x = tissue, y = mean, fill = group)) +
#   geom_col(position = position_dodge(width = 0.75), width = 0.7) +
#   facet_wrap(feature ~mark, ncol=6) +
#   scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +  # meth fractions
#   labs(x = NULL, y = "Mean methylation", fill = NULL) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     strip.text = element_text(size = 8),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )




