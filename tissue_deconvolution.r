library(forcats)
library(scales)
library(ggpubr)
setwd("/gpfs2/well/ludwig/users/cfo155/tissueMap/methcalls/")
# dir.create("figs/fig7")

#### compare TAPS CAPS on paired data ####
infile_taps <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.all.txt.all_CpG.rC.TAPS.20snprm.selN100.tissue_contribution.txt"
cfDNA_TAPS <- fread(paste0("cfDNA_deconvolution/",infile_taps)) %>% as.data.frame()
infile_caps <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.all.txt.all_CpG.rC.CAPS.20snprm.selN100.tissue_contribution.txt"
cfDNA_CAPS<- fread(paste0("cfDNA_deconvolution/",infile_caps)) %>% as.data.frame()

cfDNA_TAPS <- cfDNA_TAPS %>% 
  pivot_longer(
    cols = contains("_"),
    values_to = "value"
  )%>%
  mutate(
    group=str_remove(name,"_.*")
  )
cfDNA_TAPS_mean <- cfDNA_TAPS %>%
  group_by(tissue,group) %>%
  summarise(across(value, mean, na.rm = TRUE)) %>% as.data.frame()


cfDNA_CAPS <- cfDNA_CAPS %>% 
  pivot_longer(
    cols = contains("_"),
    values_to = "value"
  )%>%
  mutate(
    group=str_remove(name,"_.*")
  )
cfDNA_CAPS_mean <- cfDNA_CAPS %>%
  group_by(tissue,group) %>%
  summarise(across(value, mean, na.rm = TRUE)) %>% as.data.frame()

tissue_order <- cfDNA_TAPS_mean %>%
  filter(group=="healthycontrol") %>%
  arrange(value) %>%
  pull(tissue) %>%
  rev()

cfDNA_mean <- merge(cfDNA_TAPS_mean, cfDNA_CAPS_mean, by = c("tissue","group")) %>%
  dplyr::rename(mC = value.x, hmC = value.y)

##### healthy control #####
grp   <- "healthycontrol"
top_k <- 6

# union of top_k tissues from taps and caps (within the chosen group)
tissue_select <- bind_rows(
  cfDNA_mean %>% filter(group == grp) %>% slice_max(order_by = mC, n = top_k, with_ties = FALSE) %>% select(tissue),
  cfDNA_mean %>% filter(group == grp) %>% slice_max(order_by = hmC, n = top_k, with_ties = FALSE) %>% select(tissue)
) %>% distinct() %>% pull(tissue)

pie_df <- cfDNA_mean %>%
  filter(group == grp) %>%
  pivot_longer(cols = c(mC, hmC), names_to = "variable", values_to = "value") %>%
  mutate(tissue = if_else(tissue %in% tissue_select, tissue, "Other")) %>%
  group_by(variable, tissue) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(variable) %>%
  mutate(
    pct   = value / sum(value),
    label = if_else(pct >= 0.02, percent(pct, accuracy = 0.1), "")
  ) %>%
  ungroup()

# consistent legend/order across both pies; put "Other" last
level_order <- pie_df %>%
  group_by(tissue) %>%
  summarise(total = sum(value), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(tissue)

pie_df <- pie_df %>%
  mutate(tissue = factor(tissue, levels = c(level_order[level_order != "Other"], "Other")),
         variable=factor(variable, levels=c("mC","hmC")))

# plot: two pies side-by-side
lvls <- levels(pie_df$tissue)
base <- colorspace::qualitative_hcl(length(lvls) - 1, palette = "Dark 3")
pal  <- setNames(c(base, "#BDBDBD"), c(lvls[lvls != "Other"], "Other"))

pdf("figs/fig7/top6_tissue_contribution.heatlhy_control.pdf",width = 8, height = 5)
p <- ggplot(pie_df, aes(x = "", y = value, fill = tissue)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ variable, nrow = 1) +
  theme_void() +
  labs(
    title = paste0("Union of top ", top_k, " tissues across mC & hmC in ", grp),
    fill  = "Tissue"
  ) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = pal, breaks = lvls)
print(p)
dev.off()

cfDNA_merged <- merge(cfDNA_TAPS,cfDNA_CAPS,by=c("tissue","name","group")) %>%
  dplyr::rename(taps = value.x, caps = value.y)
p1 <- cfDNA_merged %>%
  filter(group=="healthycontrol") %>%
  ggplot(aes(x=taps,y=caps,color=tissue)) +
  geom_point() +
  theme_classic()

p2 <- cfDNA_merged %>%
  group_by(group,name) %>%
  summarise(
    n            = sum(complete.cases(taps, caps)),
    cor_pearson  = cor(taps, caps, use = "complete.obs"),
    cor_spearman = cor(taps, caps, method = "spearman", use = "complete.obs"),
    p_value      = cor.test(taps, caps, method = "pearson")$p.value
  ) %>%
  ungroup() %>%
  filter(group=="healthycontrol") %>%
  ggplot(aes(x=group,y=cor_pearson)) +
  geom_violin() +
  geom_jitter(width = 0.15,height = 0) +
  ylim(0,1) +
  theme_classic()
p <- cowplot::plot_grid(p1,p2,rel_widths = c(6,1))

pdf("figs/fig7/correlation_tissue_contribution.heatlhy_control.pdf",width = 9, height = 4)
print(p)
dev.off()



cfDNA_merged %>%
  group_by(group,name) %>%
  summarise(
    n            = sum(complete.cases(taps, caps)),
    cor_pearson  = cor(taps, caps, use = "complete.obs"),
    cor_spearman = cor(taps, caps, method = "spearman", use = "complete.obs"),
    p_value      = cor.test(taps, caps, method = "pearson")$p.value
  ) %>%
  ungroup() %>%
  group_by(group) %>%
  summarise(
    mean_cor_pearson            = mean(cor_pearson),
  ) 

# # A tibble: 4 × 2
# group          mean_cor_pearson
# <chr>                     <dbl>
#   1 HCC                       0.751
# 2 PDAC                      0.776
# 3 healthycontrol            0.769
# 4 livercontrol              0.751



##### PDAC ####


# helper for mean CI (Welch-style per group)
mean_ci <- function(x, conf = 0.95){
  x <- x[is.finite(x)]
  n <- length(x); m <- mean(x); s <- sd(x)
  if(n <= 1 || is.na(s)) return(data.frame(y=m, ymin=m, ymax=m))
  se <- s/sqrt(n); tcrit <- qt((1+conf)/2, df = n-1)
  data.frame(y = m, ymin = m - tcrit*se, ymax = m + tcrit*se)
}


p1 <- cfDNA_TAPS %>%
  filter(group %in% c("PDAC","healthycontrol") & tissue=="Liver") %>%
  ggplot(aes(x = group, y = value*100)) +
  # jittered raw data (zeros included)
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  # mean ± CI
  stat_summary(fun.data = mean_ci, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "two.sided"),
    comparisons = list(c("healthycontrol", "PDAC")),
    label = "p.format",   # or "p.signif"
    hide.ns = TRUE,
    label.y = max(cfDNA_TAPS$value*100, na.rm = TRUE) * 1.08
  )+
  expand_limits(y = max(cfDNA_TAPS$value*100, na.rm = TRUE)*1.1) +
  theme_classic() +
  ylab("Liver contribution") +
  ggtitle("TAPS decon") +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(1,3)])

p2 <-  cfDNA_CAPS %>%
  filter(group %in% c("PDAC","healthycontrol") & tissue=="Liver") %>%
  ggplot(aes(x = group, y = value*100)) +
  # jittered raw data (zeros included)
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  # mean ± CI
  stat_summary(fun.data = mean_ci, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "two.sided"),
    comparisons = list(c("healthycontrol", "PDAC")),
    label = "p.format",   # or "p.signif"
    hide.ns = TRUE,
    label.y = max(cfDNA_CAPS$value*100, na.rm = TRUE) * 1.08
  )+
  expand_limits(y = max(cfDNA_CAPS$value*100, na.rm = TRUE)*1.1) +
  theme_classic() +
  ylab("Liver contribution") +
  ggtitle("CAPS decon") +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(1,3)])

pdf("figs/fig7/liver_contribution.PDAC_healthycontrol.pdf",width = 7, height = 3)
cowplot::plot_grid(p1,p2)
dev.off() 

p1 <- cfDNA_TAPS %>%
  filter(group %in% c("PDAC","healthycontrol") & tissue=="Pancreas-Tumour") %>%
  ggplot(aes(x = group, y = value*100)) +
  # jittered raw data (zeros included)
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  # mean ± CI
  stat_summary(fun.data = mean_ci, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "two.sided"),
    comparisons = list(c("healthycontrol", "PDAC")),
    label = "p.format",   # or "p.signif"
    hide.ns = TRUE,
    label.y = max(cfDNA_TAPS$value*100, na.rm = TRUE) * 1.08
  )+
  expand_limits(y = max(cfDNA_TAPS$value*100, na.rm = TRUE)*1.1) +
  theme_classic() +
  ylab("Pancreas-Tumour contribution") +
  ggtitle("TAPS decon") +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(1,3)])

p2 <-  cfDNA_CAPS %>%
  filter(group %in% c("PDAC","healthycontrol") & tissue=="Pancreas-Tumour") %>%
  ggplot(aes(x = group, y = value*100)) +
  # jittered raw data (zeros included)
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  # mean ± CI
  stat_summary(fun.data = mean_ci, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "two.sided"),
    comparisons = list(c("healthycontrol", "PDAC")),
    label = "p.format",   # or "p.signif"
    hide.ns = TRUE,
    label.y = max(cfDNA_CAPS$value*100, na.rm = TRUE) * 1.08
  )+
  expand_limits(y = max(cfDNA_CAPS$value*100, na.rm = TRUE)*1.1) +
  theme_classic() +
  ylab("Pancreas-Tumour contribution") +
  ggtitle("CAPS decon") +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(1,3)])

pdf("figs/fig7/Pancreas_Tumour_contribution.PDAC_healthycontrol.pdf",width = 7, height = 3)
cowplot::plot_grid(p1,p2)
dev.off() 


p1 <- cfDNA_TAPS %>%
  filter(group %in% c("PDAC","healthycontrol") & tissue=="Liver-Tumour") %>%
  ggplot(aes(x = group, y = value*100)) +
  # jittered raw data (zeros included)
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  # mean ± CI
  stat_summary(fun.data = mean_ci, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "two.sided"),
    comparisons = list(c("healthycontrol", "PDAC")),
    label = "p.format",   # or "p.signif"
    hide.ns = TRUE,
    label.y = max(cfDNA_TAPS$value*100, na.rm = TRUE) * 1.08
  )+
  expand_limits(y = max(cfDNA_TAPS$value*100, na.rm = TRUE)*1.1) +
  theme_classic() +
  ylab("Liver-Tumour contribution") +
  ggtitle("TAPS decon") +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(1,3)])

p2 <-  cfDNA_CAPS %>%
  filter(group %in% c("PDAC","healthycontrol") & tissue=="Liver-Tumour") %>%
  ggplot(aes(x = group, y = value*100)) +
  # jittered raw data (zeros included)
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  # mean ± CI
  stat_summary(fun.data = mean_ci, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "two.sided"),
    comparisons = list(c("healthycontrol", "PDAC")),
    label = "p.format",   # or "p.signif"
    hide.ns = TRUE,
    label.y = max(cfDNA_CAPS$value*100, na.rm = TRUE) * 1.08
  )+
  expand_limits(y = max(cfDNA_CAPS$value*100, na.rm = TRUE)*1.1) +
  theme_classic() +
  ylab("Liver-Tumour contribution") +
  ggtitle("CAPS decon") +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(1,3)])

pdf("figs/fig7/Liver-Tumour_contribution.PDAC_healthycontrol.pdf",width = 7, height = 3)
cowplot::plot_grid(p1,p2)
dev.off() 

t.test(cfDNA_TAPS%>%filter(group=="healthycontrol"&tissue=="Liver") %>% pull(value),
       cfDNA_TAPS%>%filter(group=="PDAC"&tissue=="Liver") %>% pull(value), alternative = "less")
# Welch Two Sample t-test
# 
# data:  cfDNA_TAPS %>% filter(group == "healthycontrol" & tissue == "Liver") %>% pull(value) and cfDNA_TAPS %>% filter(group == "PDAC" & tissue == "Liver") %>% pull(value)
# t = -1.944, df = 16.363, p-value = 0.03465
# alternative hypothesis: true difference in means is less than 0
# 95 percent confidence interval:
#   -Inf -0.006400113
# sample estimates:
#   mean of x mean of y 
# 0.0372400 0.0992875 
t.test(cfDNA_CAPS%>%filter(group=="healthycontrol"&tissue=="Liver") %>% pull(value),
       cfDNA_CAPS%>%filter(group=="PDAC"&tissue=="Liver") %>% pull(value), alternative = "less")
# Welch Two Sample t-test
# 
# data:  cfDNA_CAPS %>% filter(group == "healthycontrol" & tissue == "Liver") %>% pull(value) and cfDNA_CAPS %>% filter(group == "PDAC" & tissue == "Liver") %>% pull(value)
# t = -1.9221, df = 20.276, p-value = 0.03438
# alternative hypothesis: true difference in means is less than 0
# 95 percent confidence interval:
#   -Inf -0.003982942
# sample estimates:
#   mean of x mean of y 
# 0.0160800 0.0546375 


##### HCC #####
# prepare filtered data
library(dplyr)
library(ggplot2)
library(ggpubr)        # stat_compare_means()
library(RColorBrewer)

# prepare filtered data
y_max <- 0.2*100

p3 <- cfDNA_TAPS[grep("NH125|NH134|NH140|NH159|NH160",cfDNA_TAPS$name,invert = TRUE),] %>%
# p3 <- cfDNA_TAPS %>%
  filter(group %in% c("healthycontrol","livercontrol", "HCC"),
         tissue == "Liver-Tumour") %>%
  mutate(group = factor(group, levels = c("healthycontrol","livercontrol", "HCC"))) %>%
  ggplot(aes(x = group, y = value*100)) +
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  stat_summary(fun.data = mean_cl_normal, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # add each comparison separately with its own y position
  stat_compare_means(method = "wilcox.test",
                     method.args = list(alternative = "two.sided"),
                     comparisons = list(c("livercontrol", "HCC")),
                     label = "p.format",
                     hide.ns = TRUE,
                     label.y = y_max * 1.02) +
  
  stat_compare_means(method = "wilcox.test",
                     method.args = list(alternative = "two.sided"),
                     comparisons = list(c("livercontrol", "healthycontrol")),
                     label = "p.format",
                     hide.ns = TRUE,
                     label.y = y_max * 1.10) +
  
  stat_compare_means(method = "wilcox.test",
                     method.args = list(alternative = "two.sided"),
                     comparisons = list(c("HCC", "healthycontrol")),
                     label = "p.format",
                     hide.ns = TRUE,
                     label.y = y_max * 1.18) +
  
  expand_limits(y = y_max * 1.25) +
  theme_classic() +
  ylab("Liver-Tumour contribution") +
  ggtitle("TAPS decon") +
  scale_fill_manual(values = brewer.pal(3, "Dark2"))

# print
p3

y_max <- 0.15*100
p4 <-  cfDNA_CAPS[grep("NH125|NH134|NH140|NH159|NH160",cfDNA_CAPS$name,invert = TRUE),]  %>%
# p4 <-  cfDNA_CAPS %>%
  filter(group %in% c("healthycontrol","livercontrol", "HCC"),
         tissue == "Liver-Tumour") %>%
  mutate(group = factor(group, levels = c("healthycontrol","livercontrol", "HCC"))) %>%
  ggplot(aes(x = group, y = value*100)) +
  geom_jitter(aes(fill = group),
              width = 0.08, height = 0,
              shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
  stat_summary(fun.data = mean_cl_normal, size = 0.8,
               geom = "errorbar", width = 0.12) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # add each comparison separately with its own y position
  stat_compare_means(method = "wilcox.test",
                     method.args = list(alternative = "two.sided"),
                     comparisons = list(c("livercontrol", "HCC")),
                     label = "p.format",
                     hide.ns = TRUE,
                     label.y = y_max * 1.02) +
  
  stat_compare_means(method = "wilcox.test",
                     method.args = list(alternative = "two.sided"),
                     comparisons = list(c("livercontrol", "healthycontrol")),
                     label = "p.format",
                     hide.ns = TRUE,
                     label.y = y_max * 1.10) +
  
  stat_compare_means(method = "wilcox.test",
                     method.args = list(alternative = "two.sided"),
                     comparisons = list(c("HCC", "healthycontrol")),
                     label = "p.format",
                     hide.ns = TRUE,
                     label.y = y_max * 1.18) +
  
  expand_limits(y = y_max * 1.25) +
  theme_classic() +
  ylab("Liver-Tumour contribution") +
  ggtitle("CAPS decon") +
  scale_fill_manual(values = brewer.pal(3, "Dark2"))


pdf("figs/fig7/livertumour_contribution.HCC_cirrhosis_excludetreated.pdf",width = 9, height = 3)
cowplot::plot_grid(p3,p4)
dev.off() 






  
  
  
# #### 2025 cfTAPS NC ####
# infile_taps <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.all.txt.cfTAPS_NC_2025.filtered.CpG.rC.selN200.tissue_contribution.txt"
# cfDNA_SCAN <- fread(paste0("cfDNA_deconvolution/",infile_taps)) %>% as.data.frame()
# 
# cfDNA_SCAN <- cfDNA_SCAN %>% 
#   pivot_longer(
#     cols = contains("_"),
#     values_to = "value"
#   )%>%
#   mutate(
#     group=str_remove(name,"_S[0-9].*")
#   )
# 
# 
# df <- cfDNA_SCAN %>%
#   filter(
#     tissue == "Colon-Tumour",
#     # tissue == "Colon",
#     group %in% c("CBHC_stage_Non-cancer", "SCAN_stage_Non-cancer",
#                  "CRC_stage_3", "CRC_stage_4")
#   ) %>%
#   mutate(cohort = if_else(group %in% c("CBHC_stage_Non-cancer","SCAN_stage_Non-cancer"),
#                           "Non-cancer", "CRC stage 3/4")) %>%
#   mutate(cohort = factor(cohort, levels=c("Non-cancer", "CRC stage 3/4")))
# 
# 
# 
# # helper for mean CI (Welch-style per group)
# mean_ci <- function(x, conf = 0.95){
#   x <- x[is.finite(x)]
#   n <- length(x); m <- mean(x); s <- sd(x)
#   if(n <= 1 || is.na(s)) return(data.frame(y=m, ymin=m, ymax=m))
#   se <- s/sqrt(n); tcrit <- qt((1+conf)/2, df = n-1)
#   data.frame(y = m, ymin = m - tcrit*se, ymax = m + tcrit*se)
# }
# 
# pval <- t.test(value ~ cohort, data = df, alternative = "less")$p.value
# 
# p5 <- ggplot(df, aes(x = cohort, y = value)) +
#   # jittered raw data (zeros included)
#   geom_jitter(aes(fill = cohort),
#               width = 0.08, height = 0,
#               shape = 21, color = "grey60", stroke = 0.4, size = 1.8, alpha = 0.7) +
#   # mean ± CI
#   stat_summary(fun.data = mean_ci, size = 0.8,
#                geom = "errorbar", width = 0.12) +
#   stat_summary(fun = mean, geom = "point", size = 3) +
#   stat_compare_means(
#     method = "t.test",
#     method.args = list(alternative = "two.sided"),
#     comparisons = list(c("Non-cancer", "CRC stage 3/4")),
#     label = "p.format",   # or "p.signif"
#     hide.ns = TRUE,
#     label.y = max(df$value, na.rm = TRUE) * 1.08
#   ) +
#   expand_limits(y = max(df$value, na.rm = TRUE)*1.1) +
#   theme_bw() +
#   ylab("Colon-Tumour") +
#   ggtitle("TAPS decon") +
#   scale_fill_manual(values=brewer.pal(6,"Dark2")[c(5,6)])
# 
# pdf("figs/fig7/colontumour_contribution.CRC.pdf",width = 3.5, height = 3)
# print(p5)
# dev.off() 


#### SCAN samples ####
infile_taps <- "all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.all.txt.SCAN_cfTAPS.samples.filtered.CpG.rC.selN200.tissue_contribution.txt"
cfDNA_SCAN <- fread(paste0("cfDNA_deconvolution/",infile_taps)) %>% as.data.frame()
tissue_order <- rev(cfDNA_SCAN$tissue)
pdf("figs/fig7/tissue_contribution.SCAN.pdf",width = 10, height = 6)
cfDNA_SCAN %>%
  melt(id.vars = c("tissue")) %>%
  mutate(tissue = factor(tissue, tissue_order)) %>%
  ggplot(aes(variable, tissue, fill = value*100)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Blues", type = "div", direction = 1,            # <-- RdBu (blue low, red high)
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
cfDNA_SCAN %>%
  melt(id.vars = c("tissue")) %>%
  mutate(tissue = factor(tissue, tissue_order)) %>%
  ggplot(aes(variable, tissue, fill = value*100)) +
  geom_tile() +
  geom_text(aes(label = round(value*100,1)), size = 3,color="black") +
  scale_fill_distiller(
    palette = "Blues", type = "div", direction = 1,            # <-- RdBu (blue low, red high)
    oob = scales::squish,
    na.value = "black"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



dev.off() 


infile_taps <- "Atlas.U250.l4.hg38.full.tsv.SCAN_cfTAPS.samples.filtered.CpG.rC.tissue_contribution.txt"
cfDNA_SCAN <- fread(paste0("cfDNA_deconvolution/",infile_taps)) %>% as.data.frame()
pdf("figs/fig7/tissue_contribution.SCAN_grail.Atlas.U250.pdf",width = 10, height = 6)
cfDNA_SCAN %>%
  melt(id.vars = c("tissue")) %>%
  ggplot(aes(variable, tissue, fill = value*100)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Blues", type = "div", direction = 1,            # <-- RdBu (blue low, red high)
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
cfDNA_SCAN %>%
  melt(id.vars = c("tissue")) %>%
  ggplot(aes(variable, tissue, fill = value*100)) +
  geom_tile() +
  geom_text(aes(label = round(value*100,1)), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "Blues", type = "div", direction = 1,            # <-- RdBu (blue low, red high)
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off() 

infile_taps <- "Atlas.U25.l4.hg38.full.tsv.SCAN_cfTAPS.samples.filtered.CpG.rC.tissue_contribution.txt"
cfDNA_SCAN <- fread(paste0("cfDNA_deconvolution/",infile_taps)) %>% as.data.frame()
pdf("figs/fig7/tissue_contribution.SCAN_grail.Atlas.U25.pdf",width = 10, height = 6)
cfDNA_SCAN %>%
  melt(id.vars = c("tissue")) %>%
  ggplot(aes(variable, tissue, fill = value*100)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Blues", type = "div", direction = 1,            # <-- RdBu (blue low, red high)
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
cfDNA_SCAN %>%
  melt(id.vars = c("tissue")) %>%
  ggplot(aes(variable, tissue, fill = value*100)) +
  geom_tile() +
  geom_text(aes(label = round(value*100,1)), size = 3,color="grey") +
  scale_fill_distiller(
    palette = "Blues", type = "div", direction = 1,            # <-- RdBu (blue low, red high)
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off() 
#### 2021 cfTAPS Sci Adv. ####





