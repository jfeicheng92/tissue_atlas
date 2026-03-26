library(data.table)
library(ggplot2)
library(dplyr)

setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/")
meth <- fread("meth/all_sample.merged.mlml.mincov10_common.summary.txt")
meth_wide <- meth %>%
  tidyr::separate(
    col = V1,              
    into = c("ID", "tissue", "mark"),
    sep = "_"
  ) %>%
  tidyr::pivot_wider(
    names_from = mark,     # column to become new column names
    values_from = V2       # column with the values
  )
colnames(meth_wide) <- c("ID","tissue","ncpg","5mC","5hmC","5umC")
dat <- meth_wide
dat$tissue <- gsub("CD34-erythroblasts","Erythroid-precursors",dat$tissue)
dat$tissue <- gsub("CD34-megakaryocytes","Megakaryocytes",dat$tissue)
dat$tumour <-"no"; dat$tumour[grep("Tumour",dat$tissue)] <- "yes";dat$tumour[grep("-Pancreatitis|-Cirrhosis",dat$tissue)] <- "pre"
dat$tissue <- gsub("-Tumour|-Pancreatitis|-Cirrhosis","",dat$tissue)
tissue_order <- c( "Brain","Breast", "Heart","Kidney","Liver", "Lung","Ovary", "Pancreas", "Prostate","Colon","Stomach","Esophagus","Spleen", 
                   "CD4-T-cells", "CD8-T-cells", "NK-cells", "B-cells", "Neutrophils", "Eosinophils", "Monocytes", "Erythroid-precursors", "Megakaryocytes")
dat$tissue <- factor(dat$tissue,levels=tissue_order)
dat <- dat[order(dat$tissue),]
dat_long <- dat %>%
  tidyr::pivot_longer(
    cols = c("5mC","5hmC","5umC"),       # columns to reshape
    names_to = "type",        # new column with the old column names
    values_to = "value",          # new column with the values
    values_drop_na = TRUE         # optional: drop NAs
  )
dat_long$type <- factor(dat_long$type, levels=c("5mC","5hmC","5umC"))

pdf("figs/fig1/all_sample.merged.mlml.mincov10_common.summary.all.pdf", width = 10, height = 6)
ggplot(dat_long,aes(x=tissue,y=value, group=tumour, color=tumour)) + 
  geom_point(aes(shape=tumour), position=position_dodge(width=0.85)) +
  facet_grid(type~.,scale="free_y") +
  scale_shape_manual(values=c(16, 15, 8)) +
  scale_color_manual(values=c("#473183","#829fd5","#e8908f"))+
  theme_bw()+
  theme(legend.position = NULL)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color="black")) +
  ylab("genome methylation %") + xlab("")+
  stat_summary(fun = "mean", 
               geom = "errorbar", color="black",
               aes(ymax = ..y.., ymin = ..y..), 
               position = position_dodge(width = 0.85), 
               width = 0.5)
ggplot(dat_long,aes(x=tissue,y=value, group=tumour, color=tumour)) + 
  geom_point(aes(shape=tumour), position=position_dodge(width=0.85)) +
  facet_grid(type~.,scale="free_y") +
  scale_shape_manual(values=c(16, 15, 8)) +
  scale_color_manual(values=c("#473183","#829fd5","#e8908f"))+
  theme_bw()+
  theme(legend.position = NULL)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color="black")) +
  ylab("genome methylation %") + xlab("")+
  stat_summary(fun = "mean", colour = "black", size = 1.8,
               geom = "text", aes(label = round(after_stat(y), 1)),
               position = position_dodge(width=0.85), vjust = 1.5) +
  stat_summary(fun = "mean", 
               geom = "errorbar", color="black",
               aes(ymax = ..y.., ymin = ..y..), 
               position = position_dodge(width = 0.85), 
               width = 0.5)


dev.off()

# lifespan <- read.table("resource/lifespan_tissue.txt",sep="\t",header = TRUE)
# select_mean <- dat_long %>%
#   group_by(tissue, tumour, type) %>%
#   summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
#   filter(type=="5hmC") %>%
#   filter(tumour=="no")%>%
#   as.data.frame() %>%
#   merge(.,lifespan,by.x=c("tissue"),by.y="tissue_atlas")
# 
# ggplot(select_mean,aes(x=mean,y=mean_lifespan_in_human_d)) +geom_point()


# tissues to include
sel_tissues <- c("Brain","Breast","Colon","Kidney","Liver","Lung","Ovary","Pancreas","Prostate","Stomach")

# the methylation marks to summarise
marks <- c("5mC","5hmC","5umC")

# helper to safely run t-tests only when variance/length permit
.safe_ttest_p <- function(x, y){
  if(length(x) > 1 && length(y) > 1 &&
     sd(x, na.rm = TRUE) > 0 && sd(y, na.rm = TRUE) > 0) {
    t.test(x, y)$p.value
  } else NA_real_
}

# build summary table: mean_normal, mean_tumour, delta, p-values for each tissue × mark
summary_df <- purrr::map_dfr(marks, function(m) {
  
  dat_m <- dat %>%
    filter(tissue %in% sel_tissues) %>%
    transmute(
      tissue,
      tumour,
      value = .data[[m]]   # avoids backticks for names like "5mC"
    ) %>%
    filter(!is.na(value))
  
  dat_m %>%
    group_split(tissue) %>%
    purrr::map_dfr(function(df) {
      n_vals <- df %>% filter(tumour == "no")  %>% pull(value)
      t_vals <- df %>% filter(tumour == "yes") %>% pull(value)
      
      mean_normal <- mean(n_vals, na.rm = TRUE)
      mean_tumour <- mean(t_vals, na.rm = TRUE)
      
      tibble(
        tissue      = first(df$tissue),
        mark        = m,
        n_normal    = length(n_vals),
        n_tumour    = length(t_vals),
        mean_normal = mean_normal,
        mean_tumour = mean_tumour,
        delta       = mean_tumour - mean_normal,
        p_wilcox    = if(length(n_vals) > 0 && length(t_vals) > 0)
          wilcox.test(t_vals, n_vals, exact = FALSE)$p.value else NA_real_,
        p_ttest     = .safe_ttest_p(t_vals, n_vals)
      )
    })
})


delta_wide <- summary_df %>%
  select(tissue, mark, delta) %>%
  pivot_wider(names_from = mark, values_from = delta)

delta_wide
summary_df <- summary_df %>%
  mutate(sig = case_when(
    p_ttest < 0.001 ~ "***",
    p_ttest < 0.01  ~ "**",
    p_ttest < 0.05  ~ "*",
    TRUE ~ ""
  ))
summary_df$mark <- factor(summary_df$mark, levels=c("5mC","5hmC","5umC"))
pdf("figs/fig1/tumour_delta_ttest.pdf",width = 5,height = 8)
ggplot(summary_df, aes(x = mark, y = tissue, fill = delta)) +
  geom_tile() +
  geom_text(aes(label = sig), color = "black", size = 5) +   # adds "*" symbols
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "tumour - normal") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
write.csv(summary_df, "figs/fig1/tumour_delta_ttest.csv", row.names = FALSE)


##### bash command #####
# WORKDIR=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/meth/
# zcat all_sample.merged.mlml.mincov10_common.gz |\
#     awk 'NR==1 {            # first line (header)
#     for (i=4; i<=NF; i++) header[i]=$i
#     next
# }
# {
#     for (i=4; i<=NF; i++) {
#         sum[i]+=$i
#         count[i]++
#     }
# }
# END {
#     for (i=4; i<=NF; i++) {
#         avg = (count[i] ? sum[i]/count[i] : 0)
#         printf "%s\t%.2f\t%d\n", header[i], avg, count[i]
#     }
# }' >all_sample.merged.mlml.mincov10_common.summary.txt



#### diagram show tri-level methylation ####

df <- data.frame(
  category = c("mC", "hmC", "umC"),
  value = c(68.5, 7.8, 100 - 68.5 - 7.8)
)

df$pos <- cumsum(df$value) - df$value/2

ggplot(df, aes(x = "", y = value, fill = category)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(y = pos, label = paste0(category, " (", round(value, 1), "%)"))) +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B", "#76B7B2")) +
  labs(title = "Composition (%)") +
  theme_void()
