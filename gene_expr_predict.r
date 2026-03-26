res_hmC <- fread("/well/ludwig/users/wdu564/gene_expression/hmC_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("hmC_",c("r2","mse","pearsonr"))))
res_mC <- fread("/well/ludwig/users/wdu564/gene_expression/mC_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("mC_",c("r2","mse","pearsonr"))))
res_uC <- fread("/well/ludwig/users/wdu564/gene_expression/umC_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("uC_",c("r2","mse","pearsonr"))))
res_cmb <- fread("/well/ludwig/users/wdu564/gene_expression/combined_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("cmb_",c("r2","mse","pearsonr"))))

res_all_cnn <- merge(res_hmC, res_mC, by = c("fold","tissues")) %>%
  merge(., res_uC, by = c("fold","tissues")) %>%
  merge(., res_cmb, by = c("fold","tissues")) %>%
  mutate(
    type = "solid",
    type = ifelse(tissues %in% c(
      "CD4-T-cells","CD8-T-cells","NK-cells","B-cells",
      "Neutrophils","Eosinophils","Monocytes",
      "Erythroid-precursors","Megakaryocytes"
    ), "blood", type),
    type = factor(type, levels = c("blood","solid")) # optional order
  )
p1 <- res_all_cnn %>%
  select(tissues, contains("pearsonr"), type) %>%
  pivot_longer(
    cols = -c(type, tissues),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(variable = factor(variable, levels = c("cmb_pearsonr","mC_pearsonr","uC_pearsonr", "hmC_pearsonr"))) %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = round(..y.., 2)),
    vjust = -6,
    size = 3.5
  ) +
  geom_jitter(aes(color = type), width = 0.1) +
  scale_color_manual(values = c(blood = "#6BAED6", solid = "#E34A33")) +
  ylim(0, 1) +
  theme_classic() +
  labs(
    y = "pearsonr",
    x = "feature",
    title = "cnn",
    fill = "type"
  )

# > mean(res_all_cnn$cmb_pearsonr)
# [1] 0.8101778
# > mean(res_all_cnn$hmC_pearsonr)
# [1] 0.6605703
# > mean(res_all_cnn$mC_pearsonr)
# [1] 0.7552283
# > mean(res_all_cnn$uC_pearsonr)
# [1] 0.7401655
# > t.test(res_all_cnn$cmb_pearsonr, res_all_cnn$hmC_pearsonr)
# 
# Welch Two Sample t-test
# 
# data:  res_all_cnn$cmb_pearsonr and res_all_cnn$hmC_pearsonr
# t = 13.215, df = 33.525, p-value = 7.566e-15
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1265879 0.1726270
# sample estimates:
#   mean of x mean of y 
# 0.8101778 0.6605703 
# 
# > t.test(res_all_cnn$cmb_pearsonr, res_all_cnn$mC_pearsonr)
# 
# Welch Two Sample t-test
# 
# data:  res_all_cnn$cmb_pearsonr and res_all_cnn$mC_pearsonr
# t = 4.3021, df = 36, p-value = 0.0001237
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.02904517 0.08085390
# sample estimates:
#   mean of x mean of y 
# 0.8101778 0.7552283 
# 
# > t.test(res_all_cnn$cmb_pearsonr, res_all_cnn$uC_pearsonr)
# 
# Welch Two Sample t-test
# 
# data:  res_all_cnn$cmb_pearsonr and res_all_cnn$uC_pearsonr
# t = 5.3244, df = 35.882, p-value = 5.609e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.04334101 0.09668366
# sample estimates:
#   mean of x mean of y 
# 0.8101778 0.7401655 
res_hmC <- fread("/well/ludwig/users/wdu564/gene_expression/rf_hmC_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("hmC_",c("r2","mse","pearsonr"))))
res_mC <- fread("/well/ludwig/users/wdu564/gene_expression/rf_mC_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("mC_",c("r2","mse","pearsonr"))))
res_uC <- fread("/well/ludwig/users/wdu564/gene_expression/rf_umC_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("uC_",c("r2","mse","pearsonr"))))
res_cmb <- fread("/well/ludwig/users/wdu564/gene_expression/rf_combined_gene_expression_prediction_metrics.txt", col.names = c("fold","tissues", paste0("cmb_",c("r2","mse","pearsonr"))))
res_all_rf <- merge(res_hmC, res_mC, by = c("fold","tissues")) %>%
  merge(., res_uC, by = c("fold","tissues")) %>%
  merge(., res_cmb, by = c("fold","tissues")) %>%
  mutate(
    type = "solid",
    type = ifelse(tissues %in% c(
      "CD4-T-cells","CD8-T-cells","NK-cells","B-cells",
      "Neutrophils","Eosinophils","Monocytes",
      "Erythroid-precursors","Megakaryocytes"
    ), "blood", type),
    type = factor(type, levels = c("blood","solid")) # optional order
  )
p2 <- res_all_rf %>%
  select(tissues, contains("pearsonr"), type) %>%
  pivot_longer(
    cols = -c(type, tissues),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(variable = factor(variable, levels = c("cmb_pearsonr","mC_pearsonr","uC_pearsonr", "hmC_pearsonr"))) %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = round(..y.., 2)),
    vjust = -6,
    size = 3.5
  ) +
  geom_jitter(aes(color = type), width = 0.1) +
  ylim(0, 1) +
  scale_color_manual(values = c(blood = "#6BAED6", solid = "#E34A33")) +
  theme_classic() +
  labs(
    y = "pearsonr",
    x = "feature",
    title = "rf",
    fill = "type"
  )



setwd("/users/ludwig/cfo155/cfo155/tissueMap/methcalls/")
pdf("figs/fig5/gene_expr_prediction_pearsonr.pdf", width = 7, height = 3)
cowplot::plot_grid(p1,p2)
dev.off()
# > mean(res_all_rf$cmb_pearsonr)
# [1] 0.7719984
# > mean(res_all_rf$hmC_pearsonr)
# [1] 0.6875054
# > mean(res_all_rf$mC_pearsonr)
# [1] 0.7327585
# > mean(res_all_rf$uC_pearsonr)
# [1] 0.7231747
# t.test(res_all_rf$cmb_pearsonr, res_all_rf$hmC_pearsonr)
# 
# Welch Two Sample t-test
# 
# data:  res_all_rf$cmb_pearsonr and res_all_rf$hmC_pearsonr
# t = 9.2079, df = 35.244, p-value = 6.571e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.06586904 0.10311693
# sample estimates:
#   mean of x mean of y 
# 0.7719984 0.6875054 
# 
# > t.test(res_all_rf$cmb_pearsonr, res_all_rf$mC_pearsonr)
# 
# Welch Two Sample t-test
# 
# data:  res_all_rf$cmb_pearsonr and res_all_rf$mC_pearsonr
# t = 2.9704, df = 30.009, p-value = 0.005805
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01226124 0.06621846
# sample estimates:
#   mean of x mean of y 
# 0.7719984 0.7327585 
# 
# > t.test(res_all_rf$cmb_pearsonr, res_all_rf$uC_pearsonr)
# 
# Welch Two Sample t-test
# 
# data:  res_all_rf$cmb_pearsonr and res_all_rf$uC_pearsonr
# t = 4.4011, df = 34.402, p-value = 9.92e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.02628882 0.07135860
# sample estimates:
#   mean of x mean of y 
# 0.7719984 0.7231747 
