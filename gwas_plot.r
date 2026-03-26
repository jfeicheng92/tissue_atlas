setwd("/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls")

dat <- fread("gwas/gwas_catalog_v1.0_vs_dmr.fisher.txt")
tumour_order <- c("breast","colon","kidney","lung","ovary","pancreas","prostate","stomach")
plot_dat <- dat %>%
  # split "dmr" into tumour_type and feature
  separate(dmr, sep = "\\.", into = c("tumour_type", "feature"), remove = FALSE, extra="merge") %>%
  # try to get p-value from several possible column names; result is "pval"
  mutate(
    pval = coalesce(
      # backtick name safe lookup where hyphen exists
      `two-tail` = if ("two-tail" %in% names(.)) .[["two-tail"]] else NULL,
    )
  ) %>%
  # safety: ensure ratio is numeric
  mutate(ratio = as.numeric(ratio)) %>%
  # significance label
  mutate(
    sig = case_when(
      is.na(pval)           ~ "ns",
      pval < 0.001          ~ "***",
      pval < 0.01           ~ "**",
      pval < 0.05           ~ "*",
      TRUE                  ~ ""
    ),
    p_label = if_else(is.na(pval), "p=NA", 
                      if_else(pval < 0.001, "p<0.001", paste0("p=", signif(pval, 3))))
  ) %>%
  filter(tumour_type %in% tumour_order) %>%
  mutate(tumour_type=factor(tumour_type, levels=rev(tumour_order)))

# ---- compute y positions for annotations ----

pos_dat <- plot_dat %>%
  group_by(tumour_type) %>%
  summarize(group_max = max(ratio, na.rm = TRUE)) %>%
  ungroup()

plot_dat <- plot_dat %>%
  left_join(pos_dat, by = "tumour_type") %>%
  mutate(
    # offset depends on overall range of ratio in that tumour group
    y_pos = group_max + 0.05 * (group_max + 1)  # tweak multiplier if needed
  )

# ---- plotting ----
p <- ggplot(plot_dat, aes(x = tumour_type, y = ratio, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "longdash", color = "gray40") +
  # significance stars above bars (dodged)
  geom_text(
    aes(label = sig, y = y_pos, group = feature),
    position = position_dodge(width = 0.9),
    vjust = 0.5,
    size = 5,
    color = "black"
  ) +
  scale_fill_manual(values=c("#c03b2f","#373c95")) +
  labs(x = "Tumour type", y = "Ratio", fill = "Feature") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  coord_flip() +
  labs(
    x = "Tumour type",
    y = "Ratio",
    fill = "Feature",
    title = "Enrichment of tumour DMRs in GWAS-defined cancer risk loci"
  )

ggsave("figs/fig6/gwas_catalog_v1.0_vs_dmr.fisher.pdf",p, height = 5, width = 5 )
