# A Tri-level methylation atlas of normal and tumour tissues
Authors: Jingfei Cheng<sup>1,2,15</sup> , Masato Inoue<sup>1,2,15</sup> , Jinfeng Chen<sup>3,4,5,15</sup> , Ella Mi<sup>1,2,15</sup> , Felix Jackson<sup>1,2,6,15</sup> , Haiqi Xu<sup>1,2</sup> , Beibei Wang<sup>7</sup> , Yanchun Peng<sup>7,8</sup> , Rory Peters<sup>9</sup> , Sakineh Hussainy<sup>1,2</sup> , Natalie J. Jooss<sup>10</sup> , Bob Amess<sup>1</sup> , Yibin Liu<sup>11,12</sup> , Benjamin Schuster-Böckler<sup>1</sup> , Bethan Psaila<sup>10</sup> , Shivan Sivakumar<sup>13</sup> , Eleanor Barnes<sup>9</sup> , Brian D Nicholson<sup>14</sup> , Tao Dong<sup>7,8</sup> , Chun-Xiao Song<sup>1,2,†</sup> 

Affiliations:  
<sup>1</sup> Ludwig Institute for Cancer Research, Nufﬁeld Department of Medicine, University of Oxford, Oxford, UK  
<sup>2</sup> Target Discovery Institute, Nufﬁeld Department of Medicine, University of Oxford, Oxford, UK 
<sup>3</sup> China National Center for Bioinformation, Beijing 100101, China 
<sup>4</sup> Beijing Institute of Genomics, Chinese Academy of Sciences, Beijing 100101, China 
<sup>5</sup> University of Chinese Academy of Sciences, Beijing 100049, China 
<sup>6</sup> Department of Computer Science, University of Oxford, Oxford, UK 
<sup>7</sup> Chinese Academy of Medical Sciences (CAMS) Oxford Institute (COI), University of Oxford, Oxford, UK  
<sup>8</sup> MRC Translational Immune Discovery Unit, MRC Weatherall Institute of Molecular Medicine, University of Oxford, Oxford, UK 
<sup>9</sup> Oxford University Hospital, Oxford, United Kingdom 
<sup>10</sup> MRC Weatherall Institute of Molecular Medicine, Radcliffe Department of Medicine and National Institute of Health Research, Oxford Biomedical Research Centre, University of Oxford, Oxford, UK 
<sup>11</sup> College of Chemistry and Molecular Sciences, Wuhan University, Wuhan, China 
<sup>12</sup> Taikang Center for Life and Medical Sciences, Wuhan University, Wuhan, China 
<sup>13</sup> Department of Immunology and Immunotherapy, School of Infection, Inflammation and Immunology, College of Medicine and Health, University of Birmingham, Birmingham B15 2TT, UK. 
<sup>14</sup>Nuffield Department of Primary Care Health Sciences, University of Oxford, Oxford OX2 6GG, UK 
<sup>15</sup>These authors contributed equally to this work. 
<sup>†</sup>Corresponding author. E-mail: chunxiao.song@ludwig.ox.ac.uk  
 
## data pre-processing
https://bitbucket.org/bsblabludwig/nxf_workflows/src/master/ramess/CAPS_tissue_map/
https://bitbucket.org/bsblabludwig/nxf_workflows/src/master/ramess/TAPSbeta_tissue_map/

Steps:
* Trim reads with Trim Galore
* Align reads with bwa-mem2
* Mark duplicate reads with Picard MarkDuplicates
* Call methylation with MethylDackel
* MLML was used to integrated 5mC and 5hmC to obtain tri-level cytosine estimates (5umC, 5mC, 5hmC)

## DMB calling
dmr_call.py
```bash
methratio=all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.bed
for c in mC umC hmC
do
python dmr_call.py --in_file $methratio  --tg_quant 0.25 --bg_quant_hypo 0.05 --bg_quant_hyper 0.05 --bg_quant_mode groups --top_n 2000 --libs $c >dmr_call_bg_quant0.05.${c}.log 2>&1 
done
```

## Code for figures
**Fig. 1.** Tri-level DNA Methylation Atlas.  
A-B. Schematic plot (BioRender)  
C. Whole-genome average levels of 5mC, 5hmC, and 5umC across all samples 
`ternary_plot.r`  
D. Differences in whole-genome cytosine modification levels between tumour - normal pairs 
`ternary_plot.r`  

**Fig. 2.** Genomic distribution of 5uC, 5hmC and 5mC.  
A. Feature enrichment analysis of genomic regions with maximal modification levels (top 0.5%)  
`high_mod.r`  
B. Average tri-level methylation profiles (5umC, 5hmC, 5mC) in 1kb bins at example tissue-specific genes  
`stacked_bar_igv.r`  
C. Hierarchical clustering of normal samples based on the top 0.5% most variable 
`top_var.r`  

**Fig. 3.** Tissue specific DMRs marked spatially distinct regions.  
A. Pairwise correlations between 5umC, 5hmC, and 5mC  
B. Proportion of hyper- and hypo-DMRs  
C. Heatmap of the top 200 tissue-specific DMRs  
D. Venn diagrams showing overlap between the top 1,000 tissue-specific DMRs  
E. Boxplots showing the percentage of top 1,000 tissue-specific DMRs in another context  
`tissue_specific_dmr_1kb.r`  


**Fig. 4.**  5(u)mC and 5hmC marked tissue specific enhancers.  
A-B. ChIP-seq signals aligned with tissue specific DMRs  
`tissue_specific_dmr_1kb_histone.r`  
C. ATAC-seq within hypo-5mC and hyper-5hmC DMRs  
`tissue_specific_dmr_1kb_atac-seq.r`  
D. Enriched sequence motifs
`tissue_specific_dmr_1kb_motif.r`  
E. Schematic plot (BioRender)  

**Fig. 5.** Relationship between tissue-specific DMRs and gene expression.  
A. Distribution of DNA modifications across gene bodies and ±10 kb flanking regions  
`meth_vs_gene_expression_bins_quantiles.r`  
B. Enrichment of genes proximal to tissue-specific DMRs among tissue-specific gene sets  
`tissue_specific_dmr_1kb_gene_expr.r`  
C. Heatmap showing scaled average expression of genes associated with tissue-specific DMRs
`tissue_specific_dmr_1kb_gene_expr.r`  
D. Spatial distribution of DMRs  
`tumour_specific_dmr_1kb_gene_expr_test.r`
E. Pearson correlation between predicted and observed gene expression levels   
`Gene_expression_prediction.ipynb`
`gene_expr_predict.r`

**Fig. 6.** Integrated analysis of tumour type-specific DMRs. 
A. Overall methylation difference between tumour - normal pairs  
`tumour_meth_1kb.r`  
B. Heatmap of the top 200 tumour-specific DMRs  
`tumour_specific_dmr_1kb.r`
C. Enrichment analysis of tumour-specific DMRs against tumour type-specific genes
`tumour_gene.r`
`tumour_specific_dmr_1kb_gene_expr_test.r`
D. ATAC-seq within hypo-5mC and hyper-5hmC DMRs
`tissue_specific_dmr_1kb_atac-seq.r` 
E. Motif enrichment analysis of tumour-specific DMRs
`tissue_specific_dmr_1kb_motif.r`
F. Enrichment analysis of tumour DMRs at cancer GWAS loci
`gwas.sh`
`gwas_plot.r`

**Fig. 7.** Comparison of tissue deconvolution using TAPS and CAPS+ Atlases
`tissue_deconvolution_CAPS.ipynb`
`tissue_deconvolution_TAPS.ipynb`
`tissue_deconvolution.r`


