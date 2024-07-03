# A methylation and hydroxymethylation atlas of normal and tumour tissues 
Authors: Masato Inoue<sup>1,2,9</sup>, Jingfei Cheng<sup>1,2,9</sup>, Felix Jackson<sup>1,2,3,9</sup>, Jinfeng Chen<sup>1,2,9,10</sup>, Haiqi Xu<sup>1,2</sup>, Beibei Wang<sup>4</sup>, Yanchun Peng<sup>4,5</sup>, Natalie J. Jooss<sup>6</sup>, Bob Amess<sup>1</sup>, Yibin Liu<sup>7,8</sup>, Benjamin Schuster-Böckler<sup>1</sup>, Bethan Psaila<sup>6</sup>, Tao Dong<sup>4,5</sup>, Chun-Xiao Song<sup>1,2,†</sup>  
Affiliations:  
<sup>1</sup>Ludwig Institute for Cancer Research, Nufﬁeld Department of Medicine, University of Oxford, Oxford, UK  
<sup>2</sup>Target Discovery Institute, Nufﬁeld Department of Medicine, University of Oxford, Oxford, UK  
<sup>3</sup>Department of Computer Science, University of Oxford, Oxford, UK  
<sup>4</sup>Chinese Academy of Medical Sciences (CAMS) Oxford Institute (COI), University of Oxford, Oxford, UK  
<sup>5</sup>MRC Translational Immune Discovery Unit, MRC Weatherall Institute of Molecular Medicine, University of Oxford, Oxford, UK  
<sup>6</sup>MRC Weatherall Institute of Molecular Medicine, Radcliffe Department of Medicine and National Institute of Health Research, Oxford Biomedical Research Centre, University of Oxford, Oxford, UK  
<sup>7</sup>College of Chemistry and Molecular Sciences, Wuhan University, Wuhan, China  
<sup>8</sup>Taikang Center for Life and Medical Sciences, Wuhan University, Wuhan, China  
<sup>9</sup>These authors contributed equally to this work  
<sup>10</sup>Present address: CAS Key Laboratory of Genome Sciences and Information, Beijing Institute of Genomics, Chinese Academy of Sciences and China National Center for Bioinformation, Beijing 100101, China  
<sup>†</sup>Corresponding author. E-mail: chunxiao.song@ludwig.ox.ac.uk
 
## data preprocessing
https://bitbucket.org/bsblabludwig/nxf_workflows/src/master/ramess/CAPS_tissue_map/
https://bitbucket.org/bsblabludwig/nxf_workflows/src/master/ramess/TAPSbeta_tissue_map/

Steps:
* Trim reads with Trim Galore
* Align reads with bwa-mem2
* Mark duplicate reads with Picard MarkDuplicates
* Call methylation with MethylDackel

## Genome segmentation

## DMB calling

## DMB analysis
dmr15 - change back to delta_means=0.01 for all, including cirrhosis + pancreatitis

## Gene expression prediction

## Code for figures
Fig. 1.

Fig. 2. Differentially modified blocks of 5mC and 5hmC are associated with tissue-specific gene expression.  
A. Heatmap  
`dmb_downstream_analysis.r`  
B. DMB number  
`dmb_downstream_analysis.r`  
C. DMB vs. gene enrich   
`dmb_downstream_analysis.r`  

Fig. 3. D(h)MBs of normal tissues and blood cells mark regulatory regions.  
A-B. Histone  
`heatlhy_dmb_validation.sh`  
`blocks_analysis_histone_deeptools.r`  
C. ATAC-seq  
`atac_seq.sh`  
`dmb_downstream_analysis.r`  
D. Motif
``
``

Fig. 4. Tumours gain unique signatures of 5mC and 5hmC.
A. Heatmap 
`dmb_downstream_analysis.r`  
B. DMB vs. gene enrich   
`dmb_downstream_analysis.r`  
C. ATAC-seq  
`atac_seq.sh`  
`dmb_downstream_analysis.r` 

Fig. 5. 5mC and 5hmC levels predict gene expression. 


Fig S1  
Fig S2  

Fig S5
Heatmap  
`dmb_downstream_analysis.r`  

Fig S14.
Heatmap 
`dmb_downstream_analysis.r`  

Fig S16.
Heatmap  
`dmb_downstream_analysis.r`  
