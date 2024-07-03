#!/bin/bash
#SBATCH --job-name=compute_matrix-array
#SBATCH -p long
#SBATCH --array 0-11:1
#SBATCH --requeue 


# Define array of dmb and histone
dmb=(TAPSbeta_Hypo CAPS_Hyper)
his=(H3K4me1 H3K27ac H3K4me3 H3K36me3 H3K27me3 H3K9me3)
dmb_idx=$((SLURM_ARRAY_TASK_ID % ${#dmb[@]}))
his_idx=$((SLURM_ARRAY_TASK_ID / ${#dmb[@]}))
sel_dmb="${dmb[$dmb_idx]}"
sel_his="${his[$his_idx]}"

# H3K4me1 72 # H3K4me1 is predominantly enriched at enhancers
# H3K27ac 73 # H3K27ac localizes to both active gene promoters and enhancer regions
# H3K4me3 148 # H3K4me3 and H3K9ac are associated with transcriptionally active gene promoter regions
# H3K9ac 5
# H3K36me3 122 # H3K36me3 are localized to gene bodies of actively transcribing genes
# H3K27me3 127
# H3K9me3 130 # silent genes have high levels of H3K9me2/3 and H3K27me3 widespread across the gene region





# tapsb_dmb=TAPSbeta.healthy_high_cov_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
# caps_dmb=CAPS.healthy_high_cov_samples.Tissue_Group_DMBs.top300.CAPS_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
# gsize=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.fai
# cat $tapsb_dmb|awk 'BEGIN{OFS="\t"}{print $1,($(NF-1)),($NF)}'|awk '$2=="Hypo"' |sed 's/_/\t/g' |bedtools sort -i - -faidx <(grep -w chr[0-9,X,Y,M]* $gsize)|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4"_"$5}' >TAPSbeta_Hypo_DMB.bed
# cat $caps_dmb |awk 'BEGIN{OFS="\t"}{print $1,($(NF-1)),($NF)}'|awk '$2=="Hyper"'|sed 's/_/\t/g' |bedtools sort -i - -faidx <(grep -w chr[0-9,X,Y,M]* $gsize)|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4"_"$5}' >CAPS_Hyper_DMB.bed

#### Histone ####
module purge
module load deepTools/3.3.1-foss-2018b-Python-3.6.6
module load plotly.py/4.4.1-foss-2018b-Python-3.6.6


# computeMatrix reference-point --referencePoint center -S `ls /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/ihec/*.${sel_his}.signal_unstranded.bigWig` \
#                 -R ${sel_dmb}_DMB.bed \
#                 -a 15000 -b 15000 -bs 100 --missingDataAsZero \
#                 -o ${sel_dmb}_DMB.${sel_his}.gz -p 4

computeMatrix reference-point --referencePoint center -S `ls /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/histone/rename_tracks/*.${sel_his}*bigWig` \
                -R ${sel_dmb}_DMB.bed \
                -a 15000 -b 15000 -bs 100 --missingDataAsZero \
                -o ${sel_dmb}_DMB.histone_${sel_his}.gz -p 4

computeMatrix reference-point --referencePoint center -S `ls /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/histone/rename_tracks/*.${sel_his}*bigWig` \
                -R ${sel_dmb}_DMB.bed \
                -a 15000 -b 15000 -bs 100 --missingDataAsZero \
                -o ${sel_dmb}_DMB.histone_${sel_his}.remove_blacklist.gz -p 4 --blackListFileName /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/hg38_blacklist_centromer.bed

#   # #### Motif ####
#   # wget -c http://homer.ucsd.edu/homer/configureHomer.pl
#   # perl configureHomer.pl -install homer
#   # perl /users/ludwig/cfo155/cfo155/tissueMap/methcalls/homer/configureHomer.pl -install hg38
#   # PATH=$PATH:/users/ludwig/cfo155/cfo155/tissueMap/methcalls/homer/bin/
#   # tapsb_dmb=TAPSbeta.healthy_high_cov_samples.Tissue_Group_DMBs.top1000.TAPSbeta_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   # caps_dmb=CAPS.healthy_high_cov_samples.Tissue_Group_DMBs.top1000.CAPS_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   # 
#   # mkdir dmbs
#   # cat $tapsb_dmb|awk 'BEGIN{OFS="\t"}{print $1,($(NF-1)),($NF)}'|sed 's/_/\t/g' |sort -k1,1 -k2,2n|grep -v start |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR,".","+" >"dmbs/TAPSbeta_"$4"_"$5".bed"}'
#   # cat $caps_dmb |awk 'BEGIN{OFS="\t"}{print $1,($(NF-1)),($NF)}'|sed 's/_/\t/g' |sort -k1,1 -k2,2n|grep -v start |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR,".","+" >"dmbs/CAPS_"$4"_"$5".bed"}'
#   # for i in `ls dmbs/CAPS_Hyper*bed dmbs/TAPSbeta_Hypo*bed`
#   # do
#   # findMotifsGenome.pl  <(cut -f1-3,6 $i) hg38 ${i/.bed/} -size 250 -len 8 
#   # done
#   
#   for i in `ls */knownResults.txt`
#   do
#       awk -v i=`dirname $i` 'BEGIN{FS="\t";OFS="\t"}{if($5<0.01)print i,$0}' $i
#   done |cat <(echo -e "bed\t`head -1 CAPS_Hyper_B-cells/knownResults.txt`") - >all_motif_q0.01.txt
#   
#   for i in `ls */knownResults.txt`
#   do
#       awk -v i=`dirname $i` 'BEGIN{FS="\t";OFS="\t"}{if(NR>1)print i,$0}' $i
#   done |cat <(echo -e "bed\t`head -1 CAPS_Hyper_B-cells/knownResults.txt`") - >all_motif.txt
#   
#   https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.24.tgz
#   https://github.com/vierstralab/motif-clustering/tree/master
#   
#   
#   
#   dir=/users/ludwig/cfo155/cfo155/tissueMap/methcalls/homer/data/knownTFs
#   meme=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/meme/meme-5.5.5/meme
#   $meme/bin/tomtom \
#   	-dist kullback \
#   	-motif-pseudo 0.1 \
#   	-text \
#   	-min-overlap 1 \
#   	all_motif.meme all_motif.meme > all_motif_tomtom.all.txt
#   $meme/libexec/meme-5.5.5/taipale2meme all_motif.taiple >all_motif.meme
#   cat *fmt >all_motif.taiple
#   motif_clustering.ipynb
#   
#   output all_motif_tomtom.all.cluster.txt
#   for i in `ls *motif`;do echo -e "$i\t`head -1 $i|cut -f2`";done|\
#       sed 's/.motif//g'|sort -k1,1 |\
#       join -1 1 - -2 1 <(sort -k1,1 all_motif_tomtom.all.cluster.txt ) -t$'\t'|\
#       sort -k3,3n >all_motif_tomtom.all.cluster.name.txt
#   
#   
#   #### Feature enrichment ####
#   
#   features=(CGIshelves.bed
#   CGIshore.bed
#   cpgIsland.bed
#   GRCh38.Regulatory_Build.CTCF_binding_site.bed
#   GRCh38.Regulatory_Build.TF_binding_site.bed
#   GRCh38.Regulatory_Build.enhancer.bed
#   GRCh38.Regulatory_Build.open_chromatin_region.bed
#   GRCh38.Regulatory_Build.promoter.bed
#   GRCh38.Regulatory_Build.promoter_flanking_region.bed
#   MANE.GRCh38.v1.0.refseq_genomic.gene.bed
#   MANE.GRCh38.v1.0.refseq_genomic.gene.down1k.bed
#   MANE.GRCh38.v1.0.refseq_genomic.gene.up1k.bed
#   MANE.GRCh38.v1.0.refseq_genomic.exon.bed
#   MANE.GRCh38.v1.0.refseq_genomic.intron.bed
#   hg38.repeatmasker.DNA_TcMar.bed
#   hg38.repeatmasker.DNA_hAT.bed
#   hg38.repeatmasker.LINE_L1.bed
#   hg38.repeatmasker.LINE_L2.bed
#   hg38.repeatmasker.LTR_ERV.bed
#   hg38.repeatmasker.SINE_Alu.bed
#   hg38.repeatmasker.SINE_MIR.bed
#   hg38.repeatmasker.Satellite_centromeric.bed
#   hg38lift_genome_100_segments/Acet.bed
#   hg38lift_genome_100_segments/BivProm.bed
#   hg38lift_genome_100_segments/DNase.bed
#   hg38lift_genome_100_segments/EnhA.bed
#   hg38lift_genome_100_segments/EnhWk.bed
#   hg38lift_genome_100_segments/GapArtf.bed
#   hg38lift_genome_100_segments/HET.bed
#   hg38lift_genome_100_segments/PromF.bed
#   hg38lift_genome_100_segments/Quies.bed
#   hg38lift_genome_100_segments/ReprPC.bed
#   hg38lift_genome_100_segments/TSS.bed
#   hg38lift_genome_100_segments/Tx.bed
#   hg38lift_genome_100_segments/TxEnh.bed
#   hg38lift_genome_100_segments/TxEx.bed
#   hg38lift_genome_100_segments/TxWk.bed
#   hg38lift_genome_100_segments/znf.bed
#   PMD_coordinates_hg38.commonPMD.bed
#   )
#   
#   for i in ${features[@]}
#   do
#       awk -v n=${i/.bed/} 'BEGIN{OFS="\t"}{print $1,$2,$3,n}' $i
#   done|sort -k1,1 -k2,2n >all_feature.bed
#   # normal 
#   #/users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr9
#   tapsb_dmb=TAPSbeta.healthy_high_cov_samples.Tissue_Group_DMBs.top1000.TAPSbeta_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   caps_dmb=CAPS.healthy_high_cov_samples.Tissue_Group_DMBs.top1000.CAPS_hmm.healthy_high_cov_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   cut -f1,88- $tapsb_dmb|grep Hypo |cut -f1|sed 's/_/\t/g' |grep -v start >${tapsb_dmb/.csv/}.hypo.bed
#   cut -f1,88- $caps_dmb|grep Hyper |cut -f1|sed 's/_/\t/g' |grep -v start >${caps_dmb/.csv/}.hyper.bed
#   # /users/ludwig/cfo155/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr10
#   tapsb_dmb=TAPSbeta.all_samples.Tissue_Group_DMBs.top300.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   caps_dmb=CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   cut -f1,122- $tapsb_dmb|grep Tumour|grep Hypo |cut -f1|sed 's/_/\t/g' |grep -v start >${tapsb_dmb/.csv/}.hypo.bed
#   cut -f1,122- $caps_dmb|grep Tumour|grep Hyper |cut -f1|sed 's/_/\t/g' |grep -v start >${caps_dmb/.csv/}.hyper.bed
#   
#   bed=(${tapsb_dmb/.csv/}.hypo.bed ${caps_dmb/.csv/}.hyper.bed)
#   # bed=(TAPSbeta_shared_dmrs.tumour_vs_healthy.Hyper.bed TAPSbeta_shared_dmrs.tumour_vs_healthy.mixed.bed TAPSbeta_shared_dmrs.tumour_vs_healthy.Hypo.bed)
#   faidx=/gpfs3/well/ludwig/users/cfo155/tissueMap/others/mC/resource/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.fai
#   # #### shuffle ####
#   # shuffle using genome as background
#   for j in ${bed[@]}
#   do
#       exp=`echo $j |cut -d '_' -f1`
#       temp=${j/.bed/}.shuffle.genome.bed
#       for seed in `seq 1 1 100`
#       do
#       bedtools shuffle -i $j -g <(grep -w chr[0-9]* $faidx) -seed $seed |\
#           awk -v s=$seed 'BEGIN{OFS="\t"}{print $0,s}'
#       done >$temp
#       awk -v s=0 'BEGIN{OFS="\t"}{print $0,s}' $j >>$temp
#   done
#   
#   featuredir=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/
#   # cat hg38lift_genome_100_segments.bed |awk 'BEGIN{OFS="\t"}{gsub(".*_|[0-9]*","",$4);print $1,$2,$3 >"hg38lift_genome_100_segments/"$4".bed"}'
#   for j in ${bed[@]}
#   do
#       temp=${j/.bed/}.shuffle.genome.bed
#       for i in ${features[@]}
#       do
#        echo -e "$i\t$j\t\
#       `bedtools intersect -a <(cut -f1-3 $featuredir/$i|grep -w chr[0-9]*|bedtools sort -i - -g $faidx) -b <(cat $temp|grep -w chr[0-9]*|bedtools sort -i - -g $faidx) -wa -wb -g <(grep -w chr[0-9]* $faidx ) -sorted |
#       cut -f4-|sort -u|cut -f4|sort |uniq -c |sed 's/^ \+//g;s/ /\t/g'|sort -k2,2|join -1 1 <(seq 0 1 100|sort -k1,1) -2 2 - -t$'\t' -a1|awk 'BEGIN{OFS="\t"}{if($2=="")print $1,0;else print $1,$2}' |cut -f2|tr '\n' '\t'`"
#       done
#   done >feature_enrichment.shuffle.genome.txt
#   
#   
#   ln -s /well/ludwig/users/qwb368/taps_tissue_atlas/dmb/TAPSbeta.all_samples.Tissue_Group_DMBs.all_delta_quants0.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.1.csv
#   ln -s /well/ludwig/users/qwb368/taps_tissue_atlas/dmb/CAPS.all_samples.Tissue_Group_DMBs.all_delta_quants0.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.csv
#   ln -s /well/ludwig/users/dyp502/tissue_atlas_v3/dmr13/CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   
#   csv=CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.csv
#   cut -f1,127- $csv|grep Hyper |grep Tumour|cut -f1|sed 's/_/\t/g' |grep -v start |sed 's/"//g'>${csv/.csv/}.hyper.bed
#   cut -f1,127- $csv|grep Hypo |grep Tumour|cut -f1|sed 's/_/\t/g' |grep -v start |sed 's/"//g'>${csv/.csv/}.hypo.bed
#   cut -f1,127- $csv|grep Tumour|sed 's/_/\t/g' |grep -v start |sed 's/"//g'|cut -f1-5 |awk -v n=${csv/.csv/} 'BEGIN{OFS="\t"}{print $1,$2,$3>n"."$4"."$5".bed"}'
#   csv=CAPS.all_samples.Tissue_Group_DMBs.all_delta_quants0.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.csv
#   cut -f1,127- $csv|grep Hyper |grep Tumour|cut -f1|sed 's/_/\t/g' |grep -v start |sed 's/"//g'>${csv/.csv/}.hyper.bed
#   cut -f1,127- $csv|grep Hypo |grep Tumour|cut -f1|sed 's/_/\t/g' |grep -v start |sed 's/"//g'>${csv/.csv/}.hypo.bed
#   cut -f1,127- $csv|grep Tumour|sed 's/_/\t/g' |grep -v start |sed 's/"//g'|cut -f1-5 |awk -v n=${csv/.csv/} 'BEGIN{OFS="\t"}{print $1,$2,$3>n"."$4"."$5".bed"}'
#   csv=TAPSbeta.all_samples.Tissue_Group_DMBs.all_delta_quants0.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.1.csv
#   cut -f1,127- $csv|grep Hyper |grep Tumour|cut -f1|sed 's/_/\t/g' |grep -v start |sed 's/"//g'>${csv/.csv/}.hyper.bed
#   cut -f1,127- $csv|grep Hypo |grep Tumour|cut -f1|sed 's/_/\t/g' |grep -v start |sed 's/"//g'>${csv/.csv/}.hypo.bed
#   cut -f1,127- $csv|grep Tumour|sed 's/_/\t/g' |grep -v start |sed 's/"//g'|cut -f1-5 |awk -v n=${csv/.csv/} 'BEGIN{OFS="\t"}{print $1,$2,$3>n"."$4"."$5".bed"}'
#   
#   bed=(CAPS.all_samples.Tissue_Group_DMBs.all_delta_quants0.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.hyper.bed CAPS.all_samples.Tissue_Group_DMBs.all_delta_quants0.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.hypo.bed CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.hyper.bed CAPS.all_samples.Tissue_Group_DMBs.top300.CAPS_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg.hypo.bed TAPSbeta.all_samples.Tissue_Group_DMBs.all_delta_quants0.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.1.hyper.bed TAPSbeta.all_samples.Tissue_Group_DMBs.all_delta_quants0.TAPSbeta_hmm.all_samples.test7_2hiddenstates_mincov3_snp_reindex.Filtered_blocks.large_blocks_4cpg_delta_healthy_tumour_0.1.hypo.bed)
#   
#   faidx=/gpfs3/well/ludwig/users/cfo155/tissueMap/others/mC/resource/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.fai
#   gene=/gpfs3/well/ludwig/users/cfo155/tissueMap/others/mC/resource/MANE.GRCh38.v1.0.refseq_genomic.gene.bed
#   for csv in `ls *csv`
#   do
#   cut -f1,127- $csv|grep Tumour|sed 's/_/\t/g' |grep -v start |sed 's/"//g'|cut -f1-5 |\
#       bedtools sort -i - -g $faidx |\
#       bedtools closest -a - -b <(awk 'BEGIN{OFS="\t"}{if($6=="+")print $1,$2,$2+1,$4;else print $1,$3-1,$3,$4 }' $gene|\
#       bedtools sort -i - -g $faidx) -g $faidx -wa -wb -d|\
#       tee ${caps_dmb}.all_genelist|\
#       awk '$10 < 200000'|cut -f4,5,9|sort -u|awk -v n=${csv/.csv/} 'BEGIN{OFS="\t"}{print $3>"genelist/"n"_"$1"_"$2"_genelist"}'
#   done