#!/bin/bash
#SBATCH --job-name=atac-seq
#SBATCH -p short
#SBATCH --array  0-467:1 #0-1771:1 #
#SBATCH --requeue 

#### intersect atac-seq and DMB for tumour
# bws=($(find /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/atac-seq/tumour/ -name "*bw"))
# beds=($(ls /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr14/*tumour.bed))
# 
# bw_idx=$((SLURM_ARRAY_TASK_ID % ${#bws[@]}))
# bed_idx=$((SLURM_ARRAY_TASK_ID / ${#bws[@]}))
# bw="${bws[$bw_idx]}"
# bed="${beds[$bed_idx]}"
# 
# output=$(echo `basename $bed|sed 's/bed$//g;s/.all_samples.Tissue_Group_DMBs.top300.*.large_blocks_4cpg//g'` `basename $bw|sed 's/.bw$/.bed/g'` | sed 's/ //g')
# echo "$output"
# module purge
# module load BEDTools
# bigWigToBedGraph=/users/ludwig/cfo155/miniconda2/bin/bigWigToWig
# bedtools map \
#          -a <(sort -k1,1 -k2,2n $bed) \
#          -b <($bigWigToBedGraph  $bw stdout|sort -k1,1 -k2,2n) \
#          -c 4 \
#          -o mean \
#          -null 0 > $output

#### download atac-seq from encode for healthy tissues
# for i in `tail -n +2 encode.atac_seq.txt|cut -d '/' -f5|tr '\n' ' '`
# do
#     a=`curl -s https://www.encodeproject.org/files/${i}/|sed 's/<[^>]*>/\n/g'|grep -v ^$|grep -P "^Dataset|Output type" -A1|paste - - - - - |cut -f2,5`
#     b=`echo "$a"|cut -f1`
#     c=`curl -s https://www.encodeproject.org/experiments/${b}/|sed 's/<[^>]*>/\n/g'|grep -v ^$|grep -P "Biosample summary" -A2 |paste - - - |cut -f2,3`
#     echo -e "$i;$a;$c"
# done >encode.atac_seq.sampleinfo.txt
# 
# tail -n +2 encode.atac_seq.txt |while read line ;do wget -c $line ;done

# cut -d ';' -f1,3 encode.atac_seq.sampleinfo.txt |\
#     tail -n +2|sed 's/;Homo sapiens//g;s/tissue.*//g;s/female.*//g;s/male.*//g'|\
#     cut -f1,2|sed 's/\t /\t/g;s/ $//g;s/ /_/g'|\
#     awk 'BEGIN{OFS="\t"}{print "ln -s "$1".bigWig",$1"."$2".bigWig"}'

#### intersect atac-seq and DMB for normal ####
bws=($(ls /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/atac-seq/encode/*.*.bigWig ))
beds=($(ls /gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/dmb_call/tissue_atlas_v3_dmr14/*normal.bed ))

bw_idx=$((SLURM_ARRAY_TASK_ID % ${#bws[@]}))
bed_idx=$((SLURM_ARRAY_TASK_ID / ${#bws[@]}))
bw="${bws[$bw_idx]}"
bed="${beds[$bed_idx]}"

output=$(echo `basename $bed|sed 's/bed$//g;s/.all_samples.Tissue_Group_DMBs.top300.*.large_blocks_4cpg//g'` `basename $bw|sed 's/.bigWig$/.bed/g'` | sed 's/ //g')
echo "$output"
module purge
module load BEDTools
bigWigToBedGraph=/users/ludwig/cfo155/miniconda2/bin/bigWigToWig
bedtools map \
         -a <(sort -k1,1 -k2,2n $bed) \
         -b <($bigWigToBedGraph  $bw stdout|sort -k1,1 -k2,2n) \
         -c 4 \
         -o mean \
         -null 0 > $output