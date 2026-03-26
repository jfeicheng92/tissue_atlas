# https://www.ebi.ac.uk/gwas/docs/file-downloads gwas_catalog_v1.0-associations_e115_r2025-12-03_full.zip
gwas=gwas-catalog-download-associations-v1.0-full.tsv # gwas_catalog_v1.0-associations_e115_r2025-12-03_full.zip
cat $gwas|awk 'BEGIN{FS="\t";OFS="\t"}{print $12,$13,$22,$8,$6}'|grep -iP "cancer|carcinoma"  > ${gwas}.cancer.txt
cat ${gwas}.cancer.txt|awk 'BEGIN{FS="\t"}{if($1~/[0-9]/ && $1!~/[;,x]/)print $0}'|sed 's/^/chr/g' >${gwas}.cancer.select.txt
cat ${gwas}.cancer.txt|awk 'BEGIN{FS="\t"}{if($1~/[;]/)print $0}'|\
awk -F "\t" '{
    n = split($1, a, ";");
    split($2, b, ";");
    split($3, c, ";");
    for (i = 1; i <= n; i++)
        printf "%s\t%s\t%s\t%s\n", a[i], b[i], c[i], $4, $5
  }'|sed 's/^/chr/g'  >>${gwas}.cancer.select.txt
cat ${gwas}.cancer.txt|awk 'BEGIN{FS="\t"}{if($1~/x/)print $0}'|\
awk -F "\t" '{
    n = split($1, a, " x ");
    split($2, b, " x ");
    split($3, c, " x ");
    for (i = 1; i <= n; i++)
        printf "%s\t%s\t%s\t%s\n", a[i], b[i], c[i], $4, $5
  }'|sed 's/^/chr/g' >>${gwas}.cancer.select.txt
# if the first colum is empty($1==""), Variant does not map to the genome




cut -f4 gwas-catalog-download-associations-v1.0-full.tsv.cancer.select.txt |sort |uniq -c|sed 's/^ *//g;s/ /\t/'|awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1}' >cancer_summary.txt
gwas=gwas-catalog-download-associations-v1.0-full.tsv # gwas_catalog_v1.0-associations_e115_r2025-12-03_full.zip
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Breast cancer")print $0,"breast"}' ${gwas}.cancer.select.txt >breast_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Kidney cancer")print $0,"kidney"}' ${gwas}.cancer.select.txt >kidney_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Hepatic cancer"||$4=="Hepatocellular carcinoma")print $0,"liver"}' ${gwas}.cancer.select.txt >liver_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Lung cancer")print $0,"lung"}' ${gwas}.cancer.select.txt >lung_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Epithelial ovarian cancer") print $0,"ovary"}' ${gwas}.cancer.select.txt >ovary_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Pancreatic cancer")print $0,"pancreas"}' ${gwas}.cancer.select.txt >pancreas_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Prostate cancer")print $0,"prostate"}' ${gwas}.cancer.select.txt >prostate_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Colorectal cancer" && $5=="www.ncbi.nlm.nih.gov/pubmed/36539618")print $0,"colon"}' ${gwas}.cancer.select.txt >colon_${gwas}.cancer.select.txt
awk 'BEGIN{OFS="\t";FS="\t"}{if($4=="Gastric cancer")print $0,"stomach"}' ${gwas}.cancer.select.txt >stomach_${gwas}.cancer.select.txt



tumours=(breast colon kidney liver lung ovary pancreas prostate stomach)
genome=/gpfs3/well/ludwig/users/cfo155/tissueMap/methcalls/resource/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.fai


# DMR files (full paths). Order: umC then hmC (you can add/remove files)
dmr_files=(
  "../dmr_call_new1/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.umC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop2000.txt"
  "../dmr_call_new1/all_sample.merged.mlml.mincov10_common.groupby.hg38_ws1000.s500.hmC.bgQ0.05-0.05.tgQ0.25.hyper_dmrs.bg_quant_modegroups.alltop2000.txt"
)
flanking=20000
# Loop over tumour names
for t in "${tumours[@]}"; do
  # GWAS file expected to be named like: breast_gwas_catalog_...cancer.txt
  gwas="${t}_gwas-catalog-download-associations-v1.0-full.tsv.cancer.select.txt"


  # Build the tumour label used inside the DMR file (capitalize first letter + "-Tumour")
  # e.g. "kidney" -> "Kidney-Tumour"
  tumour_label="${t^}-Tumour"

  for dmr in "${dmr_files[@]}"; do
    if [[ ! -f "$dmr" ]]; then
      echo "Warning: DMR file not found: $dmr  — skipping" >&2
      continue
    fi

    # infer dmr type string for output filename
    if [[ "$dmr" == *umC* ]]; then
      dmr_type="umC"
    elif [[ "$dmr" == *hmC* ]]; then
      dmr_type="hmC"
    else
      dmr_type="dmr"
    fi

    out="${t}.${dmr_type}.gwas_catalog_v1.0_vs_dmr.fisher.txt"
    echo "Running: tumour=${tumour_label}, gwas=${gwas}, dmr=${dmr}, out=${out}"

    # Prepare GWAS windows: (skip header if present) create 3-col BED: chr, start-flanking, end+flanking
    # NOTE: adjust NR>1 filter depending on your GWAS file header structure
    awk -v f=$flanking 'BEGIN{OFS="\t"}{print $1, ($2-f)>=0?($2-f):0, $2+f}' "$gwas" \
      | awk '$1~/chr[0-9].*/' |awk '$1!~/;/' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      | bedtools fisher -a - \
        -b <( tail -n +2 "$dmr" \
               | awk -v t="$tumour_label" 'BEGIN{OFS="\t"} $6==t {print $1, $2, $3}' \
               | sort -k1,1 -k2,2n ) \
        -g <( grep -E '^chr([0-9]+|X)\b' "$genome" | sort -k1,1 -k2,2n ) \
      > "$out"
    awk -v f=$flanking 'BEGIN{OFS="\t"}{print $1, ($2-f)>=0?($2-f):0, $2+f,$0}' "$gwas" \
      | awk '$1~/chr[0-9].*/' |awk '$1!~/;/' \
      | sort -k1,1 -k2,2n \
      | bedtools intersect -a - \
        -b <( tail -n +2 "$dmr" \
               | awk -v t="$tumour_label" 'BEGIN{OFS="\t"} $6==t {print $1, $2, $3}' \
               | sort -k1,1 -k2,2n ) \
      -wa -wb > ${out/.txt/}.intersection.txt
    echo "Wrote: $out"
  done
done


tail -n 1 *.gwas_catalog_v1.0_vs_dmr.fisher.txt|paste - - -|\
  sed 's/.gwas_catalog_v1.0_vs_dmr.fisher.txt <==//g;s/==> //g'|\
  cut -f1,4,5|awk 'BEGIN{OFS="\t"}{print $1,$2,$3}'|cat <(echo -e "dmr\ttwo-tail\tratio" ) -  >gwas_catalog_v1.0_vs_dmr.fisher.txt


