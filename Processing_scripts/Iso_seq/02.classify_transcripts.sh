#https://isoseq.how/classification/workflow.html

data_dir=/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data/Iso_seq
sample_dir=$1
sample_prefix=$2
sample=$3

ref=/home/nterekha/Software/Refs/GRCh37-lite_wchr/GRCh37-lite_wchr.fa
gtf=/home/nterekha/Software/Refs/gencode.v19.annotation/gencode.v19.annotation.sorted.gtf

mkdir $sample/out/classify_out

# 1. Collapse into unique isoforms
isoseq collapse --do-not-collapse-extra-5exons $sample/out/$sample_prefix.hifi.polished.hq.aligned.GRCh37-lite.bam $sample/out/$sample_prefix.hifi.flnc.bam $sample/out/$sample_prefix.collapsed.gff

# 2. Prepare input transcript GFF
pigeon prepare $sample/out/$sample_prefix.collapsed.gff

# Classify Isoforms
# 3. Transcript classification
#pigeon classify --out-dir $sample/out/classify_out $sample/out/$sample_prefix.collapsed.sorted.gff $gtf $ref

# 3. Adding FLNC counts to classification output
pigeon classify --fl $sample/out/$sample_prefix.collapsed.flnc_count.txt --out-dir $sample/out/classify_out $sample/out/$sample_prefix.collapsed.sorted.gff $gtf $ref


# 4. Filter isoforms
pigeon filter --isoforms $sample/out/$sample_prefix.collapsed.sorted.gff $sample/out/classify_out/$sample_prefix\_classification.txt


# 5. Report gene saturation
pigeon report $sample/out/classify_out/$sample_prefix\_classification.filtered_lite_classification.txt $sample/out/classify_out/$sample_prefix\_scisoseq_saturation.txt
