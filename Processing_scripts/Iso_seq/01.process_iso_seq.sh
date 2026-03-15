
data_dir=/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data/Iso_seq
sample_dir=$1
sample_prefix=$2
sample=$3

ref_path=/research/rgs01/applications/hpcf/authorized_apps/cab/Automation/REF/Homo_sapiens/NCBI/GRCh37-lite/FASTA/GRCh37-lite_wchr.fa

#1. Primer removal and demultiplexing:
lima --isoseq --dump-clips --peek-guess -j 24 $data_dir/$sample_dir/$sample_prefix.hifi_reads.bam $data_dir/primers.fasta.fasta $sample/out/$sample_prefix.hifi.demult.bam

#2. Trimming polyA tails and contacamer removal
isoseq refine --require-polya -j 24 $sample/out/$sample_prefix.hifi.demult.NEB_5p--NEB_Clontech_3p.bam $data_dir/primers.fasta.fasta $sample/out/$sample_prefix.hifi.flnc.bam

#3. Isoforms
isoseq cluster -j 24 $sample/out/$sample_prefix.hifi.flnc.bam $sample/out/$sample_prefix.hifi.polished.bam --verbose --use-qvs

#4. Map to genome
pbmm2 align $ref_path $sample/out/$sample_prefix.hifi.polished.hq.bam $sample/out/$sample_prefix.hifi.polished.hq.aligned.bam -j 24 --preset ISOSEQ --sort --log-level INFO
