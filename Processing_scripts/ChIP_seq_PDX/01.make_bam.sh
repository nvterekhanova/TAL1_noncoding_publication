sample=$1
reference=$2
out_dir=$3

# trimming 
trim_galore --fastqc --cores 8 FASTQ/$sample.fq
mv $sample\_trimmed* $out_dir/.
mv $sample\.*_trimming_report* $out_dir/.

# mapping, use bwa mem
bwa mem $reference $out_dir/$sample\_trimmed.fq -t 10 | samtools view -b - > $out_dir/$sample.unmark.bam

# mark duplicates
cat $out_dir/$sample.unmark.bam | bamsormadup > $out_dir/$sample.bam

# remove duplicates and filter by quality
samtools view -F 1024 -b -q 1 $out_dir/$sample.bam > $out_dir/BAM/$sample.rmdupq1.bam
