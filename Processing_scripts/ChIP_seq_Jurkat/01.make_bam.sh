sample=$1
out_dir=$2
index=$3

# trimming
trim_galore --fastqc --cores 8 FASTQ/$sample.fastq.gz
mv $sample\_trimmed* $out_dir/.
mv $sample\.*_trimming_report* $out_dir/.

# mapping
bowtie2 --rfg 1,1 -p 10 -k 1 -q -x $index -U $out_dir/$sample\_trimmed.fq.gz -S $out_dir/$sample.unmark.sam

# convert .sam to .bam
samtools view -bS $out_dir/$sample.unmark.sam > $out_dir/$sample.unmark.bam

# mark duplicates
cat $out_dir/$sample.unmark.bam | bamsormadup > $out_dir/$sample.bam

# remove duplicates and filter by quality
samtools view -F 1024 -b -q 1 $out_dir/$sample.bam > $out_dir/BAM/$sample.rmdupq1.bam
