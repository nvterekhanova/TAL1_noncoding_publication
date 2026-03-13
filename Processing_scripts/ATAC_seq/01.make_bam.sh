file=$1
reference=$2
out_dir=$3
sample=$4

# trimming
trim_galore --fastqc --cores 8 --output_dir $out_dir/. --paired FASTQ/$sample/$file\_R1.fq FASTQ/$sample/$file\_R2.fq  

# mapping, use bwa mem
bwa mem $reference $out_dir/$file\_R1_val_1.fq $out_dir/$file\_R2_val_2.fq -t 10 | samtools view -b - > $out_dir/$file.unmark.bam

# mark duplicates
cat $out_dir/$file.unmark.bam | bamsormadup > $out_dir/$file.bam

# remove duplicates, filter by quality and remove singletons
samtools view -F 1804 -b -q 1 $out_dir/$file.bam > $out_dir/BAM/$file.rmdupq1.bam
