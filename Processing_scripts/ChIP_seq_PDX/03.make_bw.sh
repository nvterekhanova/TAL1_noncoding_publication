module load bedtools/2.31.0
module load ucsc/051223

sample=$1
genome_sizes=$2
out_dir=$3

frag_size=`awk 'NR==1{print $3} ' $out_dir/QC/CC/phantomPeakStatsReps.tab.$sample | awk -F ',' '{print $1}' `

echo $frag_size

bedtools bamtobed -i $out_dir/BAM/$sample.rmdupq1.bam | awk 'BEGIN{OFS="\t"}{print "chr"$0}' | grep -v '^chrGL' | grep -v '^chrMT' | awk '{if($6 ~ /-/){$2=$3-1;}else{$3=$2+1}; OFS="\t"; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | bedtools slop -s -l 0 -r $frag_size -g $genome_sizes -i - > $out_dir/BED/$sample.CHIP.bed


samtools flagstat $out_dir/BAM/$sample.rmdupq1.bam > $out_dir/QC/$sample.flagstat
scalefactor=`awk 'NR==5{print $1} ' $out_dir/QC/$sample.flagstat | awk '{print 100000000/$1}' `

echo $scalefactor

bedtools bedtobam -g $genome_sizes -i $out_dir/BED/$sample.CHIP.bed > $out_dir/BAM/$sample.CHIP.bam

samtools sort -o $out_dir/BAM/$sample.CHIP.sort.bam  $out_dir/BAM/$sample.CHIP.bam

bedtools genomecov -bga -ibam $out_dir/BAM/$sample.CHIP.sort.bam -scale $scalefactor > $out_dir/Coverage/$sample.CHIP.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n $out_dir/Coverage/$sample.CHIP.bedGraph > $out_dir/Coverage/$sample.CHIP.sorted.bedGraph

bedGraphToBigWig $out_dir/Coverage/$sample.CHIP.sorted.bedGraph $genome_sizes $out_dir/Coverage/$sample.bw

