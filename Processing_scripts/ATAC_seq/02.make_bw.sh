sample=$1
genome_sizes=$2
out_dir=$3

samtools flagstat $out_dir/BAM/$sample.rmdupq1.sorted.bam > $out_dir/QC/$sample.flagstat

scalefactor=`awk 'NR==5{print $1} ' $out_dir/QC/$sample.flagstat | awk '{print 1000000000/$1}' `

echo $sample
echo $scalefactor


# nucleosome free
cat $out_dir/BED/$sample.rmdupq1.sorted.bedpe | awk 'BEGIN{OFS="\t";} $NF == "-" && $(NF-1) == "+" && $1 == $4 && $6 - $2 < 109 && $2 + $6 > 80 && $6 > $2 {mid=int(($2+$6)/2); print "chr"$1,mid-40,mid+41,$7,$8,"+"}' | grep -v '^chrGL' | grep -v '^chrMT' > $out_dir/BED/$sample.rmdupq1.free.bed # 109 == 100 + Tn5 adjustment

sort-bed --max-mem 20G $out_dir/BED/$sample.rmdupq1.free.bed > $out_dir/BED/$sample.rmdupq1.sorted.free.bed

bedtools genomecov -bg -i $out_dir/BED/$sample.rmdupq1.sorted.free.bed -g $genome_sizes -scale $scalefactor > $out_dir/Coverage/$sample.rmdupq1.free.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n $out_dir/Coverage/$sample.rmdupq1.free.bedGraph > $out_dir/Coverage/$sample.rmdupq1.free.sorted.bedGraph

bedGraphToBigWig $out_dir/Coverage/$sample.rmdupq1.free.sorted.bedGraph $genome_sizes $out_dir/Coverage/$sample.nuc_free.bw

# nucleosome free and mono/di/tri-nucleosomes
cat $out_dir/BED/$sample.rmdupq1.sorted.bedpe | awk 'BEGIN{OFS="\t";} $NF == "-" && $(NF-1) == "+" && $1 == $4  && $2 + $6 > 80 && $6 > $2 {mid=int(($2+$6)/2); print "chr"$1,mid-40,mid+41,$7,$8,"+"}' | grep -v '^chrGL' | grep -v '^chrMT' > $out_dir/BED/$sample.rmdupq1.all.bed

sort-bed --max-mem 20G $out_dir/BED/$sample.rmdupq1.all.bed > $out_dir/BED/$sample.rmdupq1.sorted.all.bed

bedtools genomecov -bg -i $out_dir/BED/$sample.rmdupq1.sorted.all.bed -g $genome_sizes -scale $scalefactor > $out_dir/Coverage/$sample.rmdupq1.all.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n $out_dir/Coverage/$sample.rmdupq1.all.bedGraph > $out_dir/Coverage/$sample.rmdupq1.all.sorted.bedGraph

bedGraphToBigWig $out_dir/Coverage/$sample.rmdupq1.all.sorted.bedGraph $genome_sizes $out_dir/Coverage/$sample.all.bw


