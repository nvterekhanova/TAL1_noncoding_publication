sample=$1
chr_sizes=$2


samtools flagstat out/$sample.bam > out/QC/$sample.flagstat

scalefactor=`awk 'NR==5{print $1} ' out/QC/$sample.flagstat | awk '{print 100000000/$1}' `

echo $sample
echo $scalefactor

bedtools genomecov -bga -ibam out/$sample.bam -scale $scalefactor -split > out/coverage/$sample.bedgraph

awk 'BEGIN{OFS="\t"}{print "chr"$0}' out/coverage/$sample.bedgraph |grep -v chrGL |grep -v chrMT > out/coverage/$sample.chr.bedgraph

sort -k1,1 -k2,2n out/coverage/$sample.chr.bedgraph > out/coverage/$sample.chr.sorted.bedgraph

bedGraphToBigWig out/coverage/$sample.chr.sorted.bedgraph $chr_sizes out/BigWig/$sample.bw

