# Merge BAM files across replicates

module load samtools/1.9

for sample in SJTALL002048 SJALL015708

do
	
    for file in input/$sample/*.bam; do echo $file >> input/Lists/$sample\_list.tsv; done; 

    bsub -P PCGP -q compbio -J merge_bam -n 1 -R 'rusage[mem=20000]' -eo logs/01.merge_bam.$sample.err -oo logs/01.merge_bam.$sample.out "samtools merge -b input/Lists/$sample\_list.tsv out/$sample\_X.bam; samtools index out/$sample\_X.bam"

done
