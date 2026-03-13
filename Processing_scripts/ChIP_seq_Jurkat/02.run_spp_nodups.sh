sample=$1
out_dir=$2

# run QC using phantompeakqualtools
Rscript /home/nterekha/Software/phantompeakqualtools/run_spp_nodups.R -c=$out_dir/BAM/$sample.rmdupq1.bam -odir=$out_dir/QC/CC -savp -rf -out=$out_dir/QC/CC/phantomPeakStatsReps.tab.$sample
