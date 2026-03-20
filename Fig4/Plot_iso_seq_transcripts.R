library(ggplot2)
library(ggtranscript)
library(dplyr)


samples <- c('Jurkat','SJTALL002048','SJALL015708')
all_class <- NULL
all_tab <- NULL
for (sample in samples){
class <- read.table(paste('../03.make_summary/input/',sample,'.classification.txt',sep=''),sep='\t',header=T)
class <- class[class$associated_gene=='TAL1',]
class_s <- class[,c('isoform','length','exons','fl_assoc')]
class_s$Sample <- sample
all_class <- rbind(all_class, class_s)

tab <- read.table(paste('../03.make_summary/input/',sample,'.collapsed.sorted.filtered_lite.gff',sep=''),sep='\t',header=F)
colnames(tab) <- c('chr','source','type','start','end','info','strand','info_2','transcript_name')
tab$Transcr_id <- gsub('^.*transcript_id (PB.*);$','\\1',tab$transcript_name)
tab_s <- tab[tab$type=='exon' & tab$Transcr_id %in% class$isoform,]
tab_s$Sample <- sample
all_tab <- rbind(all_tab,tab_s)
print(sample)
}

all_class <- all_class[order(-all_class$fl_assoc),]

write.table(all_tab, 'out/All_transcripts_exons_3_samples.20251231.tsv',sep='\t',quote=F,row.names=F)
write.table(all_class, 'out/All_transcripts_with_read_support_3_samples.20251231.tsv',sep='\t',quote=F,row.names=F)

# Start from here:
all_tab <- read.table('out/All_transcripts_exons_3_samples.20251231.tsv',sep='\t',header=T)
p <- ggplot(all_tab,aes(xstart = start, xend = end, y = transcript_name))
p <- p + geom_range(fill='blue')
p <- p + geom_intron(data = to_intron(all_tab, "transcript_name"),aes(strand = strand))
pdf(paste("plots/Iso-seq_transcripts_All_samples.ggtranscript.20251231.pdf",sep=""),width=8,height=6,useDingbats=F)
p
dev.off()

# 2026-01-01:
class <- read.table('out/All_transcripts_with_read_support_3_samples.20251231.tsv',sep='\t',header=T)
all_tab <- read.table('out/All_transcripts_exons_3_samples.20251231.tsv',sep='\t',header=T)
all_tab$Transcr_id <- factor(all_tab$Transcr_id, levels=rev(class$isoform))
all_tab$Sample <- factor(all_tab$Sample,levels=c('Jurkat','SJTALL002048','SJALL015708'))
all_tab$Type <- ifelse(all_tab$Sample=='Jurkat','Jurkat','PDX')
cols <- c('Jurkat'='#552D8D','PDX'='#931A1D')

p <- ggplot(all_tab,aes(xstart = start, xend = end, y = Transcr_id))
p <- p + geom_range(aes(fill=Type))
p <- p + geom_intron(data = to_intron(all_tab, "Transcr_id"),aes(strand = strand))
p <- p + scale_fill_manual(values=cols) 
p <- p + facet_grid(Sample~.,space='free',scales='free') + theme_minimal() + theme(legend.position='none')
pdf(paste("plots/Iso-seq_transcripts_All_samples.ggtranscript.20260101.pdf",sep=""),width=6,height=9.5,useDingbats=F)
p
dev.off()

# now also plot selected isoforms:
sel_iso <- c('PB.416.1','PB.416.12','PB.416.14','PB.657.2','PB.657.37','PB.683.6','PB.683.11')
all_tab_s <- all_tab[all_tab$Transcr_id %in% sel_iso,]
all_tab_s$Transcr_id <- factor(all_tab_s$Transcr_id,levels=rev(sel_iso))
p <- ggplot(all_tab_s,aes(xstart = start, xend = end, y = Transcr_id))
p <- p + geom_range(aes(fill=Type))
p <- p + geom_intron(data = to_intron(all_tab_s, "Transcr_id"),aes(strand = strand))
p <- p + scale_fill_manual(values=cols) + scale_x_continuous(breaks=seq(from=47682000, to=47698000, by=4000))
p <- p + facet_grid(Sample~.,space='free',scales='free') + theme_minimal() + theme(legend.position='none')
pdf(paste("plots/Iso-seq_Sel_transcripts_All_samples.ggtranscript.20260101.pdf",sep=""),width=6,height=1.2,useDingbats=F)
p
dev.off()
