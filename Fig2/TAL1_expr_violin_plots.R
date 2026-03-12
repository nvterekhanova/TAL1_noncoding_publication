library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(data.table)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(readxl)

# read annotation:
tab <- read.table('out/TAL1_alterations_Allsamples_with_TAL1_WT.annotated.20250623.tsv',sep='\t',header=T)
tab_s <- tab[,c('USI','SJ_ID','Group')]
tab_s <- tab_s[!is.na(tab_s$Group) & tab_s$Group!='TAL1 WT',]

# read gene expression:
data_dir='/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data'

expr <- read.table(paste(data_dir,'/RNAseq/Polonen_2024_samples/From_synapse/TALL_X01_tpm.tsv',sep=''),sep='\t',header=T)
gene_annot <- fread(paste(data_dir,'/RNAseq/Polonen_2024_samples/From_synapse/TALL_X01_gene_annotations_unfiltered.tsv',sep=''))
gene_annot <- as.data.frame(gene_annot)
gene_annot_s <- gene_annot[gene_annot$bioType %in% c('protein_coding','TR_C_gene'),]
gene_annot_s <- gene_annot_s[,1:2]

expr_2 <- reshape2::melt(expr, id="X")
colnames(expr_2) <- c('geneID','Sample','Expr')
expr_2$Log2_expr <- log2(expr_2$Expr + 1)

# keep only protein-coding & TR_C genes
expr_2 <- expr_2[expr_2$geneID %in% gene_annot_s$geneID,]
expr_3 <- merge(expr_2,gene_annot_s)

expr_4 <- expr_3[expr_3$geneSymbol=='TAL1',]
colnames(expr_4)[2] <- 'USI'

alt_order <- c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh','TAL1 Downstream Enh I','TAL1 Downstream Enh II')
cols_id <- brewer.pal(name='Dark2',n=6)[c(1,6,5,2,4,3)]
names(cols_id) <- alt_order

to_plot <- merge(tab_s,expr_4)
to_plot$Group <- factor(to_plot$Group, levels=rev(alt_order))


# sum(table(to_plot$Group))==nrow(to_plot)


theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

options(ggrepel.max.overlaps = Inf)

write.table(to_plot, 'out/For_violin_plot_TAL1_expr_byGroup.20250623.tsv',sep='\t',quote=F,row.names=F)

p <- ggplot(data=to_plot, aes(x=Group,y=Log2_expr))
p <- p + geom_violin(width=0.7,aes(color=Group), trim=T)
p <- p + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width=0.8)
p <- p + geom_jitter(position=position_jitter(0.2),size=0.5,fill='black',pch=19, alpha=0.5)

p <- p + geom_text_repel(aes(y=Log2_expr,x=Group,label=ifelse(USI %in% c('PAUAZV','PARSNX'), as.character(SJ_ID),NA)),size=2,
segment.size=0.1, box.padding = 0.8, min.segment.length = 0)

p <- p + scale_color_manual(values=cols_id) + guides(fill=guide_legend(override.aes=list(size=5,pch=21),title='Group'))

p <- p + theme_bw() + theme_nogrid() + labs(title="",x="",y="TAL1 expression, TPM")

p <- p + theme(axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10)) 

p <- p + scale_y_continuous(breaks = seq(from=0, to=6, by=1), labels = c(0, 1, 3, 7, 15, 31, 63))

p <- p + coord_flip()

p <- p + theme(legend.position='none') +ggtitle("")

pdf('plots/TAL1_expr_groupedByTAL1_alter20250623.pdf',width=4,height=3.5,useDingbats=F)
print(p)
dev.off()
