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

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())


##########################################################################
#Make a barplot showing subtype distribution across different TAL1 groups#
##########################################################################

alt_order <- c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh','TAL1 Downstream Enh I','TAL1 Downstream Enh II')
cols_id <- brewer.pal(name='Dark2',n=6)[c(1,6,5,2,4,3)]
names(cols_id) <- alt_order

subt_cols <- brewer.pal(name='Set1',n=4)[c(2,4)]
names(subt_cols) <- c('TAL1 DP-like','TAL1 alpha/beta-like')

tab <- read.table('out/TAL1_alterations_Allsamples_with_TAL1_WT.annotated.20250623.tsv',sep='\t',header=T)
tab_s <- tab[,c('USI','SJ_ID','Group','Reviewed.subtype')]
tab_s <- tab_s[!is.na(tab_s$Group) & tab_s$Group!='TAL1 WT',]

tab_s$Count <- 1
tab_s1 <- aggregate(tab_s$Count, by=list(tab_s$Group,tab_s$Reviewed.subtype), FUN='sum')
tab_s2 <- aggregate(tab_s$Count, by=list(tab_s$Group), FUN='sum')

colnames(tab_s1) <- c('Group','Reviewed.subtype','Count')
colnames(tab_s2) <- c('Group','Total_count')

to_plot <- merge(tab_s1,tab_s2)
to_plot$Fraction <- to_plot$Count/to_plot$Total_count


to_plot$Group <- factor(to_plot$Group, levels=alt_order)

p <- ggplot(to_plot, aes(x=Reviewed.subtype, y=Fraction*100, fill=Reviewed.subtype))+ geom_bar(stat='identity',color='black')

p <- p + scale_fill_manual(values=subt_cols)

p <- p + theme_classic() + ylab('% of cases') + xlab('')

p <- p + facet_grid(.~Group, space='free',scales='free') + ylim(0,100)

p <- p + theme_bw() + theme_nogrid() + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
p <- p + theme(strip.text.x.top = element_text(angle = 90, hjust = 0.05, vjust = 0.5, size=13))
p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
axis.text.y = element_text(colour="black", size=12),axis.ticks.x = element_blank(), legend.position='bottom')
pdf(paste("plots/TAL1_subt_across_groups_barplot.20250721.pdf",sep=""),width=6,height=4,useDingbats=F)
p
dev.off()

write.table(to_plot, 'out/Samples_bySubtype_byTAL1group_forBarPlot.20250721.tsv',sep='\t',quote=F,row.names=F)
