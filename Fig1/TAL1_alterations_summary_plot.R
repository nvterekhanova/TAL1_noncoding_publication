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


# read annotation:
tab <- read.table('out/TAL1_alterations_Allsamples_with_TAL1_WT.annotated.20250623.tsv',sep='\t',header=T)
tab_s <- tab[,c('USI','SJ_ID','Group')]
tab_s <- tab_s[!is.na(tab_s$Group) & tab_s$Group!='TAL1 WT',]


alt_order <- c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh','TAL1 Downstream Enh I','TAL1 Downstream Enh II')
cols_id <- brewer.pal(name='Dark2',n=6)[c(1,6,5,2,4,3)]
names(cols_id) <- alt_order


####################################################
#make a simple barplot listing all TAL1 alterations#
####################################################

to_plot <- tab_s

to_plot$Group <- factor(to_plot$Group, levels=alt_order)

to_plot$Count <- 1

to_plot_2 <- aggregate(to_plot$Count, by=list(to_plot$Group), FUN='sum')
colnames(to_plot_2) <- c('Group','Count')

to_plot_2 <- to_plot_2[order(to_plot_2$Count),]
to_plot_2$Group <- factor(to_plot_2$Group,levels=to_plot_2$Group)

p <- ggplot(to_plot_2, aes(x=Group, y=Count, fill=Group))+ geom_bar(stat='identity',color='black')

p <- p + scale_fill_manual(values=cols_id)

p <- p + theme_classic() + ylab('N of cases with TAL1 alteration') + xlab('')

p <- p + coord_flip()

p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=12),
axis.text.y = element_text(colour="black", size=12))
p <- p + theme(legend.position='none')

pdf(paste("plots/TAL1_alter_summary_barplot.v2.20250623.pdf",sep=""),width=5,height=2.2,useDingbats=F)
p
dev.off()
