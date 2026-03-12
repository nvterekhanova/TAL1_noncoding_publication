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
library(ggbreak)


vers <- '20260302'

data_dir='/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data'


#First, plot all ITDs of actual length:

tab <- read.table('out/Downstream_enhancer_I_alterations_and_TFBSs_for_barplot.20260301.tsv',sep='\t',header=T)
tab_2 <- tab[tab$Type %in% c('ITD','TF'),]

# remove the first GATA3 site, as it is too far away
tab_2 <- tab_2[!(tab_2$Type=='TF' & tab_2$hg19_start==47669198),]

tab_2$SJ_ID_2 <- factor(tab_2$SJ_ID_2, levels=unique(tab_2$SJ_ID_2))

to_plot <- tab_2
p1 <- ggplot()
p1 <- p1 + geom_rect(data=to_plot, aes(xmin=hg19_start,xmax=hg19_end,ymin=0,ymax=1, color=Outline, fill=Type))
p1 <- p1 + geom_vline(xintercept=unique(to_plot$MYB_start), linetype="dashed",color="grey")
p1 <- p1 + scale_fill_manual(values=c('ITD'='blue','Indel'='blue','TF'='#EF4136'))
p1 <- p1 + facet_grid(Type + SJ_ID_2~.,  scales='free', space='free', switch='y')
p1 <- p1 + theme_minimal() + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines")) + ylab('Sample')
p1 <- p1 + scale_color_manual(values=c('Short'='black','Long'=NA, 'TF'='grey'), na.value="transparent") 
p1 <- p1 + theme(axis.text.y=element_blank()) +ggtitle('Downstream Enhancer | samples') + theme(strip.placement = "outside",legend.position='bottom',
panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),axis.ticks.y=element_blank(), axis.ticks.x = element_line(linewidth = 1,
colour = "black"))
p1 <- p1 + theme(strip.text.y.left = element_text(angle = 0, hjust = 0.05, vjust = 0.5)) + theme(legend.position='none')

pdf(paste('plots/Summary_enh_gain_I_ITD_samples.hg19_coords.Full_view.',vers,'.pdf',sep=''),width=5,height=5,useDingbats=F)
p1
dev.off()

#Now plot zoomed in region, and map all TFBS sites:
tab <- read.table('out/Downstream_enhancer_I_alterations_and_TFBSs_for_barplot.20260301.tsv',sep='\t',header=T)
tab_2 <- tab[tab$Type %in% c('ITD','TF'),]
# remove the first GATA3 site, as it is too far away
tab_2 <- tab_2[!(tab_2$Type=='TF' & tab_2$hg19_start==47669198),]

tab_2$new_start <- ifelse(tab_2$hg19_start <= 47669198, 47669198, tab_2$hg19_start)
tab_2$new_end <- ifelse(tab_2$hg19_end >= 47669493, 47669493, tab_2$hg19_end)
tab_2$MYB_end_2 <- 47669393

tab_2$SJ_ID_2 <- factor(tab_2$SJ_ID_2, levels=unique(tab_2$SJ_ID_2))
to_plot <- tab_2

p1 <- ggplot()
p1 <- p1 + geom_rect(data=to_plot, aes(xmin=new_start,xmax=new_end,ymin=0,ymax=1, color=Outline, fill=Type))
p1 <- p1 + geom_vline(xintercept=unique(to_plot$MYB_start), linetype="dashed",color="grey")
p1 <- p1 + geom_vline(xintercept=unique(to_plot$MYB_end_2), linetype="dashed",color="grey")
p1 <- p1 + scale_fill_manual(values=c('ITD'='blue','Indel'='blue','TF'='#EF4136'))
p1 <- p1 + facet_grid(Type + SJ_ID_2~.,  scales='free', space='free', switch='y')
p1 <- p1 + theme_minimal() + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines")) + ylab('Sample')
p1 <- p1 + scale_color_manual(values=c('Short'='black','Long'=NA, 'TF'='grey'), na.value="transparent")
p1 <- p1 + theme(axis.text.y=element_blank()) +ggtitle('Downstream Enhancer | samples') + theme(strip.placement = "outside",legend.position='bottom',
panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),axis.ticks.y=element_blank(), axis.ticks.x = element_line(linewidth = 1,
colour = "black"))
p1 <- p1 + theme(strip.text.y.left = element_text(angle = 0, hjust = 0.05, vjust = 0.5)) + theme(legend.position='none')

pdf(paste('plots/Summary_enh_gain_I_ITD_samples.hg19_coords.Zoomed_in.',vers,'.pdf',sep=''),width=5,height=5,useDingbats=F)
p1
dev.off()
