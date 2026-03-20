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
library(ggnewscale)

data_dir='/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data'
theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())


res_st <- read.table('out/Fisher_test_results_for_gene_mut_enr_acrossGroups.20251021.tsv',sep='\t',header=T)

pathw <- read_excel('../../../Data/Pathways/JZ_format_pathway_genes_pancancer_genes_9_4_2025.xlsx',sheet=1)
pathw <- as.data.frame(pathw)
pathw_s <- pathw[pathw$disease_code=='TALL',]

# Do separately for oncogenes and tumor suppressors:
pathw_s2 <- pathw_s[pathw_s$gene_symbol %in% res_st$Targeted.gene,]
pathw_s2 <- pathw_s2[,c('gene_symbol','LOFofFunctionLOForGOFofFunctionG','ProposedPathway','disease_code')]

oncog <- pathw_s2$gene_symbol[pathw_s2$LOFofFunctionLOForGOFofFunctionG=='GOF' & !is.na(pathw_s2$LOFofFunctionLOForGOFofFunctionG)]
tsg <- pathw_s2$gene_symbol[pathw_s2$LOFofFunctionLOForGOFofFunctionG=='LOF' & !is.na(pathw_s2$LOFofFunctionLOForGOFofFunctionG)]

# try making heatmap:

# first, do it for oncogenes:
to_plot <- res_st[res_st$Fraction!=0 & res_st$Targeted.gene %in% oncog,]

to_plot$Count_group <- 1
freq_group <- aggregate(to_plot$Fraction, by=list(to_plot$Targeted.gene), FUN='sum')
colnames(freq_group) <- c('Gene','Frequency_group')
freq_group <- freq_group[order(freq_group$Frequency_group),]

to_plot$Pct <- round(to_plot$Fraction*100,1)
to_plot$Targeted.gene <- factor(to_plot$Targeted.gene, levels=freq_group$Gene)
to_plot$Group <- factor(to_plot$Group, levels=c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh',
'TAL1 Downstream Enh I','TAL1 Downstream Enh II'))

to_plot_suggest <- to_plot[to_plot$FDR>=0.05 & to_plot$FDR<0.15,]
to_plot_sign <- to_plot[to_plot$FDR<0.05,]

x <- colorRampPalette(append(brewer.pal(9,"YlOrRd"), "#FFFFFF", after=0), bias=2)(1000) 
p <- ggplot()
p <- p + geom_tile(data=to_plot,aes(x=Group, y=Targeted.gene, fill= Pct), linetype="blank") +
scale_fill_gradientn(name= "% Carriers",colours=x, na.value=NA, limit=c(0,NA))
p <- p + geom_text(data=to_plot,aes(x=Group, y=Targeted.gene, label = Pct), color="black", size=3)
p <- p + geom_tile(data=to_plot_suggest,aes(x=Group, y=Targeted.gene), color="grey",fill=NA, linewidth=1.5) 
p <- p + geom_tile(data=to_plot_sign,aes(x=Group, y=Targeted.gene), color="black",fill=NA, linewidth=1.5) 
p <- p  + theme_bw() + theme_nogrid() + scale_x_discrete(position = "top")
p <- p + theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust=0.2, hjust=0.05),
axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank()) +
theme(axis.text.x.top = element_text(vjust = 0))
pdf('plots/Alterations_top_byGroup.heatmap.oncogenes.20260209.pdf',width=6,height=4.5,useDingbats=F)
p
dev.off()

#############################################
# Now also make the same analysis for TSGs: #
#############################################

to_plot <- res_st[res_st$Fraction!=0 & res_st$Targeted.gene %in% tsg,]

to_plot$Count_group <- 1
freq_group <- aggregate(to_plot$Fraction, by=list(to_plot$Targeted.gene), FUN='sum')
colnames(freq_group) <- c('Gene','Frequency_group')
freq_group <- freq_group[order(freq_group$Frequency_group),]

to_plot$Pct <- round(to_plot$Fraction*100,1)
to_plot$Targeted.gene <- factor(to_plot$Targeted.gene, levels=freq_group$Gene)
to_plot$Group <- factor(to_plot$Group, levels=c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh',
'TAL1 Downstream Enh I','TAL1 Downstream Enh II'))

to_plot_suggest <- to_plot[to_plot$FDR>=0.05 & to_plot$FDR<0.15,]
to_plot_sign <- to_plot[to_plot$FDR<0.05,]

x <- colorRampPalette(append(brewer.pal(9,"YlGnBu"), "#FFFFFF", after=0), bias=2)(1000) 
p <- ggplot()
p <- p + geom_tile(data=to_plot,aes(x=Group, y=Targeted.gene, fill= Pct), linetype="blank") +
scale_fill_gradientn(name= "% Carriers",colours=x, na.value=NA, limit=c(0,NA))
p <- p + geom_text(data=to_plot,aes(x=Group, y=Targeted.gene, label = Pct), color="black", size=3)
p <- p + geom_tile(data=to_plot_suggest,aes(x=Group, y=Targeted.gene), color="grey",fill=NA, linewidth=1.5) 
p <- p + geom_tile(data=to_plot_sign,aes(x=Group, y=Targeted.gene), color="black",fill=NA, linewidth=1.5) 
p <- p  + theme_bw() + theme_nogrid() + scale_x_discrete(position = "top")
p <- p + theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust=0.2, hjust=0.05),
axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank()) +
theme(axis.text.x.top = element_text(vjust = 0))
pdf('plots/Alterations_top_byGroup.heatmap.tsgs.20260209.pdf',width=6,height=3.9,useDingbats=F)
p
dev.off()

