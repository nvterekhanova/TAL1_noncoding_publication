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
library(stringr)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())


# read annotation:
tab <- read.table('../5.Drivers/20250620_Cohort_overview/out/TAL1_alterations_Allsamples_with_TAL1_WT.annotated.20250623.tsv',sep='\t',header=T)
tab_s <- tab[,c('USI','SJ_ID','Group','Reviewed.subtype')]
tab_s <- tab_s[!is.na(tab_s$Group) & tab_s$Group!='TAL1 WT',]


to_plot_s2 <- fread('out/Table_summary_134_samples_with_predictions.20260304.tsv',sep='\t',header=T)
to_plot_s2 <- as.data.frame(to_plot_s2)
to_plot_s2 <- to_plot_s2[to_plot_s2$Category!='ITD',]

# also make only for samples from this study:
to_plot_s2 <- to_plot_s2[to_plot_s2$Study=='Polonen_et_al',]
to_plot_s2$SJ_ID <- gsub('^(.*_D)[0-9]','\\1',to_plot_s2$SJ_ID)
alt_order <- c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh','TAL1 Downstream Enh I','TAL1 Downstream Enh II')
cols_id <- brewer.pal(name='Dark2',n=6)[c(1,6,5,2,4,3)]
names(cols_id) <- alt_order
cat_cols <- c('SNV'='#737373','Indel, 1bp'='#984ea3','Indel, 2-3bp'='#4daf4a','Indel, 4-9bp'='#377eb8','Indel, >=10bp'='#ff7f00', 'ITD'='#e41a1c')

#this table is needed for a separate layer without Violins for single data points
to_plot_s2_filt <- to_plot_s2[to_plot_s2$Group!='TAL1 Intron Enh',]
p <- ggplot()
p <- p + geom_violin(data=to_plot_s2_filt, aes(x=Group,y=tal1_diff_in_cd34,color=Group),width=0.7, trim=T, drop=F)
p <- p + stat_summary(data=to_plot_s2_filt, aes(x=Group,y=tal1_diff_in_cd34),fun = median, fun.min = median, fun.max = median, geom = "crossbar", width=0.8)
p <- p + geom_jitter(data=to_plot_s2, aes(x=Group,y=tal1_diff_in_cd34,fill=Category,size=Var_count),
position=position_jitter(width=0.2,seed=124,height=0),color='black',pch=21, alpha=0.8)
p <- p + scale_color_manual(values=cols_id) + guides(fill=guide_legend(override.aes=list(size=5,pch=21),title='Group'))
p <- p + scale_fill_manual(values=cat_cols)
#p <- p + facet_grid(.~Study, space='free',scales='free')
p <- p + theme_bw() + theme_nogrid() + labs(title="",x="",y="Predicted TAL1 RNA-seq expr. score by AlphaGenome")
p <- p + theme(axis.text.x = element_text(colour="black", size=10, angle = 45, hjust = 1),
axis.text.y = element_text(colour="black", size=10))
p <- p + geom_text_repel(data=to_plot_s2,aes(y=tal1_diff_in_cd34,x=Group,label=ifelse(USI %in% c('PAUAZV','PARSNX'), as.character(SJ_ID),NA)),size=2,
segment.size=0.1, box.padding = 0.8, min.segment.length = 0)
p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
p <- p + theme(strip.text.x.top = element_text(angle = 90, hjust = 0.05, vjust = 0.5))
p <- p + scale_size_continuous(range  = c(0.1, 7),
                         limits = c(0, max(to_plot_s2$Var_count)),
                         breaks = c(1, 5, 10, 15, 20, 25))

pdf('plots/Violin_plot_predChange_inTAL1_expr.byIndelSize.Polonen_et_al_study.withLabs.no_ITDs.20260304.pdf',width=5,height=4,useDingbats=F)
print(p)
dev.off()




# Make plot for TAL1 expression:
to_plot <- fread('out/Table_122_samples_with_predictions.20260304.tsv',sep='\t',header=T)
to_plot <- as.data.frame(to_plot)

# make violin plot for gene expression
to_plot_s <- to_plot
to_plot_s$Category <- 'SNV'
to_plot_s$Category <- ifelse(to_plot_s$Var_size==2, 'Indel, 1bp', to_plot_s$Category)
to_plot_s$Category <- ifelse(to_plot_s$Var_size>=3 & to_plot_s$Var_size<=4,"Indel, 2-3bp",to_plot_s$Category)
to_plot_s$Category <- ifelse(to_plot_s$Var_size>=5 & to_plot_s$Var_size<=10,"Indel, 4-9bp",to_plot_s$Category)
to_plot_s$Category <- ifelse(to_plot_s$Var_size>=11,"Indel, >=10bp",to_plot_s$Category)
# all ITDs are >40bp
to_plot_s$Category <- ifelse(to_plot_s$Var_size>40, 'ITD', to_plot_s$Category)

to_plot_s$SJ_ID <- gsub('^(.*_D)[0-9]','\\1',to_plot_s$SJ_ID)
to_plot_s <- to_plot_s[to_plot_s$Category!='ITD',]

p <- ggplot(data=to_plot_s, aes(x=Group,y=Log2_expr))
p <- p + geom_violin(width=0.7,aes(color=Group), trim=T)
p <- p + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width=0.8)
p <- p + geom_jitter(position=position_jitter(width=0.2,height=0,seed=123),aes(fill=Category),size=1.5,color='black',pch=21, alpha=0.8)
p <- p + scale_color_manual(values=cols_id) + guides(fill=guide_legend(override.aes=list(size=5,pch=21),title='Group'))
p <- p + scale_fill_manual(values=cat_cols)
p <- p + theme_bw() + theme_nogrid() + labs(title="",x="",y="Observed TAL1 RNA expr. in patient samples")
p <- p + theme(axis.text.x = element_text(colour="black", size=10, angle = 45, hjust = 1),
axis.text.y = element_text(colour="black", size=10))
p <- p + geom_text_repel(aes(y=Log2_expr,x=Group,label=ifelse(USI %in% c('PAUAZV','PARSNX'), as.character(SJ_ID),NA)),size=2,
segment.size=0.1, box.padding = 0.8, min.segment.length = 0)
p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
p <- p + theme(strip.text.x.top = element_text(angle = 90, hjust = 0.05, vjust = 0.5))

pdf('plots/Violin_plot_Epression_inTAL1.byIndelSize.withLabs.no_ITDs20260304.pdf',width=4.8,height=4,useDingbats=F)
print(p)
dev.off()

