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


pr <- fread('input/Model_predictions.v5.20260302.txt')
pr <- as.data.frame(pr)
pr$USI <- pr$variant_id
pr_s <- pr[pr$output=='TRUE',]
pr_s$List <- ifelse(row.names(pr_s) %in% c(1:32),'AlphaGenome','Polonen_et_al')
pr_s$Study <- ifelse(row.names(pr_s) %in% c(1:8),'Mansour_et_al','Polonen_et_al')
pr_s$Study <- ifelse(row.names(pr_s) %in% c(9:10),'Smith_et_al',pr_s$Study)
pr_s$Study <- ifelse(row.names(pr_s) %in% c(11:12),'Liu_et_al_2020',pr_s$Study)
pr_s$Study <- ifelse(row.names(pr_s) %in% c(13:32),'Liu_et_al',pr_s$Study)


# add 28 samples with TAL1 intron enhancers -- as they have the exactly same scores as those from Liu et al. (cis-X paper),
# because same sequence change
tab_intron <- tab_s[tab_s$Group=='TAL1 Intron Enh',]
intron_var_scores <- pr_s[11:12,]
all_tab_i <- NULL
for (usi in tab_intron$USI){
    tab_i <- intron_var_scores[1,]
    tab_i$USI <- usi
    all_tab_i <- rbind(all_tab_i, tab_i)
}
all_tab_i <- as.data.frame(all_tab_i)
all_tab_i$Study <- 'Polonen_et_al'
all_tab_i$List <- 'Polonen_et_al'


pr_s2 <- pr_s[pr_s$Study!='Liu_et_al',]
pr_s2 <- rbind(pr_s2,all_tab_i)
res <- merge(pr_s2,tab_s,all.x=T)

res$Group <- ifelse(res$Study=='Mansour_et_al','TAL1 Upstream Enh',res$Group)
res$Group <- ifelse(res$Study=='Smith_et_al','TAL1 Downstream Enh II',res$Group)
res$Group <- ifelse(res$Study=='Liu_et_al_2020','TAL1 Intron Enh',res$Group)

# plot them together:
to_plot <- res
to_plot$S_1 <- gsub('chr1\\:[0-9]+\\:(.*)>(.*)$','\\1',to_plot$variant)
to_plot$S_2 <- gsub('chr1\\:[0-9]+\\:(.*)>(.*)$','\\2',to_plot$variant)
to_plot$Size_1 <- str_length(to_plot$S_1)
to_plot$Size_2 <- str_length(to_plot$S_2)
to_plot$Var_size <- ifelse(to_plot$Size_1>=to_plot$Size_2,to_plot$Size_1,to_plot$Size_2)
to_plot_s <- to_plot[!duplicated(to_plot[,c('variant','Study')]),]

#add N of samples with the same variant:
to_plot$Count <- 1
counts <- aggregate(to_plot$Count, by=list(to_plot$variant,to_plot$Study,to_plot$Group), FUN='sum')
colnames(counts) <- c('variant','Study','Group','Var_count')

to_plot_s$Category <- 'SNV'
to_plot_s$Category <- ifelse(to_plot_s$Var_size==2, 'Indel, 1bp', to_plot_s$Category) 
to_plot_s$Category <- ifelse(to_plot_s$Var_size>=3 & to_plot_s$Var_size<=4,"Indel, 2-3bp",to_plot_s$Category)
to_plot_s$Category <- ifelse(to_plot_s$Var_size>=5 & to_plot_s$Var_size<=10,"Indel, 4-9bp",to_plot_s$Category)
to_plot_s$Category <- ifelse(to_plot_s$Var_size>=11,"Indel, >=10bp",to_plot_s$Category)
# all ITDs are >40bp
to_plot_s$Category <- ifelse(to_plot_s$Var_size>40, 'ITD', to_plot_s$Category)

to_plot_s2 <- merge(to_plot_s,counts)


#save table
write.table(to_plot_s2, 'out/Table_summary_134_samples_with_predictions.20260304.tsv',quote=F,row.names=F,sep='\t')

to_plot_s2 <- fread('out/Table_summary_134_samples_with_predictions.20260304.tsv',sep='\t',header=T)
to_plot_s2 <- as.data.frame(to_plot_s2)
to_plot_s2$Study <- factor(to_plot_s2$Study, levels=c('Polonen_et_al','Mansour_et_al','Liu_et_al_2020','Smith_et_al'))

cat_cols <- c('SNV'='#737373','Indel, 1bp'='#984ea3','Indel, 2-3bp'='#4daf4a','Indel, 4-9bp'='#377eb8','Indel, >=10bp'='#ff7f00', 'ITD'='#e41a1c')

alt_order <- c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh','TAL1 Downstream Enh I','TAL1 Downstream Enh II')
cols_id <- brewer.pal(name='Dark2',n=6)[c(1,6,5,2,4,3)]
names(cols_id) <- alt_order

subt_cols <- brewer.pal(name='Set1',n=4)[c(2,4)]
names(subt_cols) <- c('TAL1 DP-like','TAL1 alpha/beta-like')


# also make only for samples from this study:
to_plot_s2 <- to_plot_s2[to_plot_s2$Study=='Polonen_et_al',]
to_plot_s2$SJ_ID <- gsub('^(.*_D)[0-9]','\\1',to_plot_s2$SJ_ID)
to_plot_s2_filt <- to_plot_s2[to_plot_s2$Group!='TAL1 Intron Enh',]

p <- ggplot()
p <- p + geom_violin(data=to_plot_s2_filt, aes(x=Group,y=tal1_diff_in_cd34,color=Group),width=0.7, trim=T, drop=F)
p <- p + stat_summary(data=to_plot_s2_filt, aes(x=Group,y=tal1_diff_in_cd34),fun = median, fun.min = median, fun.max = median, geom = "crossbar", width=0.8)
p <- p + geom_jitter(data=to_plot_s2, aes(x=Group,y=tal1_diff_in_cd34,fill=Category,size=Var_count),
position=position_jitter(width=0.2,seed=123,height=0),color='black',pch=21, alpha=0.8)
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
pdf('plots/Violin_plot_predChange_inTAL1_expr.byVarSize.Polonen_et_al_study.withLabs.20260304.pdf',width=5,height=4,useDingbats=F)
print(p)
dev.off()




# also add epxression of TAL1:
tab_var <- to_plot[to_plot$Study=='Polonen_et_al',]

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

# keep only protein-coding  and TR_C genes
expr_2 <- expr_2[expr_2$geneID %in% gene_annot_s$geneID,]
expr_3 <- merge(expr_2,gene_annot_s)

colnames(expr_3)[2] <- 'USI'

expr_4 <- merge(expr_3, tab_s)
expr_4_s <- expr_4[expr_4$geneSymbol=='TAL1',]
expr_4_s <- expr_4_s[,c('USI','Log2_expr')]

# combine two tables:
res <- merge(tab_var,expr_4_s)

write.table(res, 'out/Table_122_samples_with_predictions.20260304.tsv',quote=F,row.names=F,sep='\t')

#################
#Start from here#
#################

to_plot <- fread('out/Table_122_samples_with_predictions.20260304.tsv',sep='\t',header=T)
to_plot <- as.data.frame(to_plot)

p <- ggplot(aes(Log2_expr,  tal1_diff_in_cd34),data=to_plot) 
p <- p + geom_point(aes(color=Group))
p <- p + geom_smooth(method="lm")+theme_classic()
p <- p + facet_wrap(Group~.,scales = "free")
#p <- p + facet_grid(.~Group,scales = "free",space="free")
p <- p + ggpubr::stat_cor(method = "pearson", label.y = -1)
p <- p + scale_color_manual(values=cols_id)
p <- p + theme(legend.position='bottom')

#pdf('plots/Scatterplot_cor_predChange_vs_expr.TAL1.v4.20260302.pdf',width=8,height=3.5,useDingbats=F)
#print(p)
#dev.off()

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

pdf('plots/Violin_plot_Epression_inTAL1.byVarSize.withLabs.20260304.pdf',width=4.8,height=4,useDingbats=F)
print(p)
dev.off()

