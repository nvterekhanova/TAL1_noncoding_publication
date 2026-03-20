library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(ggnewscale)
library(ggrastr)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())


# read annotation:
tab <- read.table(paste('/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Analysis/5.Drivers/',
'20250620_Cohort_overview/out/TAL1_alterations_Allsamples_with_TAL1_WT.annotated.20250623.tsv',sep=''),sep='\t',header=T)
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

# keep only protein-coding  and TR_C genes
expr_2 <- expr_2[expr_2$geneID %in% gene_annot_s$geneID,]
expr_3 <- merge(expr_2,gene_annot_s)

colnames(expr_3)[2] <- 'USI'

expr_4 <- merge(expr_3, tab_s)

expr_4$Group_2 <- ifelse(expr_4$Group=='TAL1 Downstream Enh I','TAL1 Downstream Enh I','Other')

# add subtype information
subt_info <- read.table(paste('/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Analysis/5.Drivers/',
'20250620_Cohort_overview/out/Subtype_info_sheet_4.tsv',sep=''),sep='\t',header=T)
subt_info_s <- subt_info[,c('sample','Reviewed.subtype')]
colnames(subt_info_s)[1] <- 'USI'
subt_info_s$Reviewed.subtype <- ifelse(subt_info_s$Reviewed.subtype=='TAL1 <U+03B1><U+03B2>-like','TAL1 alpha/beta-like',
subt_info_s$Reviewed.subtype)
subt_info_s$Reviewed.subtype <- ifelse(subt_info_s$Reviewed.subtype=='LMO2 <U+03B3><U+03B4>-like','LMO2 gamma/delta-like',
subt_info_s$Reviewed.subtype)

expr_5 <- merge(expr_4, subt_info_s)

#final counts of samples:
#test <- expr_5[expr_5$geneSymbol=='TAL1',]
#TAL1 Downstream Enh I  Other
#        	    25    366


#Run diff. expr. analysis (takes ~30 min):
all_stat <- NULL
for (gene in unique(expr_5$geneID)){
    tab <- expr_5[expr_5$geneID==gene,]
    geneSymbol <- unique(tab$geneSymbol)
    test <- glm((tab$Log2_expr ~ tab$Group_2 + tab$Reviewed.subtype))
    st <- summary(test)$coefficients
    st_g <- st[2,4]
    st_s <- st[3,4]
    est_g <- st[2,1]
    est_s <- st[3,1]	  
    stat <- cbind(gene,geneSymbol,st_g,est_g,st_s,est_s)
    all_stat <- rbind(all_stat,stat)
}

all_stat <- as.data.frame(all_stat)
colnames(all_stat) <- c('geneID','geneSymbol','P_value_g', 'Est_g','P_value_subt','Est_subt')
all_stat$P_value_g <- as.numeric(as.character(unlist(all_stat$P_value_g)))
all_stat$P_value_subt <- as.numeric(as.character(unlist(all_stat$P_value_subt)))
all_stat$Est_g <- as.numeric(as.character(unlist(all_stat$Est_g)))

all_stat$FDR_g <- p.adjust(all_stat$P_value_g, method='fdr')
all_stat$FDR_subt <- p.adjust(all_stat$P_value_subt, method='fdr')
all_stat <- all_stat[order(all_stat$FDR_g),]


write.table(all_stat,'out/DEGs_DEA_1_vs_allOthers.20250625.tsv',sep='\t',quote=F,row.names=F)
write.table(expr_5, 'out/Expr_across_TAL1_alt_groups.20250625.tsv',sep='\t',quote=F,row.names=F)



# now plot the top results:

all_stat <- read.table('out/DEGs_DEA_1_vs_allOthers.20250625.tsv',sep='\t',header=T)
expr_5 <- read.table('out/Expr_across_TAL1_alt_groups.20250625.tsv',sep='\t',header=T)

alt_order <- c('STIL-TAL1','TAL1 Upstream Enh','TAL1 Translocation','TAL1 Intron Enh','TAL1 Downstream Enh I','TAL1 Downstream Enh II')
cols_id <- brewer.pal(name='Dark2',n=6)[c(1,6,5,2,4,3)]
names(cols_id) <- alt_order


subt_cols <- brewer.pal(name='Set1',n=4)[c(2,4)]
names(subt_cols) <- c('TAL1 DP-like','TAL1 alpha/beta-like')

sel_genes <- all_stat$geneSymbol[all_stat$FDR_g <0.1 & !is.na(all_stat$FDR_g) & all_stat$FDR_g < all_stat$FDR_subt & abs(all_stat$Est_g)>0.2]
sel_genes <- c(sel_genes, 'LMO1')
to_plot <- expr_5[expr_5$geneSymbol %in% sel_genes,]
to_plot$Group <- factor(to_plot$Group, levels=alt_order)

# To get two scale colors, use solution from hereL https://stackoverflow.com/questions/75164501/how-to-have-two-scale-color-manuals-in-ggplot

pdf('plots/DEA_1_vs_allOthers_markers.v3.20250625.pdf',width=4,height=3.5,useDingbats=F)
for (gene in sel_genes){
    to_plot_2 <- to_plot[to_plot$geneSymbol==gene,]
    p <- ggplot(data=to_plot_2, aes(x=Group, y=Log2_expr))

    p <- p + geom_violin(width=0.7,aes(color=Group))

    p <- p + scale_color_manual(values=cols_id)
    p <- p + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width=0.8)
    p <- p + new_scale_color()
    p <- p + geom_jitter(position=position_jitter(0.2),size=0.5, alpha=0.5, shape=19,aes(color=Reviewed.subtype))

    p <- p + scale_color_manual(values=subt_cols)

    p <- p + theme_bw() + theme_nogrid() + labs(title="",x="",y="Gene expression")
    p <- p + theme(axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(colour="black",
    size=10)) 
    p <- p + theme_bw() + theme_nogrid() + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
    p <- p + theme(strip.text.x.top = element_text(angle = 90, hjust = 0.05, vjust = 0.5))+ggtitle(gene)
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=10,  angle=45, vjust = 1,hjust=1),
axis.text.y = element_text(colour="black", size=12),axis.ticks.x = element_blank(), legend.position='bottom')
print(p)
print(gene)
}
dev.off()


# make volcano plot
# add log2(Fold Change)
expr <- read.table('out/Expr_across_TAL1_alt_groups.20250625.tsv',sep='\t',header=T)

all_st <- NULL
for (gene in unique(expr$geneID)){
    expr_s <- expr[expr$geneID==gene,]
    gene_symbol <- unique(expr_s$geneSymbol)
    enh_1_e <- mean(expr_s$Log2_expr[expr_s$Group=='TAL1 Downstream Enh I'])
    other_e <- mean(expr_s$Log2_expr[expr_s$Group!='TAL1 Downstream Enh I'])
    f_ch <- enh_1_e - other_e
    st <- c(gene,gene_symbol,enh_1_e,other_e,f_ch)
    all_st <- rbind(all_st,st)
    }

all_st <- as.data.frame(all_st)
colnames(all_st) <- c('geneID','geneSymbol','Mean_enh_I_group','Mean_other_groups','Log2_fch')

all_stat <- read.table('out/DEGs_DEA_1_vs_allOthers.20250625.tsv',sep='\t',header=T)
all_stat <- merge(all_stat,all_st)
all_stat <- all_stat[order(all_stat$FDR_g),]
write.table(all_stat,'out/DEGs_DEA_1_vs_allOthers.withFch.20250625.tsv',sep='\t',quote=F,row.names=F)

# 2025-09-10: start from here from now on:
# Make line on FDR=0.2 (suggestive)
all_stat <- read.table('out/DEGs_DEA_1_vs_allOthers.withFch.20250625.tsv',sep='\t',header=T)
all_stat$FDR <- all_stat$FDR_g
all_stat$Direction <- ifelse(all_stat$Log2_fch > 0, 'UP', 'DOWN')

# get top markers
mark <- all_stat[all_stat$FDR_g < 0.05 & !is.na(all_stat$FDR_g),]
mark <- mark[order(-abs(mark$Log2_fch)),]
top_up <- mark$geneSymbol[mark$Log2_fch > 0][1:10]
top_down <- mark$geneSymbol[mark$Log2_fch < 0][1:10]

top_up <- c(top_up, "EPHA3", "LMO1","TLX2")

# remove duplicates (keep those with highest abs(Log2_fch))
all_stat <- all_stat[order(-abs(all_stat$Log2_fch)),]
all_stat <- all_stat[!duplicated(all_stat$geneSymbol),]


all_stat$Label <- ifelse(all_stat$geneSymbol %in% c(top_up, top_down), all_stat$geneSymbol, NA)

p <- ggplot(all_stat, aes(Log2_fch, -log10(FDR))) + geom_point_rast(alpha=0.4, colour="grey", data = all_stat) 
p <- p + geom_point_rast(data=all_stat[!is.na(all_stat$Label) ,],aes(color=Direction), alpha=0.95)
p <- p + theme_minimal() +geom_hline(yintercept=-log10(0.2), alpha=0.5)+geom_vline(xintercept=0, alpha=0.5)
p <- p + ylab("-log10(FDR)") + xlab("Log2(fold change)") 
p <- p + geom_text_repel(data=all_stat[!is.na(all_stat$Label),],
                         aes(y=-log10(FDR), x=Log2_fch, label=geneSymbol),alpha=1,colour = "black",segment.size=0.1,max.overlaps = Inf,
			 box.padding = 0.5,)
p <- p + theme(legend.position='right',legend.text=element_text(size=12),axis.title.x = element_text(size=12), 
               axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), 
               axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),
               panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
p <- p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

pdf('plots/DEA_1_vs_allOthers_markers.volcano_plot.FoldChange.v4.20251129.pdf',width=6,height=5,useDingbats=F)
p
dev.off()
