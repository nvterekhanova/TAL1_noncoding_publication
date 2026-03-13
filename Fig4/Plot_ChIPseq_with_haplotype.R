library(ggplot2)
library(RColorBrewer)
library(reshape)
library(reshape2)

tab <- read.table('out/RNA_H3K27ac_H3K27me3_assays_SNPs_with_haplotype_and_ASE_annot.20251214.tsv',sep='\t',header=T)

tab_2 <- tab[tab$Sample=='SJTALL002048' & tab$Assay %in% c('H3K27ac','H3K27me3','RNA_PDX'),]

#change to Haplotype_1 VAF everywhere:
tab_2$Assay_vaf <- ifelse(tab_2$haplotype==1,tab_2$Assay_vaf,1-tab_2$Assay_vaf)

tab_h <- tab_2[!duplicated(tab_2[,c('position','haplotype')]),]
tab_h2 <- reshape2::dcast(data=tab_h, chr+position+REF+ALT~haplotype,value.var='Sample')

# for variant with denoted haplotype assign 1.2, and for other haplotype assign -0.2 (for plotting):
colnames(tab_h2)[5:6] <- c('Hap_1','Hap_2')
tab_h2$Hap_1 <- ifelse(tab_h2$Hap_1=='SJTALL002048' & !is.na(tab_h2$Hap_1),1.2,-0.2)
tab_h2$Hap_2 <- ifelse(tab_h2$Hap_2=='SJTALL002048' & !is.na(tab_h2$Hap_2),1.2,-0.2)
tab_h2$Assay_depth <- 30
tab_h2$Sample <- 'SJTALL002048'
tab_h3 <- reshape2::melt(tab_h2,id.vars=c('chr','position','REF','ALT','Sample','Assay_depth'))
colnames(tab_h3)[7:8] <- c('Assay','Assay_vaf')
tab_h3$Data <- 'Variant'

# change Hap_1 to 'ALT', and Hap_2 to REF (as we want to show where Hap_1 correspond to ALT or REF and it was before coded in their frequencies)
tab_h3$Assay <- ifelse(tab_h3$Assay=='Hap_1','ALT','REF')
tab_h3$FDR <- 1

tab_3 <- tab_2[,c('chr','position','REF','ALT','Sample','Assay_depth','Assay','Assay_vaf','FDR')]
tab_3$Data <- 'Assay'
tab_3 <- tab_3[,colnames(tab_h3)]

both <- rbind(tab_3,tab_h3)

# calculate median for the 5 clusters of significant ASEs:
tab_ase <- tab_3[tab_3$FDR < 0.05,]
cl_h3k27me3 <- median(tab_ase$Assay_vaf[tab_ase$Assay=='H3K27me3'])
cl_rna <- median(tab_ase$Assay_vaf[tab_ase$Assay=='RNA_PDX'])
# 6 SNPs (2 of those are at adjacent nts, so not seen at the plot)
cl_h3k27ac_1 <- median(tab_ase$Assay_vaf[tab_ase$Assay=='H3K27ac' & tab_ase$position < 47673000])
# 4 SNPs
cl_h3k27ac_2 <- median(tab_ase$Assay_vaf[tab_ase$Assay=='H3K27ac' & tab_ase$position > 47673000 & tab_ase$position < 47682000])
# 3 SNPs
cl_h3k27ac_3 <- median(tab_ase$Assay_vaf[tab_ase$Assay=='H3K27ac' & tab_ase$position > 47682000])



# plotting
to_plot <- both
col_intercept <- 'black'
p <- ggplot(to_plot, aes(x=position,y=Assay_vaf,fill=Assay)) + geom_hline(yintercept=0.5, color='grey', linetype = "dashed", linewidth = 0.5)
p <- p + geom_hline(yintercept=cl_rna, color=col_intercept, linewidth = 0.3) + geom_hline(yintercept=cl_h3k27me3, color=col_intercept, linewidth = 0.3)
p <- p + geom_hline(yintercept=cl_h3k27ac_1, color=col_intercept, linewidth = 0.3)
p <- p + geom_hline(yintercept=cl_h3k27ac_2, color=col_intercept, linewidth = 0.3)
p <- p + geom_hline(yintercept=cl_h3k27ac_3, color=col_intercept, linewidth = 0.3)
p <- p + geom_point(alpha=0.8, aes(size=Assay_depth,pch=Data,stroke=ifelse(FDR<0.05,1,ifelse(FDR>=0.05 & FDR<0.15,0.5,NA))),color='black') 
p <- p + theme_bw() + theme_classic()
p <- p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))
p <- p + scale_fill_manual(values=c('H3K27ac'='#E41A1C','H3K27me3'='#377EB8','RNA_PDX'='#984ea3','ALT'='red','REF'='blue'))
p <- p + scale_shape_manual(values=c('Assay'=21,'Variant'=22))
p <- p + scale_size_continuous(name='Coverage',range=c(1,8)) + scale_y_continuous(breaks=seq(0,1,by=0.2))
p <- p + theme(axis.text.x = element_text(colour="black"),
axis.text.y = element_text(colour="black", size=12))
p <- p + labs(x = "Genomic coordinate", y = "Hap1 VAF in assay")

pdf(paste("plots/Summary_ChIP_vars_withHapl_info.min5depth.dotplot.specific_window.v3.20251219.pdf",sep=''),width=14, height=4.3,useDingbats=FALSE)
p
dev.off()

