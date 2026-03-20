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

data_dir='/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data'
theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

tab <- read.table('../20250620_Cohort_overview/out/TAL1_alterations_Allsamples_with_TAL1_WT.annotated.20250623.tsv',sep='\t',header=T)
tab_s <- tab[,c('USI','SJ_ID','Group','Reviewed.subtype')]

# Read table with alterations annotations
alt_info <- read_excel(paste(data_dir,'/RNAseq/Polonen_2024_samples/41586_2024_7807_MOESM4_ESM.xlsx',sep=''),sheet=10)
alt_info <- as.data.frame(alt_info)
alt_info <- alt_info[alt_info$USI %in% tab_s$USI,]

alt_info_s <- alt_info[!duplicated(alt_info[,c('USI','Targeted.gene')]),]
alt_info_s <- alt_info_s[,c('USI','Targeted.gene','Alteration','variantID')]

# add PARSNX TAL1 SNV/Indel manually:
alt_info_s <- rbind(alt_info_s,c('PARSNX','TAL1','SNV/Indel','PARSNX__Indel__chr1_47203728'))

alt_info_s2 <- merge(alt_info_s,tab_s)

alt_info_s2$Count <- 1
stat <- aggregate(alt_info_s2$Count, by=list(alt_info_s2$Reviewed.subtype,alt_info_s2$Targeted.gene),FUN='sum')
colnames(stat) <- c('Reviewed.subtype','Targeted.gene','Count')

tab_s$Count <- 1
total <- aggregate(tab_s$Count, by=list(tab_s$Reviewed.subtype), FUN='sum')
colnames(total) <- c('Reviewed.subtype','Total')

res <- merge(stat,total,all.x=T)
res$Fraction <- res$Count/res$Total

genes_sel <- unique(res$Targeted.gene[res$Fraction>0.5])
genes_sel <- genes_sel[!(genes_sel %in% c('9q','9q34.11'))]

stat_2 <- stat[stat$Targeted.gene %in% genes_sel,]
stat_2 <- reshape2::dcast(data=stat_2,Reviewed.subtype~Targeted.gene,value.var='Count')
row.names(stat_2) <- stat_2[,1]
stat_2 <- stat_2[,-1]
stat_2[is.na(stat_2)] <- 0
stat_3 <- reshape2::melt(as.matrix(stat_2))
colnames(stat_3) <- c('Reviewed.subtype','Targeted.gene','Count')

res_2 <- merge(stat_3,total,all.x=T)
res_2$Fraction <- res_2$Count/res_2$Total

# run Fisher.test

tab <- res_2
subts <- unique(tab$Reviewed.subtype)
res_st <- NULL
for (subt in subts){
    all_st <- NULL
    for (gene in unique(tab$Targeted.gene)){
    	tab_g <- tab[tab$Targeted.gene==gene,]
	tab_s <- tab_g[tab_g$Reviewed.subtype==subt,]
        tab_other <- tab_g[tab_g$Reviewed.subtype!=subt,]
        fract <- tab_s$Fraction
        subt_mut <- tab_s$Count
        subt_no_mut <- tab_s$Total - tab_s$Count

	other_mut <- sum(tab_other$Count)
        other_no_mut <- sum(tab_other$Total) - sum(tab_other$Count)
        other_c_fr <- other_mut/(other_mut+other_no_mut)
        mat <- matrix(c(subt_mut,subt_no_mut,other_mut,other_no_mut),nrow=2)
        test <- fisher.test(mat,alternative='greater')
        st <- c(gene, subt,test$p.value,fract,other_c_fr)
        all_st <- rbind(all_st,st)
}
all_st <- as.data.frame(all_st)
colnames(all_st) <- c('Targeted.gene','Reviewed.subtype','P_value','Fraction','Other_subt_fraction')
all_st$P_value <- as.numeric(as.character(unlist(all_st$P_value)))
all_st$FDR <- p.adjust(all_st$P_value,method='fdr')
res_st <- rbind(res_st,all_st)
}
res_st$Fraction <- as.numeric(as.character(unlist(res_st$Fraction)))
write.table(res_st, 'out/Fisher_test_results_for_gene_mut_enr_acrossSubtypes.20251021.tsv',sep='\t',quote=F,row.names=F)



# Now do the analysis for samples separated by Group:
alt_info_s3 <- alt_info_s2[!is.na(alt_info_s2$Group) & alt_info_s2$Group!='TAL1 WT',]
stat <- aggregate(alt_info_s3$Count, by=list(alt_info_s3$Group,alt_info_s3$Targeted.gene),FUN='sum')
colnames(stat) <- c('Group','Targeted.gene','Count')

tab <- read.table('../20250620_Cohort_overview/out/TAL1_alterations_Allsamples_with_TAL1_WT.annotated.20250623.tsv',sep='\t',header=T)
tab_s <- tab[,c('USI','SJ_ID','Group','Reviewed.subtype')]
tab_s$Count <- 1
total <- aggregate(tab_s$Count, by=list(tab_s$Group), FUN='sum')
colnames(total) <- c('Group','Total')

res <- merge(stat,total,all.x=T)
res$Fraction <- res$Count/res$Total

genes_sel <- unique(res$Targeted.gene[res$Fraction>0.1])
genes_sel <- genes_sel[!(genes_sel %in% c('9p','6q','17q','12p','8p','8q','8p12','9q'))]

stat_2 <- stat[stat$Targeted.gene %in% genes_sel,]
stat_2 <- reshape2::dcast(data=stat_2,Group~Targeted.gene,value.var='Count')
row.names(stat_2) <- stat_2[,1]
stat_2 <- stat_2[,-1]
stat_2[is.na(stat_2)] <- 0
stat_3 <- reshape2::melt(as.matrix(stat_2))
colnames(stat_3) <- c('Group','Targeted.gene','Count')

res_2 <- merge(stat_3,total,all.x=T)
res_2$Fraction <- res_2$Count/res_2$Total

# run Fisher.test

# substitute "Reviewed.subtype" to "Group"
tab <- res_2
subts <- unique(tab$Group)
res_st <- NULL
for (subt in subts){
    all_st <- NULL
    for (gene in unique(tab$Targeted.gene)){
    	tab_g <- tab[tab$Targeted.gene==gene,]
	tab_s <- tab_g[tab_g$Group==subt,]
        tab_other <- tab_g[tab_g$Group!=subt,]
        fract <- tab_s$Fraction
        subt_mut <- tab_s$Count
        subt_no_mut <- tab_s$Total - tab_s$Count

	other_mut <- sum(tab_other$Count)
        other_no_mut <- sum(tab_other$Total) - sum(tab_other$Count)
        other_c_fr <- other_mut/(other_mut+other_no_mut)
        mat <- matrix(c(subt_mut,subt_no_mut,other_mut,other_no_mut),nrow=2)
        test <- fisher.test(mat,alternative='greater')
        st <- c(gene, subt,test$p.value,fract,other_c_fr)
        all_st <- rbind(all_st,st)
}
all_st <- as.data.frame(all_st)
colnames(all_st) <- c('Targeted.gene','Group','P_value','Fraction','Other_subt_fraction')
all_st$P_value <- as.numeric(as.character(unlist(all_st$P_value)))
all_st$FDR <- p.adjust(all_st$P_value,method='fdr')
res_st <- rbind(res_st,all_st)
}
res_st$Fraction <- as.numeric(as.character(unlist(res_st$Fraction)))
write.table(res_st, 'out/Fisher_test_results_for_gene_mut_enr_acrossGroups.20251021.tsv',sep='\t',quote=F,row.names=F)
