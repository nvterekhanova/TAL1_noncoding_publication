library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

# PDX:
l1 <- c('MYB', 43, 13)
l2 <- c('GATA3', 40, 11)
l3 <- c('H3K27ac', 14, 6)
l5 <- c('ATAC-seq', 168, 23)

to_plot <- as.data.frame(rbind(l1,l2,l3,l5))
to_plot$V2 <- as.numeric(as.character(unlist(to_plot$V2)))
to_plot$V3 <- as.numeric(as.character(unlist(to_plot$V3)))
colnames(to_plot) <- c('Data','Indel','WT')
to_plot$Sample <- 'SJTALL002048_X'
pdx_tab <- to_plot

# Jurkat:
l1 <- c('MYB', 47, 2)
l2 <- c('GATA3', 12, 1)
l3 <- c('RUNX1', 6, 0)

to_plot <- as.data.frame(rbind(l1,l2,l3))
to_plot$V2 <- as.numeric(as.character(unlist(to_plot$V2)))
to_plot$V3 <- as.numeric(as.character(unlist(to_plot$V3)))
colnames(to_plot) <- c('Data','Indel','WT')
to_plot$Sample <- 'Jurkat'
jurkat_tab <- to_plot

# combine two tables together:
to_plot <- rbind(pdx_tab,jurkat_tab)


to_plot$Percent_Indel <- to_plot$Indel*100/(to_plot$Indel+to_plot$WT)
to_plot$Percent_WT <- to_plot$WT*100/(to_plot$Indel+to_plot$WT)
to_plot <- to_plot[,c('Data','Percent_Indel','Percent_WT','Sample')]
colnames(to_plot)[2:3] <- c('Indel','WT')
to_plot <- melt(to_plot, id=c('Data','Sample'))
colnames(to_plot)[3:4] <- c('Reads','percent_reads')

to_plot$Data <- factor(to_plot$Data, levels=c('MYB','GATA3','RUNX1','ATAC-seq','H3K27ac'))
to_plot$Sample_2 <- ifelse(to_plot$Sample=='Jurkat','TAL1 Upstream Enh','TAL1 Downstream Enh I')
to_plot$Sample_2 <- factor(to_plot$Sample_2, levels=c('TAL1 Downstream Enh I','TAL1 Upstream Enh'))

# replace 0s with NAs:
to_plot$percent_reads[to_plot$percent_reads==0] <- NA

#color scale:
cols <- c('Indel'='#d73027', 'WT'='#4575b4')
p <- ggplot(to_plot, aes(x=Sample_2, y=percent_reads, fill=Reads))+ geom_bar(stat='identity',color='black', width=0.6, position='stack')
p <- p + scale_fill_manual(values=cols)
p <- p + theme_classic() + ylab('% of reads') + xlab('') 
p <- p + facet_grid(Data~., space='free', scales='free') #+ ylim(0,60)
p <- p + theme(axis.text = element_text(colour="black", size=12), legend.position='none')
p <- p + coord_flip() 
p <- p + theme(strip.text.y.right = element_text(angle = 0, hjust = 0.05, vjust = 0.5, size=10))
p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
pdf(paste("plots/Percent_reads_supporting_Insertion_and_WT_Jurkat.PDX.20260109.pdf",sep=""),width=5.5,height=3.5,useDingbats=F)
p
dev.off()




######################################################################################################
# Perform a two-sided binomial test with a null probability of success 0.5 in a Bernoulli experiment #
######################################################################################################

# PDX:
l1 <- c('MYB', 43, 13)
l2 <- c('GATA3', 40, 11)
l3 <- c('H3K27ac', 14, 6)
#l4 <- c('RUNX1', 1, 0)
l5 <- c('ATAC-seq', 168, 23)

to_plot <- as.data.frame(rbind(l1,l2,l3,l5))
to_plot$V2 <- as.numeric(as.character(unlist(to_plot$V2)))
to_plot$V3 <- as.numeric(as.character(unlist(to_plot$V3)))
colnames(to_plot) <- c('Data','Indel','WT')
to_plot$Sample <- 'SJTALL002048_X'
pdx_tab <- to_plot

# Jurkat:
l1 <- c('MYB', 47, 2)
l2 <- c('GATA3', 12, 1)
l3 <- c('RUNX1', 6, 0)

to_plot <- as.data.frame(rbind(l1,l2,l3))
to_plot$V2 <- as.numeric(as.character(unlist(to_plot$V2)))
to_plot$V3 <- as.numeric(as.character(unlist(to_plot$V3)))
colnames(to_plot) <- c('Data','Indel','WT')
to_plot$Sample <- 'Jurkat'
jurkat_tab <- to_plot

# combine two tables together:
to_plot <- rbind(pdx_tab,jurkat_tab)
to_plot$Coverage <- to_plot$Indel+to_plot$WT

all_st <- NULL
for (sample in unique(to_plot$Sample)){
    tab <- to_plot[to_plot$Sample==sample,]
    all_st_sample <- NULL
    for (i in 1:nrow(tab)){
        str <- tab[i,]
        st <- binom.test(str$Indel, str$Coverage, p = 0.5, alternative='two.sided')
        str$p_val <- st$p.value
        all_st_sample <- rbind(all_st_sample, str)
        }
    all_st_sample$FDR <- p.adjust(all_st_sample$p_val, metho='fdr')
    all_st <- rbind(all_st, all_st_sample)
}

write.table(all_st, 'out/N_reads_with_Indel_and_WT_SJTALL002048_Jurkat.binomial_test_res.20260109.tsv',sep='\t',quote=F)

# now also compare MYB counts and GATA3 counts between PDX and Jurkat samples:
t1 <- to_plot[to_plot$Data=='MYB',]
t1 <- t1[,2:3]
st1 <- fisher.test(t1)

t1 <- to_plot[to_plot$Data=='GATA3',]
t1 <- t1[,2:3]
st2 <- fisher.test(t1)


