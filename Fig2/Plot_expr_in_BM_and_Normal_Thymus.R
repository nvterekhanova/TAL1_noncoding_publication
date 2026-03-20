library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)

obj <- readRDS(paste('/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data/Public_DATA/scRNAseq_Cordes_2022/',
'GSE195812_HCACD34BM_MantonBM_CordesThymus_integrated.rds',sep=''))

table(obj$'Final Annotations')
table(obj$'Global Annotations')

sel_cell_t <- c('HSC','MPP','LMPP','CLP','DN1','DN2','DN3','ISP','DP_CD3min','DP_CD3plus','CD4','CD8')

Idents(obj) <- obj$'Global Annotations'
obj_s <- subset(obj, idents = sel_cell_t)
Idents(obj_s) <- factor(Idents(obj_s), levels=sel_cell_t)

# obj contain only 1 assay, RNA
markers <- c('TAL1','CD34','CD7','CD1B','RAG1','CD4','CD8A')

#use scale.by="size", based on this: https://bioinformatics.stackexchange.com/questions/10738/how-do-i-increase-the-minimum-dot-size-in-seurats-dotplot-function
plot <- DotPlot(obj_s, assay='RNA', features=markers, scale.by = "size") +
scale_colour_gradient2(low = "#0d0887", mid = "#cc4778", high = "#f0f921") + coord_flip() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf('plots/Markers_v1.BM_and_Thymus.20250403.pdf', width=8, height=4)
print(plot)
dev.off()

# 2025-04-10, make UMAP plot:
Idents(obj) <- obj$'Global Annotations'

cols <- c(brewer.pal(name='Paired',n=12),brewer.pal(name='Dark2',n=8))
names(cols) <-  unique(Idents(obj))
p <- DimPlot(object = obj, reduction = "UMAP",label=TRUE,label.size=6, pt.size=0.2) + scale_color_manual(values=cols)

pdf("plots/MantonBM_CordesThymus_integrated.UMAP.pdf",height=7,width=8)
print(p)
dev.off()

p <- DimPlot(object = obj, reduction = "UMAP3D",label=TRUE,label.size=6, pt.size=0.2) + scale_color_manual(values=cols)

pdf("plots/MantonBM_CordesThymus_integrated.UMAP3D.pdf",height=7,width=8)
print(p)
dev.off()

sel_cell_t <- c('HSC','MPP','LMPP','CLP','DN1','DN2','DN3','ISP','DP_CD3min','DP_CD3plus','CD4','CD8')

Idents(obj) <- obj$'Global Annotations'
obj_s <- subset(obj, idents = sel_cell_t)
Idents(obj_s) <- factor(Idents(obj_s), levels=sel_cell_t)

# obj contain only 1 assay, RNA
markers <- c('TAL1','MYB','GATA3','RUNX1','CD34','CD7','CD1B','RAG1','CD4','CD8A')

#use scale.by="size",
#based on this: https://bioinformatics.stackexchange.com/questions/10738/how-do-i-increase-the-minimum-dot-size-in-seurats-dotplot-function
plot <- DotPlot(obj_s, assay='RNA', features=markers, scale.by = "size") +
scale_colour_gradient2(low = "#0d0887", mid = "#cc4778", high = "#f0f921") + coord_flip() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf('plots/Markers_v1.BM_and_Thymus.more_markers.20250410.pdf', width=8, height=4)
print(plot)
dev.off()

cols <- brewer.pal(name='Paired',n=12)
names(cols) <-  unique(Idents(obj_s))
p <- DimPlot(object = obj_s, reduction = "UMAP",label=TRUE,label.size=6, pt.size=0.2) + scale_color_manual(values=cols)

pdf("plots/MantonBM_CordesThymus_integrated.sel_cell_types.UMAP.pdf",height=7,width=8)
print(p)
dev.off()

# save average expression per cell group:

obj <- readRDS(paste('/research_jude/rgs01_jude/groups/zhanggrp/projects/TALL_noncoding/common/nterekha/Data/Public_DATA/scRNAseq_Cordes_2022/',
'GSE195812_HCACD34BM_MantonBM_CordesThymus_integrated.rds',sep=''))

Idents(obj) <- obj$'Global Annotations'

var_genes <- rownames(obj)

expr_aver <- AverageExpression(obj, assays = 'RNA', slot ='data',features=var_genes)
tab_expr <- expr_aver$RNA


file2write <- paste0("out/AverageExpressionRNA_MantonBM_CordesThymus_integrated_SlotData_AllGenes.20250410.tsv")
write.table(tab_expr, file = file2write, quote = F, sep = "\t", row.names = T)
