library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(reshape2)
library(tidyverse)
library(data.table)

#Read gene list annotation from HGNC
genes_list=read.csv('hgnc-gene_list.txt',sep='\t')
genes_prot_c=genes_list[genes_list$Locus.type=='Gene with protein product',]


degs=read.table('out/Pearson_corr_PT_with_Regulons_expr.min500Cells.20240508.tsv',sep='\t',header=T)
degs=degs[degs$FDR<0.05,]
degs=degs[order(degs$FDR),]
#degs=degs[degs$Feature %in% genes_prot_c$Symbol,]



tab=read.table('out/Regulons_scores_data.CD8.20240508.tsv',sep='\t',header=T)
tab=reshape2::dcast(data=tab,Feature~Barcode,value.var='Score')
rownames(tab)=tab$Feature
tab=tab[,-1]

pt_tab=read.table('../20240508_Correlation_with_PT_ATAC_RNA_shared/out/All_PTs_monocle_objs.10comp.CD8.After_reOrdering.state16_root.20240313.tsv',sep='\t',header=T)

#length(intersect(colnames(tab),pt_tab$sample_name))==ncol(tab)

pt_tab=pt_tab[order(pt_tab$Pseudotime),]

tab=tab[,pt_tab$sample_name]

#Keep only top N degs:
tab_s=tab[degs$Feature[1:50],]


#Make a heatamap:
cell_type_cols=c(brewer.pal(n=8, name='Dark2'))
names(cell_type_cols)=unique(pt_tab$cell_type)

col_fun_pt=colorRamp2(c(0, 20, 40), c('#ffffcc','#fd8d3c','#800026'))

col_fun_state= c(brewer.pal(name='Dark2',n=6),brewer.pal(name='Paired',n=11))
names(col_fun_state)=unique(pt_tab$State)


column_ha=HeatmapAnnotation(Cell_type=pt_tab$cell_type, Pseudotime=pt_tab$Pseudotime, State=pt_tab$State, col=list(Cell_type=cell_type_cols, Pseudotime=col_fun_pt,State=col_fun_state))


color_scale=colorRamp2(c(-1, 0, 1), c('#0d0887','#cc4778','#f0f921'))

x=Heatmap(t(scale(t(tab_s))),bottom_annotation=column_ha,col= color_scale,show_row_names = T,show_column_names = FALSE,show_row_dend=T,show_column_dend=FALSE,name='Gene expression', cluster_rows=T, use_raster=T,cluster_columns=F)


pdf("plots/Top_50_selected_Regulons.20240508.pdf",width=10,height=9,useDingbats=FALSE)
print(x)
dev.off()


#########################
#Now try more regulons###
#########################

N_genes=100

tab=read.table('out/Regulons_scores_data.CD8.20240508.tsv',sep='\t',header=T)
tab=reshape2::dcast(data=tab,Feature~Barcode,value.var='Score')
rownames(tab)=tab$Feature
tab=tab[,-1]

pt_tab=read.table('../20240508_Correlation_with_PT_ATAC_RNA_shared/out/All_PTs_monocle_objs.10comp.CD8.After_reOrdering.state16_root.20240313.tsv',sep='\t',header=T)

#length(intersect(colnames(tab),pt_tab$sample_name))==ncol(tab)

pt_tab=pt_tab[order(pt_tab$Pseudotime),]

tab=tab[,pt_tab$sample_name]

#Keep only top N degs:
tab_s=tab[degs$Feature[1:N_genes],]



cell_type_cols=c(brewer.pal(n=8, name='Dark2'))
names(cell_type_cols)=unique(pt_tab$cell_type)

col_fun_pt=colorRamp2(c(0, 20, 40), c('#ffffcc','#fd8d3c','#800026'))

col_fun_state= c(brewer.pal(name='Dark2',n=6),brewer.pal(name='Paired',n=11))
names(col_fun_state)=unique(pt_tab$State)

column_ha=HeatmapAnnotation(Cell_type=pt_tab$cell_type, Pseudotime=pt_tab$Pseudotime, State=pt_tab$State, col=list(Cell_type=cell_type_cols, Pseudotime=col_fun_pt,State=col_fun_state))

degs_to_label=rownames(tab_s) %in% degs$Feature[1:50]

row_ha=rowAnnotation(link = anno_mark(at = which(degs_to_label), labels = rownames(tab_s[rownames(tab_s) %in% degs$Feature[1:50]]),labels_gp = gpar(fontsize = 12),side='right'))

color_scale_rna=colorRamp2(c(-1, 0, 1), c('#440154','#21918c','#fde725'))
color_scale_atac=colorRamp2(c(-1, 0, 1), c('#1A194D', '#B2519F', '#FBF05E'))

color_scale=colorRamp2(c(-1, 0, 1), c('#0d0887','#cc4778','#f0f921'))
x=Heatmap(t(scale(t(tab_s))),bottom_annotation=column_ha,right_annotation=row_ha,col= color_scale,show_row_names = F,show_column_names = FALSE,show_row_dend=T,show_column_dend=FALSE,name='Gene expression', cluster_rows=T, use_raster=T,cluster_columns=F)


pdf(paste("plots/Top_",N_genes,"_selected_Regulons.20240508.pdf",sep=''),width=10,height=10,useDingbats=FALSE)
print(x)
dev.off()






