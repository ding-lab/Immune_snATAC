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


degs=read.table('out/Pearson_corr_PT_with_RNA_expr.min500Cells.20240430.tsv',sep='\t',header=T)
degs=degs[degs$FDR<0.05,]
degs=degs[order(degs$FDR),]
degs=degs[degs$Feature %in% genes_prot_c$Symbol,]

dacrs=read.table('out/Pearson_corr_PT_with_ATAC_expr.min50Cells.forPeaks_assoc_wDEGs.20240508.tsv',sep='\t',header=T)
dacrs=dacrs[dacrs$FDR<0.05,]
dacrs=dacrs[order(dacrs$FDR),]

#Now retain the peaks with the same direction as RNA expression for the same gene:
degs_s=degs[,c('Feature','Estimate','FDR')]
colnames(degs_s)=c('Gene','RNA_estimate','RNA_FDR')

dacrs_2=merge(dacrs,degs_s,all.x=T)

dacrs_3=dacrs_2[(dacrs_2$Estimate>0 & dacrs_2$RNA_estimate>0) | (dacrs_2$Estimate<0 & dacrs_2$RNA_estimate<0),]

dacrs_3=dacrs_3[!duplicated(dacrs_3$Gene),]

#Also read selected Genes for RNA heatmap:
rna_rds=readRDS('for_plot/300_Genes_order_inRNA_heatmap.rds')

dacrs_3=dacrs_3[dacrs_3$Gene %in% rna_rds,]
rownames(dacrs_3)=dacrs_3$Gene
dacrs_3=dacrs_3[rna_rds,]

tab=read.table('out/ATAC_CD8_correlated_matrix.associated_wDEGs.20240508.tsv',sep='\t',header=T)
colnames(tab)<-gsub('\\.','\\-',colnames(tab))

pt_tab=read.table('out/All_PTs_monocle_objs.10comp.CD8.After_reOrdering.state16_root.20240313.tsv',sep='\t',header=T)

#length(intersect(colnames(tab),pt_tab$sample_name))==ncol(tab)

pt_tab=pt_tab[order(pt_tab$Pseudotime),]

tab=tab[,pt_tab$sample_name]

#Keep only top N degs:
tab_s=tab[dacrs_3$Feature[1:50],]


#Make a heatamap:
cell_type_cols=c(brewer.pal(n=8, name='Dark2'))
names(cell_type_cols)=unique(pt_tab$cell_type)

col_fun_pt=colorRamp2(c(0, 20, 40), c('#ffffcc','#fd8d3c','#800026'))

col_fun_state= c(brewer.pal(name='Dark2',n=6),brewer.pal(name='Paired',n=11))
names(col_fun_state)=unique(pt_tab$State)


column_ha=HeatmapAnnotation(Cell_type=pt_tab$cell_type, Pseudotime=pt_tab$Pseudotime, State=pt_tab$State, col=list(Cell_type=cell_type_cols, Pseudotime=col_fun_pt,State=col_fun_state))
row_ha=rowAnnotation(ID=anno_text(dacrs_3$Gene[1:50],gp = gpar(fontsize = 10)))


#color_scale=colorRamp2(c(-1, 0, 1), c('#1A194D', '#B2519F', '#FBF05E'))
#color_scale=colorRamp2(c(-1, 0, 1), c('#440154','#21918c','#fde725'))

color_scale_rna=colorRamp2(c(-1, 0, 1), c('#2C2B52', '#3AB69F', '#F1EB47'))
color_scale_atac=colorRamp2(c(-1, 0, 1), c('#0d0887','#cc4778','#f0f921'))
x=Heatmap(t(scale(t(tab_s))),bottom_annotation=column_ha, right_annotation=row_ha, col= color_scale_atac,show_row_names = FALSE, show_column_names = FALSE,show_row_dend=T,show_column_dend=FALSE,name='Peak accessibility', cluster_rows=T, use_raster=T,cluster_columns=F)


pdf("plots/Top_50_selected_ATAC.min500Cells.shared_withRNA.20240508.pdf",width=10,height=9,useDingbats=FALSE)
print(x)
dev.off()


######################
#Now try more DACRs###
######################

#N_genes=nrow(dacrs_3)
N_genes=300

tab=read.table('out/ATAC_CD8_correlated_matrix.associated_wDEGs.20240508.tsv',sep='\t',header=T)
colnames(tab)<-gsub('\\.','\\-',colnames(tab))

pt_tab=read.table('out/All_PTs_monocle_objs.10comp.CD8.After_reOrdering.state16_root.20240313.tsv',sep='\t',header=T)

#length(intersect(colnames(tab),pt_tab$sample_name))==ncol(tab)

pt_tab=pt_tab[order(pt_tab$Pseudotime),]

tab=tab[,pt_tab$sample_name]

#Keep only top N degs:
tab_s=tab[dacrs_3$Feature[1:N_genes],]

top_degs=dacrs_3[order(dacrs_3$RNA_FDR),]

cell_type_cols=c(brewer.pal(n=8, name='Dark2'))
names(cell_type_cols)=unique(pt_tab$cell_type)

col_fun_pt=colorRamp2(c(0, 20, 40), c('#ffffcc','#fd8d3c','#800026'))

col_fun_state= c(brewer.pal(name='Dark2',n=6),brewer.pal(name='Paired',n=11))
names(col_fun_state)=unique(pt_tab$State)

column_ha=HeatmapAnnotation(Cell_type=pt_tab$cell_type, Pseudotime=pt_tab$Pseudotime, State=pt_tab$State, col=list(Cell_type=cell_type_cols, Pseudotime=col_fun_pt,State=col_fun_state))

dacrs_to_label=rownames(tab_s) %in% top_degs$Feature[1:50]

row_ha=rowAnnotation(link = anno_mark(at = which(dacrs_to_label), labels = dacrs_3$Gene[rownames(tab_s) %in% top_degs$Feature[1:50]],labels_gp = gpar(fontsize = 12),side='right'))

color_scale_rna=colorRamp2(c(-1, 0, 1), c('#2C2B52', '#3AB69F', '#F1EB47'))
color_scale_atac=colorRamp2(c(-1, 0, 1), c('#0d0887','#cc4778','#f0f921'))
x=Heatmap(t(scale(t(tab_s))),bottom_annotation=column_ha,right_annotation=row_ha,col= color_scale_atac,show_row_names = FALSE, show_column_names = FALSE,show_row_dend=T,show_column_dend=FALSE,name='Gene expression', cluster_rows=F, use_raster=T,cluster_columns=F)


pdf(paste("plots/Top_",N_genes,"_selected_ATACs.min500Cells.shared_withRNA.20240508.pdf",sep=''),width=10,height=10,useDingbats=FALSE)
print(x)
dev.off()

#The heatmap with results looks quite sparse, maybe need imputation here...






