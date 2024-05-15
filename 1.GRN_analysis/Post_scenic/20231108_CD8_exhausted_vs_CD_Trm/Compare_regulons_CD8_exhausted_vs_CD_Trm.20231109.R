#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(Matrix)
library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(GSEABase)
library(AUCell)


##############################################
#Do comparison for groups of interest#########
##############################################

library(data.table)

wdir='/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/Analysis/1.GRN_analysis/Run.v.20230912'

mat=fread(paste(wdir, '/Post_scenic/20230913_Filter_regulons/Update_regByFreq/out/UpdatedByFreq.0.8_regulons_CellsAUC.20230914.tsv',sep=''))

mat=as.data.frame(mat)
rownames(mat)=mat[,1]
mat=mat[,-1]

res=mat

annot=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/DATA/data_freeze_v1/',
'PanImmune_integrated_allRNA_by_chemistry_v8.0_df_T_cells_NK_no_doublets_metadata_09112023_res1.4.tsv',sep=''),sep='\t',header=T)
annot_s=annot[,c('X','cell_type_v8.4_rna')]
colnames(annot_s)[1]='Barcode'

info=read.table(paste(wdir, '/Object/T_cell_PancanObj_1500cellsPer_seurat_cluster.20230912.tsv',sep=''), sep='\t', header=T)
info=merge(info,annot_s)
info=info[,c(1:2,4)]

rownames(info)=info$Barcode
info=info[colnames(res),]
info$Cohort=gsub('(.*)_(.*)_(.*)','\\1',info$Barcode)
info$Cohort=gsub('(.*)_(.*)','\\1',info$Cohort)


reg_annot=read.table(paste(wdir, '/Post_scenic/20230913_Filter_regulons/Update_regByFreq/out/Regulons_new_annot.20230914.tsv',
sep=''),sep='\t',header=T)
sel_regs=reg_annot$Regulon[reg_annot$Genes_N>=20]
#cell_types=unique(info$cell_type_v8.4_rna)

cell_t1='CD8 T-cells exhausted GZMK ENTPD1'
cell_t2='CD8 Trm exhausted ENTPD1'
final_wilcoxon_stat=NULL
#for (cell_t1 in cell_types){
    print(cell_t1)
    res_1=res[,colnames(res) %in% rownames(info)[info$cell_type_v8.4_rna==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(info)[info$cell_type_v8.4_rna==cell_t2]]
    all_wilcoxon_stat=NULL
    for (motif in 1:nrow(res)){
    	   tf=rownames(res)[motif]
	   if (tf %in% sel_regs){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
	}
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','Regulon')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
#}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat$cell_t2=cell_t2
write.table(final_wilcoxon_stat,'out/Regulons_AUC_difference.CD8_exhausted_vs_CD_Trm.20231108.tsv',sep='\t',quote=F,row.names=F)


################################################
####Also select 200 cells per group and save:###
################################################
info$ID=info$cell_type_v8.4_rna

n_cells=200
all_cells_s=NULL
for (id in unique(info$ID)){
    info_1=info[info$ID==id,]
    cells_s=sample(rownames(info_1),min(n_cells,nrow(info_1)))
    all_cells_s=c(all_cells_s,cells_s)
}


###Also save AUC-scores for the same cells:
res_s=res[,all_cells_s]
info_sel=info[info$Barcode %in% colnames(res_s),]

write.table(res_s, 'out/AUC_Pancan_Immune_200_RandomCellsPerCellType.20230915.tsv',sep='\t',quote=F)
write.table(info_sel, 'out/Annot_cells_AUC_Pancan_Immune_200_RandomCellsPerCellType.20230915.tsv',sep='\t',quote=F)


##########################
###2023-09-27#############
##########################
# Now also do comparison for each cluster vs others

##############################################
#Do comparison for each cluster vs others)####
##############################################

library(data.table)
mat=fread('out/UpdatedByFreq.0.8_regulons_CellsAUC.20230914.tsv')
mat=as.data.frame(mat)
rownames(mat)=mat[,1]
mat=mat[,-1]

res=mat

annot=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/DATA/data_freeze_v1/',
'PanImmune_integrated_allRNA_by_chemistry_v8.0_df_T_cells_NK_no_doublets_metadata_09112023_res1.4.tsv',sep=''),sep='\t',header=T)
annot_s=annot[,c('X','cell_type_v8.4_rna')]
colnames(annot_s)[1]='Barcode'

info=read.table('../../../Object/T_cell_PancanObj_1500cellsPer_seurat_cluster.20230912.tsv',sep='\t',header=T)
info=merge(info,annot_s)
info=info[,c(1:2,4)]

rownames(info)=info$Barcode
info=info[colnames(res),]
info$Cohort=gsub('(.*)_(.*)_(.*)','\\1',info$Barcode)
info$Cohort=gsub('(.*)_(.*)','\\1',info$Cohort)


reg_annot=read.table('out/Regulons_new_annot.20230914.tsv',sep='\t',header=T)
sel_regs=reg_annot$Regulon[reg_annot$Genes_N>=20]

cell_types=unique(info$Seurat_clusters)

final_wilcoxon_stat=NULL
for (cell_t1 in cell_types){
    print(cell_t1)
    res_1=res[,colnames(res) %in% rownames(info)[info$Seurat_clusters==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(info)[info$Seurat_clusters!=cell_t1]]
    all_wilcoxon_stat=NULL
    for (motif in 1:nrow(res)){
    	   tf=rownames(res)[motif]
	   if (tf %in% sel_regs){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
	}
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','Regulon')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat$cell_t2=cell_t2
write.table(final_wilcoxon_stat,'out/Regulons_AUC_difference.Each_cluster_vs_others.20230927.tsv',sep='\t',quote=F,row.names=F)

