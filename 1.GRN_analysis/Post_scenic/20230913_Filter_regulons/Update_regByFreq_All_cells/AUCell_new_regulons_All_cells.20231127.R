#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(Matrix)
library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(GSEABase)
library(AUCell)



###https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/importing_pySCENIC.html
###Import from the loom file

#We don't need filter genes, as they were filtered before inferring the regulons

wdir='/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/Analysis/1.GRN_analysis/Run.v.20230912'
loom=open_loom(paste(wdir,'/Object/data/T_cell_PancanObj_All_cells.loom',sep=''),mode='r')

###Follow the tutorial here http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
#Check AUCell here:
#https://www.bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html

###Make a new regulons-array using updated regulons:
tfs=read.table('../TF_frequency.10Iter.20230914.tsv',sep='\t',header=T)
tfs_s=tfs[tfs$Freq>=8,]
targ=read.table('../Target_frequency_byTF.10Iter.20230914.tsv',sep='\t',header=T)
targ_s=targ[targ$Freq>=8,]

regulons_new=vector(mode='list',length=nrow(tfs_s))
for (i in 1: length(tfs_s$all_tfs)){
    tf=tfs_s$all_tfs[i]
    genes=targ_s$Gene[targ_s$TF==tf]
    regulons_new[[i]]=genes
    names(regulons_new)[[i]]=tf
}
exprMatrix=get_dgem(loom)

# Calculate enrichment scores
#this works only in the latest version (1.17):
#Install:devtools::install_github("aertslab/AUCell")

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=F, nCores=40)

#Use the version with multi-processing
cells_AUC=AUCell_calcAUC(
  regulons_new,
  cells_rankings,
  nCores = 40,
  normAUC = TRUE,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),
  verbose = TRUE
)
mat=getAUC(cells_AUC)
write.table(mat, "out/UpdatedByFreq.0.8_regulons_CellsAUC_All_cells.20231127.tsv",sep='\t',quote=F,row.names=T)

saveRDS(regulons_new, "out/UpdatedByFreq.0.8_regulons_All_cells.20231127.rds")

###Annotate regulons_new:
all_st=NULL
for (i in 1:length(regulons_new)){
    reg=names(regulons_new)[[i]]
    genes_n=length(regulons_new[[i]])
    st=cbind(reg,genes_n)
    all_st=rbind(all_st,st)
}
all_st=as.data.frame(all_st)
colnames(all_st)=c('Regulon','Genes_N')
write.table(all_st,'out/Regulons_new_annot.20230914.tsv',sep='\t',quote=F,row.names=F)



#######################

  