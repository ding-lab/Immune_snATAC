###Downsample object, need to use Signac v.1.8.

library(Signac)
library(Seurat)

n_components=10

#data_dir='/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/DATA/data_freeze_v1/allATAC_new_peaks'

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

#Read RDS object:

obj=readRDS('/diskmnt/Projects/snATAC_analysis/immune/obj/v8.0/multiome_int_cancer_atac/PanImmune_int_RNA_ATAC_Tcell_multiome_int_atac_res1.6.rds')

Idents(obj) <- obj$cell_type_v8.6_multi

DefaultAssay(obj)<-'ATAC_immune'

obj$cell_type = obj$cell_type_v8.6_multi

obj_s=subset(obj, cell_type %in% c('CD8 CTL ZNF683','CD8 T-cells exhausted GZMK','CD8 T-cells exhausted NR4A1',
'CD8 T-cells GZMK','CD8 T-cells GZMK TCF7 LEF1','CD8 T-cells naive','CD8 Trm','CD8 Trm exhausted ENTPD1'))
obj_s=subset(obj_s, Cancer %in% c('BRCA','ccRCC','CESC','CRC','GBM','HNSCC','OV','PDAC','SKCM','UCEC'))

Idents(obj_s)=obj_s$cell_type


#Downsample to reduce the calculation:
obj_downs=subset(x = obj_s, downsample=500) 

saveRDS(obj_downs, 'PanImmune_RNA_ATAC.v8.0.CD8.500cells_perCellType.rds')