#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(monocle)
library(Signac)
library(Seurat)
library(RColorBrewer)
library(ggplot2)


n_components=10

#data_dir='/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/DATA/data_freeze_v1/allATAC_new_peaks'
run_date='20240313'

#Read RDS object:
obj_downs=readRDS('PanImmune_RNA_ATAC.v8.0.CD4.500cells_perCellType.rds')
DefaultAssay(obj_downs)<-'ATAC_immune'
Idents(obj_downs)=obj_downs$cell_type

###Use slot counts:
x=obj_downs@assays$ATAC_immune@counts

###Change from rowSums(x!=0):
peaks_cov=rowSums(x)

peaks_cov=peaks_cov[order(-peaks_cov)]

###Select the top covered across cells peaks 
peaks_sel=names(peaks_cov[1:50000])

x1=x[rownames(x) %in% peaks_sel,]
obj_downs[['for_m']]<-CreateAssayObject(x1)


#Load Seurat object
seurat_object <- obj_downs

#Extract data, phenotype data, and feature data from the SeuratObject
#Use slot counts instead of data:
data <- as(as.matrix(seurat_object@assays$for_m@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)



#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())


monocle_cds<-estimateSizeFactors(monocle_cds)
monocle_cds<-estimateDispersions(monocle_cds)



monocle_cds <- reduceDimension(monocle_cds, method = 'DDRTree', max_components =n_components)
monocle_cds <- orderCells(monocle_cds)

monocle_cds$ID=monocle_cds$cell_type
monocle_cds$ID=factor(monocle_cds$cell_type)



cell_type_cols=c(brewer.pal(n=6, name='Dark2'))

p1=plot_cell_trajectory(monocle_cds, color_by = "ID", cell_size = 1) +
scale_color_manual(values=cell_type_cols)

pdf(paste('plots/CD4_obj_ATAC_500_cells_perCelltype_', n_components,
'_maxComp.',run_date,'.pdf',sep=''),width=9,height=7,useDingbats=F)
print(p1)
dev.off()

saveRDS(monocle_cds,paste('Monocle_RDS/CD4_obj_ATAC_500_cells_perCelltype_10_maxComp.',run_date,'.rds',sep=''))

#Now try to color by cancer type:
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_immune_ATAC_data_freeze/v8.0/Colors_panatac_v3.0.rds')
dis_cols=colors$Cancer

monocle_cds=readRDS(paste('Monocle_RDS/CD4_obj_ATAC_500_cells_perCelltype_10_maxComp.',run_date,'.rds',sep=''))

p1=plot_cell_trajectory(monocle_cds, color_by = "Cancer", cell_size = 1) +
scale_color_manual(values=dis_cols)

pdf(paste('plots/CD4_obj_ATAC_500_cells_perCelltype_', n_components,
'_maxComp_byCancer.',run_date,'.pdf',sep=''),width=9,height=7,useDingbats=F)
print(p1)
dev.off()

#Get barcodes and their PT:
monocle_cds=readRDS(paste('Monocle_RDS/CD4_obj_ATAC_500_cells_perCelltype_10_maxComp.',run_date,'.rds',sep=''))


p1=plot_cell_trajectory(monocle_cds, color_by = "cell_type", cell_size = 1)
tab <- ggplot_build(p1)[["plot"]][["data"]]
tab_s=tab[,c('sample_name','data_dim_1', 'data_dim_2','Pseudotime','cell_type')]

#tab=cbind(sampleNames(monocle_cds),monocle_cds$Pseudotime,monocle_cds$State,monocle_cds$cell_type)

write.table(tab_s,paste('Monocle_RDS/All_PTs_monocle_objs.10comp.CD4.',run_date,'.tsv',sep=''),sep='\t',quote=F,row.names=F)