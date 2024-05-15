library(Signac)
library(Seurat)
library(Matrix)
library(SCENIC)
library(SCopeLoomR)


wdir='/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/Analysis/1.GRN_analysis/Run.v.20230912/Object/'

data_dir='/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/DATA/data_freeze_v1'

rna=readRDS(paste(data_dir,'/PanImmune_integrated_allRNA_by_cancer_chemistry_v8.0_df_T_cells_NK_no_doublets_res1.4.rds',sep=''))

###Downsample cells:

Idents(rna)=rna$seurat_clusters
rna_downs=subset(x = rna, downsample = 1500)

assay=GetAssay(rna_downs,assay='RNA')
data=assay@data
genes_list=read.csv(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/9.GRN_analysis/',
'Data/hgnc-gene_list.txt',sep=''),sep='\t')
genes_prot_c=genes_list[genes_list$Locus.type=='Gene with protein product',]


data=data[rownames(data) %in% genes_prot_c$Symbol,]

###2022-01-23, remove RPL/RPS/MT- genes:
genes=rownames(data)
genes_r=genes[grepl('^MT-|^RPL|^RPS', genes)]
data=data[!(rownames(data) %in% genes_r),]

cellInfo_1 <- data.frame(Seurat_clusters=rna_downs$seurat_clusters, cell_type_v8.0_Tcell=rna_downs$cell_type_v8.0_Tcell)
exprMatrix=data

dir.create(paste(wdir,"/data",sep=''))

loom <- build_loom(paste(wdir,"/data/T_cell_PancanObj.loom",sep=''), dgem=exprMatrix)

loom <- add_cell_annotation(loom, cellInfo_1)
close_loom(loom)

cellInfo_1$Barcode=row.names(cellInfo_1)

write.table(cellInfo_1,paste(wdir,'/T_cell_PancanObj_1500cellsPer_seurat_cluster.20230912.tsv',sep=''),
sep='\t',quote=F,row.names=F)
