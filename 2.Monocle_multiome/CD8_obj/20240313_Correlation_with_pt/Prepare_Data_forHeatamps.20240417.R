library(Signac)
library(Seurat)
library(reshape)
library(reshape2)
library(monocle)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

#data_dir='/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/DATA/data_freeze_v1/allATAC_new_peaks'

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp


obj=readRDS('../PanImmune_RNA_ATAC.v8.0.CD8.500cells_perCellType.rds')

RNA <- GetAssayData(object=obj, assay='SCT', slot='data')

degs=read.table('out/Pearson_corr_PT_with_RNA_expr.min500Cells.20240430.tsv',sep='\t',header=T)
degs=degs[degs$FDR<0.05,]

RNA_s=RNA[rownames(RNA) %in% degs$Feature,]

write.table(RNA_s, 'for_heatmap/RNA_CD8_correlated_genes_matrix.20240501.tsv',sep='\t',quote=F,row.names=T,col.names=T)


#Now also extract data for ATAC assay:
ATAC <- GetAssayData(object=obj, assay='ATAC_immune', slot='data')

dacrs=read.table('out/Pearson_corr_PT_with_ATAC_expr.min500Cells.20240430.tsv',sep='\t',header=T)
dacrs=dacrs[dacrs$FDR<0.05,]

ATAC_s=ATAC[rownames(ATAC) %in% dacrs$Feature,]

write.table(ATAC_s, 'for_heatmap/ATAC_CD8_correlated_genes_matrix.20240501.tsv',sep='\t',quote=F,row.names=T,col.names=T)


#Now also annotate peaks:

peaks_1=StringToGRanges(rownames(ATAC), sep = c("-", "-"))

#Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")

anno=as.data.frame(peakAnno)
anno$peak_ID=rownames(ATAC)

write.table(anno, "for_heatmap/Peaks_immune_snATAC.Annotated.20240430.tsv",sep='\t',row.names=F,quote=F)


#Now also extract data for the latest list of peaks (associated with DEGs):

obj=readRDS('../PanImmune_RNA_ATAC.v8.0.CD8.500cells_perCellType.rds')

ATAC <- GetAssayData(object=obj, assay='ATAC_immune', slot='data')

dacrs=read.table('out/Pearson_corr_PT_with_ATAC_expr.min50Cells.forPeaks_assoc_wDEGs.20240508.tsv',sep='\t',header=T)
#dacrs=dacrs[dacrs$FDR<0.05,]

ATAC_s=ATAC[rownames(ATAC) %in% dacrs$Feature,]

write.table(ATAC_s, 'for_heatmap/ATAC_CD8_correlated_matrix.associated_wDEGs.20240508.tsv',sep='\t',quote=F,row.names=T,col.names=T)
