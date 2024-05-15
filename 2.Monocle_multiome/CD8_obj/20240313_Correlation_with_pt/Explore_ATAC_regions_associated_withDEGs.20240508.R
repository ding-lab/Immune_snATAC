library(Signac)
library(Seurat)
library(reshape)
library(reshape2)
library(tidyverse)

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp


ATAC_m=read.table('out/ATAC_immune_slot_data.CD8.20240312.tsv',sep='\t',header=T)

annot=read_delim('for_heatmap/Peaks_immune_snATAC.Annotated.20240430.tsv',delim='\t')
annot=as.data.frame(annot)
annot_s=annot[,c('peak_ID','SYMBOL')]
colnames(annot_s)=c('Feature','Gene')

degs=read.table('out/Pearson_corr_PT_with_RNA_expr.min500Cells.20240430.tsv',sep='\t',header=T)
degs_s=degs[degs$FDR<0.05,]

annot_s2=annot_s[annot_s$Gene %in% degs_s$Feature,]

ATAC_m=ATAC_m[ATAC_m$Feature %in% annot_s2$Feature,]

#The length of peaks in ATAC_m matrix is smaller, because we removed peaks with coverage in <= 10 cells.

pt_data=read.table('../Monocle_RDS/All_PTs_monocle_objs.10comp.CD8.After_reOrdering.state16_root.20240313.tsv',sep='\t',header=T)
pt_data=pt_data[,c('sample_name','Pseudotime','cell_type')]
colnames(pt_data)[1]<-'Barcode'

library(doParallel)
registerDoParallel(cores=50)

#Estimate number that should be processed by a single core, when using 70 cores -- 70 is too much memory, try 50:
all_unique_ids=unique(ATAC_m$Feature)
n_per_core=round(length(all_unique_ids)/50)+1

all_st=NULL
all_st<-foreach(ids_set_n=c(1:50)) %dopar% {
    first_id=(ids_set_n-1)*n_per_core+1
    last_id=min((ids_set_n)*n_per_core,length(all_unique_ids))
    sel_ids=all_unique_ids[first_id:last_id]
    ATAC_m2=ATAC_m[ATAC_m$Feature %in% sel_ids,]

all_st_per_core=NULL
    for (peak in unique(ATAC_m2$Feature)){
       atac_s=ATAC_m2[ATAC_m2$Feature==peak,]
           length_non_zero=nrow(atac_s[atac_s$Expression!=0,])
           if (length_non_zero>=50){
              atac_s=merge(atac_s,pt_data)
              st=cor.test(atac_s$Pseudotime,atac_s$Expression)
              stat=c(peak,st$p.value,st$estimate,"ATAC")
              all_st_per_core=rbind(all_st_per_core,stat)
              }
           }
           return(all_st_per_core)

}

all_st_f=NULL
for (i in 1:50){
    all_st_1=as.data.frame(all_st[[i]])
    all_st_f=rbind(all_st_f,all_st_1)
}

all_st_f=as.data.frame(all_st_f)
colnames(all_st_f)<-c('Feature','P_value','Estimate','Data')
all_st_f$FDR=p.adjust(all_st_f$P_value, method='fdr')

all_st_f2=merge(all_st_f,annot_s2,all.x=T)
write.table(all_st_f2, 'out/Pearson_corr_PT_with_ATAC_expr.min50Cells.forPeaks_assoc_wDEGs.20240508.tsv',sep='\t',quote=F,row.names=F)
