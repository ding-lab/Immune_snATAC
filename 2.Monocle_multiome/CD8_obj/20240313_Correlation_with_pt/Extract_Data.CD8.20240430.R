library(Signac)
library(Seurat)
library(reshape)
library(reshape2)

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

#Read RDS object:

#obj=readRDS('/diskmnt/Projects/snATAC_analysis/immune/obj/v8.0/multiome_int_cancer_atac/PanImmune_int_RNA_ATAC_Tcell_multiome_int_atac_res1.6.rds')
obj=readRDS('../PanImmune_RNA_ATAC.v8.0.CD8.500cells_perCellType.rds')

RNA <- GetAssayData(object=obj, assay='SCT', slot='data')
genes_cov=apply(RNA,1,function(x) length(x[x>0]))
genes_cov_s=genes_cov[genes_cov>10]

RNA_2=RNA[rownames(RNA) %in% names(genes_cov_s),]
RNA_m=reshape2::melt(as.matrix(RNA_2))
colnames(RNA_m)=c('Feature','Barcode','Expression')
RNA_m$Data='RNA'

ATAC <- GetAssayData(object=obj, assay='ATAC_immune', slot='data')
peaks_cov=apply(ATAC,1,function(x) length(x[x>0]))
peaks_cov_s=peaks_cov[peaks_cov>10]

ATAC_2=ATAC[rownames(ATAC) %in% names(peaks_cov_s),]
ATAC_m=reshape2::melt(as.matrix(ATAC_2))
colnames(ATAC_m)=c('Feature','Barcode','Expression')
ATAC_m$Data='ATAC'

write.table(RNA_m, 'out/RNA_SCT_slot_data.CD8.20240312.tsv',sep='\t',row.names=F, quote=F)
write.table(ATAC_m, 'out/ATAC_immune_slot_data.CD8.20240312.tsv',sep='\t',row.names=F, quote=F)

#Now perform correlation analysis between the new PT and gene expr/ATAC accessibilities:
RNA_m=read.table('out/RNA_SCT_slot_data.CD8.20240312.tsv',sep='\t',header=T)
pt_data=read.table('../Monocle_RDS/All_PTs_monocle_objs.10comp.CD8.After_reOrdering.state16_root.20240313.tsv',sep='\t',header=T)
pt_data=pt_data[,c('sample_name','Pseudotime','cell_type')]
colnames(pt_data)[1]<-'Barcode'

library(doParallel)
registerDoParallel(cores=50)

#Estimate number that should be processed by a single core, when using 70 cores -- 70 is too much memory, try 50:
all_unique_ids=unique(RNA_m$Feature)
n_per_core=round(length(all_unique_ids)/50)+1

all_st=NULL
all_st<-foreach(ids_set_n=c(1:50)) %dopar% {
    first_id=(ids_set_n-1)*n_per_core+1
    last_id=min((ids_set_n)*n_per_core,length(all_unique_ids))
    sel_ids=all_unique_ids[first_id:last_id]
    RNA_m2=RNA_m[RNA_m$Feature %in% sel_ids,]

all_st_per_core=NULL
    for (gene in unique(RNA_m2$Feature)){
    rna_s=RNA_m2[RNA_m2$Feature==gene,]
    length_non_zero=nrow(rna_s[rna_s$Expression!=0,])
    if (length_non_zero>=500){
        rna_s=merge(rna_s,pt_data)
    	st=cor.test(rna_s$Pseudotime,rna_s$Expression)
    	stat=c(gene,st$p.value,st$estimate,"RNA")
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
write.table(all_st_f, 'out/Pearson_corr_PT_with_RNA_expr.min500Cells.20240430.tsv',sep='\t',quote=F,row.names=F)

#Now do the same with ATAC data:
#Now perform correlation analysis between the new PT and gene expr/ATAC accessibilities:

ATAC_m=read.table('out/ATAC_immune_slot_data.CD8.20240312.tsv',sep='\t',header=T)
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
    	   if (length_non_zero>=500){
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
write.table(all_st_f, 'out/Pearson_corr_PT_with_ATAC_expr.min500Cells.20240430.tsv',sep='\t',quote=F,row.names=F)


