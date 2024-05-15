while getopts ":i:" opt; do
  case ${opt} in
    i )
      iteration=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      ;;
  esac
done
shift $((OPTIND -1))

#STEP 1: correlation analysis to infer potential TF targets based on gene expression
#in bash


disease='T_cell_PancanObj'
wdir="/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/Analysis/1.GRN_analysis/Run.v.20230912/Object"
wdir_2="/diskmnt/Projects/Users/nvterekhanova/Immune_snATAC/Analysis/1.GRN_analysis/Run.v.20231027/Object"

obj_path=$wdir"/data/T_cell_PancanObj.v2.loom"
data_path="/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/9.GRN_analysis/Data"

#pyscenic grn $obj_path $data_path/JAPSAR2020_unique_motifs.txt -o $wdir/$iteration/adj.csv --num_workers 80 --sparse



#STEP 2-3: Regulon prediction aka cisTarget from CLI

#https://github.com/aertslab/pySCENIC/issues/275
pyscenic ctx $wdir/$iteration/adj.csv \
    $data_path/cisTarget_databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather $data_path/cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
    --annotations_fname $data_path/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname $obj_path \
    --output $wdir_2/$iteration/reg.csv \
    --mask_dropouts \
    --num_workers 80 \
    --all_modules

