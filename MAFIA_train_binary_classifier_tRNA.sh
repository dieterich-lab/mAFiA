#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --mem=120GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_train_binary_classifier_tRNA.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

set -e -u

cd ${HOME}/git/MAFIA
WORKSPACE=/beegfs/prj/tRNA_Berlin/newBatchDec2022_Spombe
OUTDIR=${WORKSPACE}/achan/parasail_filtered_classifier_models

python3 tRNA_train_binary_classifier.py \
--wt_bam_file ${WORKSPACE}/achan/mapping/parasail/tRNA_Q.bam.filtered.sorted \
--wt_fast5_dir ${WORKSPACE}/agg_fast5_pass/tRNA_Q \
--ivt_bam_file ${WORKSPACE}/achan/mapping/parasail/tRNA_IVT.bam.filtered.sorted \
--ivt_fast5_dir ${WORKSPACE}/agg_fast5_pass/tRNA_IVT \
--ref_file ${WORKSPACE}/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta \
--mod_file ${WORKSPACE}/S_pombe_Q_only.tsv \
--max_num_reads -1 \
--min_coverage 0 \
--model_path ${HOME}/pytorch_models/tRNA_IVT/tRNA_IVT-epoch29.torch \
--classifier logistic_regression \
--classifier_model_dir ${OUTDIR} \
--outfile ${OUTDIR}/logreg_auc.tsv
