#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --mem=192GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_test_binary_classifier_mRNA_NoNorm_PassFail.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian
DATASET=HEK293A_WT
FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all
REF=${HOME}/Data/genomes/GRCh38_96.fa
BAM=${WORKSPACE}/mapping/${DATASET}_sorted.bam
MOD_FILE=${HOME}/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv
BACKBON_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXTRACTION_LAYER=convlayers.conv21
CLASSIFIER_MODEL_DIR=${WORKSPACE}/MAFIA_classifiers/A_m6A_NoNorm_PassFail
OUTFILE=${WORKSPACE}/results/res_${DATASET}_NoNorm_PassFail.tsv

set -e -u

cd ${HOME}/git/MAFIA

python3 ${HOME}/git/MAFIA/test_binary_classifier.py \
--test_bam_file ${BAM} \
--test_fast5_dir ${FAST5_DIR} \
--ref_file ${REF} \
--mod_file ${MOD_FILE} \
--max_num_reads -1 \
--min_coverage 10 \
--backbone_model_path ${BACKBON_MODEL} \
--extraction_layer ${EXTRACTION_LAYER} \
--feature_width 0 \
--classifier logistic_regression \
--classifier_model_dir ${CLASSIFIER_MODEL_DIR} \
--use_opt_thresh False \
--outfile ${OUTFILE}
