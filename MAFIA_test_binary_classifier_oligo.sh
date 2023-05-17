#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=MAFIA_test_oligo
#SBATCH --output=/home/achan/slurm/MAFIA_test_oligo_%A.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
TRAIN_DATASET=WUE_combined
TEST_DATASET=ISA_run1_A
#TEST_DATASET=ISA_run1_m6A

BAM=${PRJ_DIR}/${TEST_DATASET}/spomelette_q80.bam
FAST5_DIR=${PRJ_DIR}/${TEST_DATASET}/fast5
REF=${PRJ_DIR}/${TEST_DATASET}/ref_recon.fa
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
CLASSIFIER_MODEL_DIR=${PRJ_DIR}/MAFIA_classifiers/${TRAIN_DATASET}
OUTFILE=${PRJ_DIR}/results/res_train_${TRAIN_DATASET}_test_${TEST_DATASET}.tsv

set -e -u -f

python3 ${HOME}/git/MAFIA/oligo_test_binary_classifier.py \
--test_bam_file ${BAM} \
--test_fast5_dir ${FAST5_DIR} \
--ref_file ${REF} \
--backbone_model_path ${BACKBONE_MODEL} \
--classifier_model_dir ${CLASSIFIER_MODEL_DIR} \
--outfile ${OUTFILE}
