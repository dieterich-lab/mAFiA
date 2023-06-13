#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=test_MAFIA_oligo
#SBATCH --output=/home/achan/slurm/test_MAFIA_oligo_%A.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

########################################################################################################################
#TRAIN_DATASET=ISA_run1
#TRAIN_DATASET=WUE_batches1+2
########################################################################################################################
#TEST_DATASET=ISA_run1_A
#TEST_DATASET=ISA_run1_m6A
#TEST_DATASET=WUE_batch1_A
########################################################################################################################
#TEST_DATASET=RL_M4_M5
#TEST_DATASET=RL_M4_M5star
#TEST_DATASET=RL_M4star_M5
#TEST_DATASET=RL_M4star_M5star
########################################################################################################################

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
BAM=${PRJ_DIR}/${TEST_DATASET}/spomelette_q70.bam
FAST5_DIR=${PRJ_DIR}/${TEST_DATASET}/fast5
LIG_REF=${PRJ_DIR}/${TEST_DATASET}/ligation_ref.fa
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
CLASSIFIER_MODEL_DIR=${PRJ_DIR}/MAFIA_classifiers/${TRAIN_DATASET}
OUTFILE=${PRJ_DIR}/results/res_train_${TRAIN_DATASET}_test_${TEST_DATASET}_q70.tsv
########################################################################################################################

set -e -f

python3 -u ${HOME}/git/MAFIA/oligo_test_binary_classifier.py \
--test_bam_file ${BAM} \
--test_fast5_dir ${FAST5_DIR} \
--ref_file ${LIG_REF} \
--backbone_model_path ${BACKBONE_MODEL} \
--classifier_model_dir ${CLASSIFIER_MODEL_DIR} \
--outfile ${OUTFILE}
