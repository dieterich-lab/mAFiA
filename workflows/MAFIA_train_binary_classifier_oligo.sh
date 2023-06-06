#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --mem=120GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_train_binary_classifier_oligo_%A.out

eval "$(conda shell.bash hook)"
#eval "$(/home/achan/miniconda3/condabin/conda shell.bash hook)"
conda activate MAFIA
set -e -f
cd ${HOME}/git/MAFIA

DATA_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian

############################################################################################
### feature settings #######################################################################
############################################################################################
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXT_LAYER=convlayers.conv21
CLASSIFIER=logistic_regression
SCALER=MaxAbs

############################################################################################
### source data ############################################################################
############################################################################################
#### Claudia batches 1 & 2 ###
#TRAIN_DATASET=WUE_batch1
#TRAIN_DATASET=WUE_batch2
#TRAIN_DATASET=WUE_batches1+2

#UNM_BAM=${DATA_DIR}/${TRAIN_DATASET}_A/spomelette_q80.bam
#UNM_FAST5=${DATA_DIR}/${TRAIN_DATASET}_A/fast5
#MOD_BAM=${DATA_DIR}/${TRAIN_DATASET}_m6A/spomelette_q80.bam
#MOD_FAST5=${DATA_DIR}/${TRAIN_DATASET}_m6A/fast5
#REF=${DATA_DIR}/reference/${TRAIN_DATASET}_ref_recon.fa
#OUTDIR=${DATA_DIR}/MAFIA_classifiers/${TRAIN_DATASET}_spomelette_${CLASSIFIER}_${SCALER}

### Isabel RL top 6 ###
TRAIN_DATASET=ISA_run1
UNM_BAM=${DATA_DIR}/${TRAIN_DATASET}_A/spomelette_q70.bam
UNM_FAST5=${DATA_DIR}/${TRAIN_DATASET}_A/fast5
MOD_BAM=${DATA_DIR}/${TRAIN_DATASET}_m6A/spomelette_q70.bam
MOD_FAST5=${DATA_DIR}/${TRAIN_DATASET}_m6A/fast5
REF=${DATA_DIR}/reference/${TRAIN_DATASET}_ligation_ref.fa
OUTDIR=${DATA_DIR}/MAFIA_classifiers/${TRAIN_DATASET}_spomelette_${CLASSIFIER}_${SCALER}

### Isabel RL Mix 1 & 3 ###
#UNM_BAM=${DATA_DIR}/RL_Mix1_A_RTA/mapped.bam
#UNM_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix1_A_RTA/*/fast5_*"
#MOD_BAM=${DATA_DIR}/RL_Mix3_m6A_RTA/mapped.bam
#MOD_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix3_m6A_RTA/*/fast5_*"
#REF=${DATA_DIR}/reference/RL_Mix1_Mix3_blocks8.fasta
#OUTDIR=${DATA_DIR}/MAFIA_classifiers/RL_run2_Mix1_Mix3_${SCALER}

### Isabel RL Mix 2 & 4 ###
#UNM_BAM=${DATA_DIR}/RL_Mix2_A_RTA/mapped.bam
#UNM_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix2_A_RTA/*/fast5_*"
#MOD_BAM=${DATA_DIR}/RL_Mix4_m6A_RTA/mapped.bam
#MOD_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix4_m6A_RTA/*/fast5_*"
#REF=${DATA_DIR}/reference/RL_Mix2_Mix4_blocks8.fasta
#OUTDIR=${DATA_DIR}/MAFIA_classifiers/RL_run2_Mix2_Mix4_${SCALER}

python3 -u oligo_train_binary_classifier.py \
--unm_bam_file ${UNM_BAM} \
--unm_fast5_dir ${UNM_FAST5} \
--mod_bam_file ${MOD_BAM} \
--mod_fast5_dir ${MOD_FAST5} \
--ref_file ${REF} \
--backbone_model_path ${BACKBONE_MODEL} \
--scaler ${SCALER} \
--min_coverage 10 \
--classifier_model_dir ${OUTDIR}
