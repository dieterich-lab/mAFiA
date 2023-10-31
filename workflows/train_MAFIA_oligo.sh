#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --mem=120GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --job-name=train_MAFIA_oligo
#SBATCH --output=/home/achan/slurm/train_MAFIA_oligo_%A.out

eval "$(conda shell.bash hook)"
#eval "$(/home/achan/miniconda3/condabin/conda shell.bash hook)"
conda activate MAFIA
set -e -f
cd ${HOME}/git/mAFiA_dev

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A
DATA_DIR=${PRJ_DIR}/oligo

############################################################################################
### Training data ##########################################################################
############################################################################################
### ISA ###
#TRAIN_DATASET=ISA_runs1-3

### WUE ###
#TRAIN_DATASET=WUE_batches1-2

### ISA-WUE ###
#TRAIN_DATASET=ISA-WUE

### ISA 12 more ###
TRAIN_DATASET=ISA_mix1

############################################################################################
### Backbone settings ######################################################################
############################################################################################
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXT_LAYER=convlayers.conv21
CLASSIFIER=logistic_regression
SCALER=MaxAbs
############################################################################################
############################################################################################

UNM_BAM=${DATA_DIR}/${TRAIN_DATASET}_A/spomelette_q70.bam
UNM_FAST5=${DATA_DIR}/${TRAIN_DATASET}_A/fast5
MOD_BAM=${DATA_DIR}/${TRAIN_DATASET}_m6A/spomelette_q70.bam
MOD_FAST5=${DATA_DIR}/${TRAIN_DATASET}_m6A/fast5
REF=${DATA_DIR}/${TRAIN_DATASET}_A_m6A/ligation_ref.fasta
OUTDIR=${PRJ_DIR}/MAFIA_classifiers/${TRAIN_DATASET}

python3 -u ${HOME}/git/mAFiA_dev/oligo/oligo_train_binary_classifier.py \
--unm_bam_file ${UNM_BAM} \
--unm_fast5_dir ${UNM_FAST5} \
--mod_bam_file ${MOD_BAM} \
--mod_fast5_dir ${MOD_FAST5} \
--ref_file ${REF} \
--backbone_model_path ${BACKBONE_MODEL} \
--scaler ${SCALER} \
--min_coverage 10 \
--classifier_model_dir ${OUTDIR}
