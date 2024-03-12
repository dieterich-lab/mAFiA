#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:turing:1
#SBATCH --mem=120GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --job-name=train_Gmorah_oligo
#SBATCH --output=/home/achan/slurm/train_Gmorah_oligo_%A.out

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

set -e -f

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm
DATA_DIR=${PRJ_DIR}/oligo

############################################################################################
### Backbone settings ######################################################################
############################################################################################
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXT_LAYER=convlayers.conv21
CLASSIFIER=logistic_regression
SCALER=MaxAbs

############################################################################################
############################################################################################
orig=NMG
#run=mix38_mix44
unm_ds=unm
mod_ds=mod

UNM_BAM=${DATA_DIR}/${orig}_${run}_${unm_ds}/spomelette_q70.bam
UNM_FAST5=${DATA_DIR}/${orig}_${run}_${unm_ds}/fast5
MOD_BAM=${DATA_DIR}/${orig}_${run}_${mod_ds}/spomelette_q70.bam
MOD_FAST5=${DATA_DIR}/${orig}_${run}_${mod_ds}/fast5
REF=${DATA_DIR}/ligation_ref/ligation_ref_${run}.fasta
OUTDIR=${PRJ_DIR}/Gmorah_classifiers/${orig}_${run}

python3 -u ${HOME}/git/mAFiA_dev/oligo/oligo_train_binary_classifier.py \
--unm_bam_file ${UNM_BAM} \
--unm_fast5_dir ${UNM_FAST5} \
--mod_bam_file ${MOD_BAM} \
--mod_fast5_dir ${MOD_FAST5} \
--ref_file ${REF} \
--backbone_model_path ${BACKBONE_MODEL} \
--scaler ${SCALER} \
--min_coverage 1 \
--classifier_model_dir ${OUTDIR}
