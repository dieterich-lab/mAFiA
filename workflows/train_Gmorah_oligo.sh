#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:turing:1
#SBATCH --mem=250GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --job-name=train_Gmorah
#SBATCH --output=/home/achan/slurm/train_Gmorah_oligo_%A.out

eval "$(/home/achan/miniconda3/condabin/conda shell.bash hook)"
conda activate MAFIA
set -e -f

prj_dir=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm
annotation=${prj_dir}/oligo_reference/oligo_mod_location.bed

############################################################################################
### Training data ##########################################################################
############################################################################################
#train_ds=A1-4
#train_ds=A5-8
train_ds=A9-12

############################################################################################
### Backbone settings ######################################################################
############################################################################################
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
ext_layer=convlayers.conv21
classifier=logistic_regression
scaler=MaxAbs
############################################################################################
############################################################################################

unm_bam=${prj_dir}/${train_ds}/spomelette_q70.bam
unm_fast5=${prj_dir}/${train_ds}/fast5_pass
mod_bam=${unm_bam}
mod_fast5=${unm_fast5}
ref=${prj_dir}/${train_ds}/ligation_ref.fasta
outdir=${prj_dir}/Gmorah_models/${train_ds}

python3 -u ${HOME}/git/Gmorah/oligo/oligo_train_binary_classifier.py \
--unm_bam_file ${unm_bam} \
--unm_fast5_dir ${unm_fast5} \
--mod_bam_file ${mod_bam} \
--mod_fast5_dir ${mod_fast5} \
--ref_file ${ref} \
--annotation ${annotation} \
--backbone_model_path ${backbone} \
--scaler ${scaler} \
--min_coverage 10 \
--classifier_model_dir ${outdir}
