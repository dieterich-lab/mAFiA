#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=160GB
#SBATCH --verbose
#SBATCH --job-name=DRACH_chrX
#SBATCH --output=/home/achan/slurm/DRACH_chrX.out

chr=X

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/100_WT_0_IVT
bam=${workspace}/chr${chr}/sorted.chr${chr}.bam
fast5_dir=${workspace}/chr${chr}/fast5
ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/GRCh38_102/GRCh38_102.chr${chr}.fa
mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/site_annotations/DRACH.GRCh38_102.chr${chr}.bed
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
classifiers=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/DRACH
output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH/chr${chr}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

test_mAFiA \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--ref_file ${ref} \
--mod_file ${mod} \
--min_coverage 10 \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--mod_prob_thresh 0.5 \
--out_dir ${output}
