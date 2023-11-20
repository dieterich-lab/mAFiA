#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500GB
#SBATCH --verbose
#SBATCH --job-name=test_mAFiA
#SBATCH --output=/home/achan/slurm/test_mAFiA_%A.out

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/100_WT_0_IVT
bam=${workspace}/chrX_q50.bam
fast5_dir=${workspace}/fast5_chrX
ref=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa
mod=/home/achan/Data/GLORI/bed_files/GLORI_chrX_ref5mer.tsv
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
classifiers=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/DRACH
output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/results/mAFiA_deep_chrX

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
