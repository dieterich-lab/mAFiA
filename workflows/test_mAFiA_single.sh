#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=test_mAFiA
#SBATCH --output=/home/achan/slurm/test_mAFiA_%A.out

workspace=/prj/TRR319_RMaP/Project_B01/Adrian/JK_HEK293_DMSO_merged
bam=${workspace}/genome_filtered_q50.bam
fast5_dir=${workspace}/fast5
ref=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa
mod=/prj/TRR319_RMaP/Project_B01/Adrian/possible_IR_sites_edge50.tsv
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
classifiers=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/release_v1/models/classifiers
output=${workspace}/mAFiA

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

test_mAFiA \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--ref_file ${ref} \
--mod_file ${mod} \
--min_coverage 20 \
--max_num_reads 1000 \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--mod_prob_thresh 0.5 \
--out_dir ${output}
