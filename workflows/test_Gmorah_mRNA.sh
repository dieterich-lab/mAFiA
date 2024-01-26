#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=250GB
#SBATCH --verbose
#SBATCH --job-name=test_Gmorah_mRNA
#SBATCH --output=/home/achan/slurm/test_Gmorah_mRNA_%A.out

#workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/HEK293
bam=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/100_WT_0_IVT/genome_filtered_q50.bam
fast5_dir=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/100_WT_0_IVT/fast5
ref=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa
mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/Gm_sites_mRNA_HeLa_HepG2_merged.bed
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
classifiers=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Gmorah_models/combined
output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/HEK293/Gmorah

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
