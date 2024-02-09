#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=180GB
#SBATCH --verbose
#SBATCH --job-name=mAFiA_multimod_chrX
#SBATCH --output=/home/achan/slurm/mAFiA_multimod_chrX.out

ds=100_WT_0_IVT
chr=X

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/${ds}
bam=${workspace}/chr${chr}/sorted.chr${chr}.bam
fast5_dir=${workspace}/chr${chr}/fast5

ref=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa
mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102/BID_GLORI.chrX.tsv

backbone=${HOME}/git/mAFiA/models/RODAN_HEK293_IVT.torch
classifiers=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/psi-co-mAFiA

output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/mAFiA_multimod

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_process_reads.py \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--ref_file ${ref} \
--mod_file ${mod} \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--out_dir ${output}
