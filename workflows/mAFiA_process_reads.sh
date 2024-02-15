#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=psico-mAFiA_HEK293_WT
#SBATCH --output=/home/achan/slurm/psico-mAFiA_HEK293_WT_chr%a.out

ds=100_WT_0_IVT
#ds=0_WT_100_IVT
#chr=X
chr=${SLURM_ARRAY_TASK_ID}

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/${ds}
bam=${workspace}/chr${chr}/sorted.chr${chr}.bam
fast5_dir=${workspace}/chr${chr}/fast5

mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102/m6A.psi.GRCh38_102.chr${chr}.bed
#mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102/BID_GLORI.chrX.tsv

backbone=${HOME}/git/mAFiA/models/RODAN_HEK293_IVT.torch
classifiers=${HOME}/git/mAFiA/models/psi-co-mAFiA

output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/psico-mAFiA/HEK293/${ds}/chr${chr}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_process_reads_parallel.py \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--mod_file ${mod} \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--num_jobs 4 \
--out_dir ${output}
