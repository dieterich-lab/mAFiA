#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=Gmorah_HepG2
#SBATCH --output=/home/achan/slurm/Gmorah_HepG2_%A_chr%a.out

ds=rep1

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
chr="X"
else
chr=${SLURM_ARRAY_TASK_ID}
fi

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HepG2/${ds}
bam=${workspace}/chr${chr}/sorted.chr${chr}.bam
fast5_dir=${workspace}/chr${chr}/fast5

mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102/Gm.GRCh38_102.chr${chr}.bed

backbone=${HOME}/git/mAFiA/models/RODAN_HEK293_IVT.torch
classifiers=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gmorah

output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/Gmorah/HepG2/${ds}/chr${chr}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_process_reads_parallel.py \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--mod_file ${mod} \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--num_jobs 4 \
--batchsize 1024 \
--out_dir ${output}
