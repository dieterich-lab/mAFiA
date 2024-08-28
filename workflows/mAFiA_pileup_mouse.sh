#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup_mouse
#SBATCH --output=/home/achan/slurm/pileup_mouse_%A_chr%a.out

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
chr="X"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 24 ]]
then
chr="MT"
else
chr=${SLURM_ARRAY_TASK_ID}
fi

workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/${ds}
bam=${workspace}/chr${chr}/mAFiA.reads.bam
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.chr${chr}.bed
output=${workspace}/chr${chr}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
--bam_file ${bam} \
--mod_file ${mod} \
--min_coverage 20 \
--out_dir ${output} \
--num_jobs 36
