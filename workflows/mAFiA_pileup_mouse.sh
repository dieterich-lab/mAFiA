#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup_mouse
#SBATCH --output=/home/achan/slurm/pileup_mouse_%A_chr%a.out

#ds=100_WT_0_IVT
#ds=0_WT_100_IVT
#chr=X
if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
chr="X"
else
chr=${SLURM_ARRAY_TASK_ID}
fi

workspace=/prj/Dewenter_TAC_Backs_lab/achan/psico-mAFiA_results/${ds}
bam=${workspace}/chr${chr}/mAFiA.reads.bam
mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.chr${chr}.bed
output=${workspace}/chr${chr}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
--bam_file ${bam} \
--mod_file ${mod} \
--min_coverage 20 \
--out_dir ${output} \
--num_jobs 36