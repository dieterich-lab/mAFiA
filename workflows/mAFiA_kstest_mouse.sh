#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=kstest_mouse
#SBATCH --output=/home/achan/slurm/kstest_mouse_%A_chr%a.out

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
chr="X"
else
chr=${SLURM_ARRAY_TASK_ID}
fi

workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/${ds}
bam1=${workspace}/SHAM_day${day}/chr${chr}/mAFiA.reads.bam
bam2=${workspace}/TAC_day${day}/chr${chr}/mAFiA.reads.bam
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.chr${chr}.bed
out_dir=${workspace}/SHAM_TAC_day${day}/chr${chr}
out_filename=KSTest_SHAM_TAC_day${day}.bed

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_KS_test_two_samples.py \
--bam_file_1 ${bam1} \
--bam_file_2 ${bam2} \
--mod_file ${mod} \
--out_dir ${out_dir} \
--num_jobs 36 \
--out_filename
