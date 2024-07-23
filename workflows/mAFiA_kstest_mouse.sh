#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=kstest_mouse
#SBATCH --output=/home/achan/slurm/kstest_mouse_%A.out

set -e -u

workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC
bam1=${workspace}/SHAM_day${day}/chrALL.mAFiA.reads.bam
bam2=${workspace}/TAC_day${day}/chrALL.mAFiA.reads.bam
bed1=${workspace}/SHAM_day${day}/chrALL.mAFiA.sites.bed
bed2=${workspace}/TAC_day${day}/chrALL.mAFiA.sites.bed
out_dir=${workspace}/KS_test
out_filename=SHAM_TAC_day${day}.bed

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_KS_test_two_samples.py \
--bam_file_1 ${bam1} \
--bam_file_2 ${bam2} \
--bed_file_1 ${bed1} \
--bed_file_2 ${bed2} \
--min_coverage 50 \
--out_dir ${out_dir} \
--num_jobs 36 \
--out_filename ${out_filename}
