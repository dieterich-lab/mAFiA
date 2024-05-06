#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup
#SBATCH --output=/home/achan/slurm/pileup_%A_chr%a.out

#ds=HEK_siCtrl_input_rep1
#ds=HEK_siMETTL3_input_rep1
#ds=HEK_siTRUB1_input_rep1

#ds=100_WT_0_IVT
#ds=0_WT_100_IVT
#ds=Mettl3-KO

#chr=X
if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
chr="X"
else
chr=${SLURM_ARRAY_TASK_ID}
fi

workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/${ds}
#workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HeLa
#workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/NanoSPA/${ds}
#workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HeLa_SRR28796313

bam=${workspace}/chr${chr}/mAFiA.reads.bam
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102/m6A.psi.GRCh38_102.chr${chr}.bed
output=${workspace}/chr${chr}

source /home/achan/git/mAFiA/mafia-venv/bin/activate

python3 /home/achan/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
--bam_file ${bam} \
--mod_file ${mod} \
--min_coverage 20 \
--out_dir ${output} \
--num_jobs 36
