#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup_Fhl1
#SBATCH --output=/home/achan/slurm/pileup_Fhl1_%A.out

workspace="/beegfs/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/TAC_day56"
out_dir=${workspace}

cd ${workspace}

chr="X"
gene="Fhl1"

tx_outer="ENSMUST00000023854"
tx_inner="ENSMUST00000114772"

range_outer="56754335-56754391"
range_inner="56779723-56779794"

samtools view -f 0 chrALL.mAFiA.reads.bam ${chr}:${range_inner} -o ${gene}.mAFiA.reads.bam
samtools index ${gene}.mAFiA.reads.bam

samtools view ${gene}.mAFiA.reads.bam ${chr}:${range_outer} -o ${gene}.${tx_outer}.mAFiA.reads.bam
samtools index ${gene}.${tx_outer}.mAFiA.reads.bam

samtools view ${gene}.${tx_outer}.mAFiA.reads.bam | cut -f1 > ids_${tx_outer}
samtools view -N ids_${tx_outer} ${gene}.mAFiA.reads.bam -U ${gene}.${tx_inner}.mAFiA.reads.bam
samtools index ${gene}.${tx_inner}.mAFiA.reads.bam

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.chr${chr}.bed

for tx in ${tx_inner} ${tx_outer}
do
  bam=${workspace}/${gene}.${tx}.mAFiA.reads.bam
  python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
    --bam_file ${bam} \
    --mod_file ${mod} \
    --min_coverage 1 \
    --out_dir ${out_dir} \
    --out_filename ${gene}.${tx}.bed \
    --num_jobs 36
done

echo "Finished"
