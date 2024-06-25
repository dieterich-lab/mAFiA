#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup_bambu_Mettl3
#SBATCH --output=/home/achan/slurm/pileup_Rcan1.out

workspace="/beegfs/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/40-33"
out_dir=${workspace}

cd ${workspace}

chr="16"

samtools view -f 16 chrALL.mAFiA.reads.bam 16:92399896-92400077 -o Rcan1.mAFiA.reads.bam
samtools index Rcan1.mAFiA.reads.bam

samtools view Rcan1.mAFiA.reads.bam 16:92465833-92466146 -o Rcan1.exon2.mAFiA.reads.bam
samtools index Rcan1.exon2.mAFiA.reads.bam

samtools view Rcan1.exon2.mAFiA.reads.bam | cut -f1 > ids_exon2
samtools view -N ids_exon2 Rcan1.mAFiA.reads.bam -U Rcan1.exon4.mAFiA.reads.bam
samtools index Rcan1.exon4.mAFiA.reads.bam

source ${HOME}/git/mAFiA/mafia-venv/bin/activate
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102/m6A.psi.GRCh38_102.chr${chr}.bed

for exon in 2 4
do
  bam="${workspace}/Rcan1.exon${exon}.mAFiA.reads.bam"
  python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
    --bam_file ${bam} \
    --mod_file ${mod} \
    --min_coverage 1 \
    --out_dir ${out_dir} \
    --out_filename "Rcan1.exon${exon}.bed" \
    --num_jobs 36
done

echo "Finished"
