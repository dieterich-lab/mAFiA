#!/usr/bin/env bash
#SBATCH --partition=long
#SBATCH --cpus-per-task=40
#SBATCH --mem=250GB
#SBATCH --verbose
#SBATCH --job-name=minimap_mESC_Mettl3_KO_merged
#SBATCH --output=/home/achan/slurm/minimap_mESC_Mettl3_KO_merged.out

#DATASET=mESC_WT_DMSO_merged
DATASET=mESC_Mettl3_KO_merged
#DATASET=mESC_WT_STM_merged


WORKSPACE=/prj/TRR319_RMaP/Project_B01/Adrian/${DATASET}
REF_GENOME=/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa
SAM_GENOME=${WORKSPACE}/genome_mapped.sam
BAM_GENOME=${WORKSPACE}/genome_filtered_q50.bam


if test -f ${WORKSPACE}/basecall_merged.fasta
then
  cat ${WORKSPACE}/part*/rodan.fasta >> ${WORKSPACE}/basecall_merged.fasta
else
  cat ${WORKSPACE}/part*/rodan.fasta > ${WORKSPACE}/basecall_merged.fasta
fi

module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF_GENOME} ${WORKSPACE}/basecall_merged.fasta > ${SAM_GENOME}

samtools view -bST ${REF_GENOME} -q50 ${SAM_GENOME} | samtools sort - > ${BAM_GENOME}
samtools index ${BAM_GENOME}

echo "minimap finsihed. Now performing QC..."

source ${HOME}/git/mAFiA/mafia-venv/bin/activate
samtools flagstats ${SAM_GENOME} > genome_qc.txt
${HOME}/git/renata/accuracy.py ${SAM_GENOME} ${REF_GENOME} >> ${WORKSPACE}/genome_qc.txt

less ${WORKSPACE}/genome_qc.txt
