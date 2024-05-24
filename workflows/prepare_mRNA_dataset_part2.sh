#!/usr/bin/env bash
#SBATCH --partition=medium
#SBATCH --cpus-per-task=40
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=prepare_mRNA_dataset_part2
#SBATCH --output=/home/achan/slurm/prepare_mRNA_dataset_part2_%A.out

shopt -s globstar

cd ${WORKSPACE}
if [ ! -f "${WORKSPACE}"/basecall_merged.fasta ]; then
  cat "${WORKSPACE}"/part*/rodan.fasta > "${WORKSPACE}"/basecall_merged.fasta
fi

########################################################################################################################
#### align to genome ###################################################################################################
########################################################################################################################
sam_genome=${WORKSPACE}/genome_mapped.sam
bam_genome=${WORKSPACE}/genome_filtered_q50.bam

echo "Mapping to ${REF_GENOME}..."

module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF_GENOME} ${WORKSPACE}/basecall_merged.fasta > ${sam_genome}

### check stats and accuracy ###
echo "Quality control on ${sam_genome}..."

samtools flagstats ${sam_genome} > genome_qc.txt
echo >> genome_qc.txt
${HOME}/git/renata/accuracy.py ${sam_genome} ${REF_GENOME} >> genome_qc.txt

#### Convert to BAM and index ###
samtools view -bST ${REF_GENOME} -q50 ${sam_genome} | samtools sort - > ${bam_genome}
samtools index ${bam_genome}

#########################################################################################################################
#### split by chromosome ################################################################################################
#########################################################################################################################
echo "Split by chromosome..."

module load ont-fast5-api/4.1.1_deb12

# shellcheck disable=SC2203
if [[ $(basename "${REF_GENOME}") == *"GRCm"* ]]
then
  chrs=$(echo {1..19} X)
else
  chrs=$(echo {1..22} X)
fi

for chr in ${chrs}
do
  mkdir -p chr${chr}
  samtools view -h genome_filtered_q50.bam ${chr} | samtools sort - > chr${chr}/sorted.chr${chr}.bam
  samtools index chr${chr}/sorted.chr${chr}.bam
done

for chr in ${chrs}
do
  samtools view chr${chr}/sorted.chr${chr}.bam | cut -f1 > chr${chr}/read_ids.txt
  mkdir -p chr${chr}/fast5
#  sbatch -c 40 --mem 80GB \
  fast5_subset -t 36 -i fast5 -s chr${chr}/fast5 -l chr${chr}/read_ids.txt
done

#########################################################################################################################
#### clean up ###########################################################################################################
#########################################################################################################################
if [ -f "${WORKSPACE}"/basecall_merged.fasta ]; then
  rm -rf ${WORKSPACE}/part*
  rm ${WORKSPACE}/fast5_paths_all
fi

echo "Finished on ${WORKSPACE}"
