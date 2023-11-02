#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --verbose
#SBATCH --job-name=prep_mouse_heart
#SBATCH --output=/home/achan/slurm/prep_mouse_heart_%A.out

shopt -s globstar

ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

########################################################################################################################
#DATASET=40-26
DATASET=40-27
#DATASET=40-28
#DATASET=40-29
#DATASET=40-30
#DATASET=40-31
#DATASET=40-32
#DATASET=40-33
#DATASET=40-34

WORKSPACE=/prj/Dewenter_TAC_Backs_lab/achan/${DATASET}
DATA_DIR=/prj/Dewenter_TAC_Backs_lab/raw_data/Nanopore_dRNA/Cologne
########################################################################################################################
#mkdir -p ${WORKSPACE}
#cd ${WORKSPACE}

#if test -f basecall_merged.fasta
#then
#  cp basecall_merged.fasta basecall_merged.fasta.backup
#  mkdir _old
#  mv basecall_merged.fasta.backup filtered_q50.bam filtered_q50.bam.bai mapped.sam _old
#fi

########################################################################################################################
### basecall large number of reads #####################################################################################
########################################################################################################################
source ${HOME}/git/renata/virtualenv/bin/activate
#
#FAST5_DIR=${DATA_DIR}/${DATASET}/*/fast5_pass
#
#FILENAME_PREFIX=fast5_paths_part
#ls -1 ${FAST5_DIR}/*.fast5 > ${WORKSPACE}/fast5_paths_all
#split -a3 -l10 -d ${WORKSPACE}/fast5_paths_all ${WORKSPACE}/${FILENAME_PREFIX}
#
#NUM_ARRAYS=""
#for f in ${WORKSPACE}/fast5_paths_part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
#NUM_ARRAYS=${NUM_ARRAYS%,*}
#sbatch --array=${NUM_ARRAYS} --export=ALL,WORKSPACE=${WORKSPACE},FILENAME_PREFIX=${FILENAME_PREFIX},ARCH=${ARCH},MODEL=${MODEL} ${HOME}/git/mAFiA_dev/workflows/array_basecaller.sh

#for f in ${WORKSPACE}/part*.fasta; do echo $f; grep '>' $f | wc -l; done

if test -f basecall_merged.fasta
then
  cat ${WORKSPACE}/part*.fasta >> ${WORKSPACE}/basecall_merged.fasta
else
  cat ${WORKSPACE}/part*.fasta > ${WORKSPACE}/basecall_merged.fasta
fi

########################################################################################################################
#### align to genome ###################################################################################################
########################################################################################################################
#REF_GENOME=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa
REF_GENOME=/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa
SAM_GENOME=${WORKSPACE}/genome_mapped.sam
BAM_GENOME=${WORKSPACE}/genome_filtered_q50.bam

module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF_GENOME} ${WORKSPACE}/basecall_merged.fasta > ${SAM_GENOME}

### check stats and accuracy ###
samtools flagstats ${SAM_GENOME} > genome_qc.txt
${HOME}/git/renata/accuracy.py ${SAM_GENOME} ${REF_GENOME} >> genome_qc.txt

#### Convert to BAM and index ###
samtools view -bST ${REF_GENOME} -q50 ${SAM_GENOME} | samtools sort - > ${BAM_GENOME}
samtools index ${BAM_GENOME}

########################################################################################################################
#### align to transcriptome ############################################################################################
########################################################################################################################
#REF_TRANSCRIPTOME=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.cdna.all.fa
#SAM_TRANSCRIPTOME=${WORKSPACE}/transcriptome_mapped.sam
#BAM_TRANSCRIPTOME=${WORKSPACE}/transcriptome_mapped.bam
#
#module purge
#module load minimap2
##minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF_TRANSCRIPTOME} ${WORKSPACE}/basecall_merged.fasta > ${SAM_TRANSCRIPTOME}
#minimap2 --secondary=no -ax map-ont -k14 -t 36 --cs ${REF_TRANSCRIPTOME} ${WORKSPACE}/basecall_merged.fasta > ${SAM_TRANSCRIPTOME}
#### check stats and accuracy ###
#samtools flagstats ${SAM_TRANSCRIPTOME} > transcriptome_qc.txt
##${HOME}/git/renata/accuracy.py ${SAM_TRANSCRIPTOME} ${REF_TRANSCRIPTOME} >> transcriptome_qc.txt
#
##### Convert to BAM and index ###
#samtools view -bST ${REF_TRANSCRIPTOME} ${SAM_TRANSCRIPTOME} | samtools sort - > ${BAM_TRANSCRIPTOME}
#samtools index ${BAM_TRANSCRIPTOME}

########################################################################################################################
### clean up ###########################################################################################################
########################################################################################################################
#rm ${WORKSPACE}/part*.fasta
#rm ${WORKSPACE}/fast5_paths_all