#!/usr/bin/env bash
shopt -s globstar

#ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

### new HEK293 #########################################################################################################
#DATASET=0_WT_100_IVT
#DATASET=25_WT_75_IVT
#DATASET=50_WT_50_IVT
#DATASET=75_WT_25_IVT
#DATASET=100_WT_0_IVT
#DATASET=P2_WT
#DATASET=Mettl3-KO

#WORKSPACE=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/${DATASET}

########################################################################################################################
#DATASET=A1_WT_CD
#DATASET=A2_WT_CD
#WORKSPACE=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/Federica_Accornero/${DATASET}

########################################################################################################################

#DATASET=JK_HEK293_DMSO_1_2_RTA
#DATASET=JK_HEK293_DMSO_3_4_RTA
#DATASET=JK_HEK293_STM2457_5_6_RTA
#DATASET=JK_HEK293_STM2457_7_8_RTA
#DATASET=JK_HEK293_DMSO_merged

#DATASET=mESC_WT_DMSO_merged
#DATASET=mESC_Mettl3_KO_merged
#DATASET=mESC_WT_STM_merged

#WORKSPACE=/prj/TRR319_RMaP/Project_B01/Adrian/${DATASET}

########################################################################################################################
#DATASET=40-26
#DATASET=40-27
#DATASET=40-28
#DATASET=40-29
#DATASET=40-30
#DATASET=40-31
#DATASET=40-32
#DATASET=40-33
#DATASET=40-34

#WORKSPACE=/prj/Dewenter_TAC_Backs_lab/achan/${DATASET}
#DATA_DIR=/prj/Dewenter_TAC_Backs_lab/raw_data/Nanopore_dRNA/Cologne
########################################################################################################################

#DATASET=col0
#DATASET=vir1
#WORKSPACE=/scratch/achan/Arabidopsis_thaliana/${DATASET}

########################################################################################################################

#DATASET=WT
#DATASET=IME4_KO
#WORKSPACE=/scratch/achan/Saccharomyces_cerevisiae/${DATASET}

########################################################################################################################
#ds=HEK_siCtrl_input_rep1
#ds=HEK_siMETTL3_input_rep1
#ds=HEK_siTRUB1_input_rep1

ds=HEK_siCtrl_input_rep2
#ds=HEK_siMETTL3_input_rep2
#ds=HEK_siTRUB1_input_rep2

WORKSPACE=/prj/TRR319_RMaP_BaseCalling/Adrian/NanoSPA/${ds}
########################################################################################################################

mkdir -p ${WORKSPACE}
cd ${WORKSPACE}

#if test -f basecall_merged.fasta
#then
#  cp basecall_merged.fasta basecall_merged.fasta.backup
#  mkdir _old
#  mv basecall_merged.fasta.backup filtered_q50.bam filtered_q50.bam.bai mapped.sam _old
#fi

########################################################################################################################
### basecall large number of reads #####################################################################################
########################################################################################################################
#source ${HOME}/git/renata/virtualenv/bin/activate
source ${HOME}/git/mAFiA/mafia-venv/bin/activate

#FAST5_DIR=${DATA_DIR}/${DATASET}/*/fast5_pass
FAST5_DIR=${WORKSPACE}/fast5

FILENAME_PREFIX=fast5_paths_part
ls -1 ${FAST5_DIR}/*.fast5 > ${WORKSPACE}/fast5_paths_all
split -a3 -l10 -d ${WORKSPACE}/fast5_paths_all ${WORKSPACE}/${FILENAME_PREFIX}

NUM_ARRAYS=""
for f in ${WORKSPACE}/fast5_paths_part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}
sbatch --array=${NUM_ARRAYS} --export=ALL,WORKSPACE=${WORKSPACE},FILENAME_PREFIX=${FILENAME_PREFIX},ARCH=${ARCH},MODEL=${MODEL} ${HOME}/git/mAFiA_dev/workflows/array_basecaller.sh

#for f in ${WORKSPACE}/part*.fasta; do echo $f; grep '>' $f | wc -l; done

if test -f basecall_merged.fasta
then
  cat ${WORKSPACE}/part*/rodan.fasta >> ${WORKSPACE}/basecall_merged.fasta
else
  cat ${WORKSPACE}/part*/rodan.fasta > ${WORKSPACE}/basecall_merged.fasta
fi

########################################################################################################################
#### align to genome ###################################################################################################
########################################################################################################################
REF_GENOME=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa
#REF_GENOME=/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa
#REF_GENOME='/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/reference/TAIR10_chr_all.fasta'
#REF_GENOME='/biodb/genomes/saccharomyces_cerevisiae/R64-1-1_96/R64-1-1_96.fa'
SAM_GENOME=${WORKSPACE}/genome_mapped.sam
BAM_GENOME=${WORKSPACE}/genome_filtered_q50.bam

module purge
module load minimap2
srun -c 40 --mem 120GB \
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF_GENOME} ${WORKSPACE}/basecall_merged.fasta > ${SAM_GENOME}

### check stats and accuracy ###
samtools flagstats ${SAM_GENOME} > genome_qc.txt
echo >> genome_qc.txt
${HOME}/git/renata/accuracy.py ${SAM_GENOME} ${REF_GENOME} >> genome_qc.txt

#### Convert to BAM and index ###
samtools view -bST ${REF_GENOME} -q50 ${SAM_GENOME} | samtools sort - > ${BAM_GENOME}
samtools index ${BAM_GENOME}

########################################################################################################################
### split by chromosome ################################################################################################
########################################################################################################################
module load ont-fast5-api
#for chr in {1..5}
#for chr in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito
#for chr in {1..22} X

for chr in {1..22} X
do
  mkdir chr${chr}
  samtools view -h genome_filtered_q50.bam ${chr} | samtools sort - > chr${chr}/sorted.chr${chr}.bam
  samtools index chr${chr}/sorted.chr${chr}.bam
done

for chr in {1..22} X
do
  samtools view chr${chr}/sorted.chr${chr}.bam | cut -f1 > chr${chr}/read_ids.txt
  mkdir chr${chr}/fast5
  srun -c 40 --mem 64GB -o ${HOME}/slurm/fast5_subset_chr${chr}.out -e ${HOME}/slurm/fast5_subset_chr${chr}.err \
  fast5_subset -t 36 -i fast5 -s chr${chr}/fast5 -l chr${chr}/read_ids.txt &
done

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
rm -rf ${WORKSPACE}/part*
rm ${WORKSPACE}/fast5_paths_all

########################################################################################################################
### filter and merge bam ###############################################################################################
########################################################################################################################

### filter ###
#for chr in {1..19} X
#do
#  mv chr${chr}/mAFiA.reads.bam chr${chr}/_mAFiA.reads.bam
#  samtools index chr${chr}/_mAFiA.reads.bam
#  samtools view -h chr${chr}/_mAFiA.reads.bam ${chr} | samtools sort - > chr${chr}/mAFiA.reads.bam
#done

### merge bam ###
#samtools merge -o chrALL.mAFiA.reads.bam chr*/mAFiA.reads.bam
#samtools index chrALL.mAFiA.reads.bam
#
#### merge bed ###
#cp chr1/mAFiA.sites.bed chrALL.mAFiA.sites.bed
#for chr in {2..22} X
#do
#  tail -n+2 chr${chr}/mAFiA.sites.bed >> chrALL.mAFiA.sites.bed
#done