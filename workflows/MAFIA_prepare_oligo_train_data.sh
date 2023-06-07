#!/usr/bin/env bash

ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
FILTER_SCORE=70

########################################################################################################################
### Würzburg ###########################################################################################################
########################################################################################################################
HOMOPOLYMER=0

WUE_BATCH=2
A_m6A=m6A

REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_oligo_ref_batch${WUE_BATCH}.fasta
TRAIN_DATASET=WUE_batch${WUE_BATCH}_${A_m6A}
FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}/fast5

########################################################################################################################
### Isabel run 1 #######################################################################################################
########################################################################################################################
#HOMOPOLYMER=1

#REF=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/ISA_run1.fasta
#
##TRAIN_DATASET=ISA_run1_A
#TRAIN_DATASET=ISA_run1_m6A
#
#FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}/fast5
#
########################################################################################################################
### Isabel run 2 #######################################################################################################
########################################################################################################################
#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/RL_Mix1_Mix3_blocks8.fasta
#TRAIN_DATASET=RL_Mix1_A
#TRAIN_DATASET=RL_Mix3_m6A

#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/RL_Mix2_Mix4_blocks8.fasta
#TRAIN_DATASET=RL_Mix2_A
#TRAIN_DATASET=RL_Mix4_m6A

#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/${TRAIN_DATASET}/*

########################################################################################################################
### Isabel hetero ######################################################################################################
########################################################################################################################
#REF=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/ISA_M4_M5.fasta
#HOMOPOLYMER=0
#
#TRAIN_DATASET=RL_M4_M5
##TRAIN_DATASET=RL_M4_M5star
##TRAIN_DATASET=RL_M4star_M5
##TRAIN_DATASET=RL_M4star_M5star
#
#FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}/fast5

########################################################################################################################

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}
mkdir -p ${WORKSPACE}
cd ${WORKSPACE}

FASTA=${WORKSPACE}/renata.fasta
SAM=${WORKSPACE}/spomelette.sam
LIGATION_REF=${WORKSPACE}/ligation_ref.fasta

#### basecall with Rodan IVT ###
#deactivate
source ${HOME}/git/renata/virtualenv/bin/activate

#srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
#python3 -u ${HOME}/git/renata/basecall_viterbi.py \
#--fast5dir ${FAST5_DIR} \
#--arch ${ARCH} \
#--model ${MODEL} \
#--batchsize 2048 \
#--decoder viterbi \
#> ${FASTA} &

#### align with minimap ###
#module purge
#module load minimap2
#minimap2 --secondary=no -ax map-ont -t 36 --cs ${REF} ${FASTA} > ${SAM}

### align with spomelette ###
python3 -u ${HOME}/git/MAFIA/spanish_omelette_alignment.py \
--ref_file ${REF} \
--query_file ${FASTA} \
--recon_ref_file ${LIGATION_REF} \
--sam_file ${SAM} \
--homopolymer ${HOMOPOLYMER} \
--write_cs

### filter by quality ###
FILTERED_SAM=${SAM/.sam/_q${FILTER_SCORE}.sam}
samtools view -h -q${FILTER_SCORE} ${SAM} > ${FILTERED_SAM}

### check read num and accuracy ###
samtools flagstats ${FILTERED_SAM}
${HOME}/git/renata/accuracy.py ${FILTERED_SAM} ${LIGATION_REF}


#### Convert to BAM ###
BAM=${FILTERED_SAM//.sam/.bam}
samtools view -bST ${LIGATION_REF} ${FILTERED_SAM} | samtools sort - > ${BAM}
samtools index ${BAM}
