#!/usr/bin/env bash

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian

ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
FILTER_SCORE=70

########################################################################################################################
### WUE batches 1-2 ####################################################################################################
########################################################################################################################
HOMOPOLYMER=0

ORIG=WUE
RUN=batch2
A_m6A=A

########################################################################################################################
### ISA runs 1-3 #######################################################################################################
########################################################################################################################
#HOMOPOLYMER=1
#
#ORIG=ISA
#RUN=run3_2
#MOD=A
#MOD=m6A

########################################################################################################################
### Isabel hetero ######################################################################################################
########################################################################################################################
#REF=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/ISA_M4_M5.fasta
#HOMOPOLYMER=0
#
#DATASET=RL_M4_M5
##DATASET=RL_M4_M5star
##DATASET=RL_M4star_M5
##DATASET=RL_M4star_M5star
#
#FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${DATASET}/fast5

########################################################################################################################

REF=${PRJ_DIR}/reference/${ORIG}_oligo_ref_${RUN}.fasta
DATASET=${ORIG}_${RUN}_${MOD}

WORKSPACE=${PRJ_DIR}/${DATASET}
cd ${WORKSPACE}

FAST5_DIR=${WORKSPACE}/fast5
FASTA=${WORKSPACE}/renata.fasta
SAM=${WORKSPACE}/spomelette.sam
LIGATION_REF=${WORKSPACE}/ligation_ref.fasta

#### basecall with Rodan IVT ###
#deactivate
source ${HOME}/git/renata/virtualenv/bin/activate

echo "Basecalling ${FAST5_DIR}"
srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
python3 -u ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${FAST5_DIR} \
--arch ${ARCH} \
--model ${MODEL} \
--batchsize 2048 \
--decoder viterbi \
> ${FASTA}

#### align with minimap ###
#module purge
#module load minimap2
#minimap2 --secondary=no -ax map-ont -t 36 --cs ${REF} ${FASTA} > ${SAM}

### align with spomelette ###
echo "Basecalling finished. Now aligning reads to ${REF}"
python3 -u ${HOME}/git/MAFIA/spanish_omelette_alignment.py \
--ref_file ${REF} \
--query_file ${FASTA} \
--recon_ref_file ${LIGATION_REF} \
--sam_file ${SAM} \
--homopolymer ${HOMOPOLYMER} \
--write_cs

### filter by quality ###
echo "Filtering and converting ${SAM}"
FILTERED_SAM=${SAM/.sam/_q${FILTER_SCORE}.sam}
samtools view -h -q${FILTER_SCORE} ${SAM} > ${FILTERED_SAM}

### check read num and accuracy ###
samtools flagstats ${FILTERED_SAM}
${HOME}/git/renata/accuracy.py ${FILTERED_SAM} ${LIGATION_REF}


#### Convert to BAM ###
echo "Converting ${BAM}"
BAM=${FILTERED_SAM//.sam/.bam}
samtools view -bST ${LIGATION_REF} ${FILTERED_SAM} | samtools sort - > ${BAM}
samtools index ${BAM}
