#!/usr/bin/env bash

########################################################################################################################
### Train data #########################################################################################################
########################################################################################################################
#TRAIN_DATASET=ISA_runs1-3
#TRAIN_DATASET=WUE_batches1-2
#TRAIN_DATASET=ISA-WUE
#TRAIN_DATASET=ISA_mixes1-4
TRAIN_DATASET=DRACH

########################################################################################################################
### Test data ##########################################################################################################
########################################################################################################################
# new HEK293
TEST_DATASET=100_WT_0_IVT
#TEST_DATASET=75_WT_25_IVT
#TEST_DATASET=50_WT_50_IVT
#TEST_DATASET=25_WT_75_IVT
#TEST_DATASET=0_WT_100_IVT

# P2
#TEST_DATASET=P2_WT
########################################################################################################################
PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A
TEST_DIR=${PRJ_DIR}/HEK293/${TEST_DATASET}
FAST5_DIR=${TEST_DIR}/fast5
BAM=${TEST_DIR}/genome_filtered_q50.bam
OUTDIR=${PRJ_DIR}/results/train_${TRAIN_DATASET}_test_${TEST_DATASET}

#FAST5_DIR=${TEST_DIR}/fast5_chrX
#BAM=${TEST_DIR}/chrX_q50.bam
#OUTDIR=${PRJ_DIR}/results/train_${TRAIN_DATASET}_test_${TEST_DATASET}_chrX
#OUTFILE=${OUTDIR}/res_train_${TRAIN_DATASET}_test_${TEST_DATASET}_chrX.tsv

REF=${HOME}/Data/genomes/GRCh38_96.fa
MOD_FILE=${HOME}/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv
#MOD_FILE=${HOME}/Data/GLORI/bed_files/GLORI_chrX.tsv
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
CLASSIFIER_MODEL_DIR=${PRJ_DIR}/MAFIA_classifiers/${TRAIN_DATASET}
########################################################################################################################
########################################################################################################################
mkdir -p ${OUTDIR}

NUM_ARRAYS=""
for f in ${MOD_FILE}.part*; do ff=${f##*part}; ff=$(echo $ff | sed 's/^0*//'); NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}

sbatch --array=${NUM_ARRAYS} \
--export=ALL,\
FAST5_DIR=${FAST5_DIR},\
REF=${REF},\
BAM=${BAM},\
MOD_FILE=${MOD_FILE},\
BACKBONE_MODEL=${BACKBONE_MODEL},\
CLASSIFIER_MODEL_DIR=${CLASSIFIER_MODEL_DIR},\
OUTDIR=${OUTDIR} \
${HOME}/git/mAFiA_dev/workflows/array_test_MAFIA_mRNA.sh

### concat output ###
#cp ${OUTFILE}.part00 ${OUTFILE}.merged
#for num in ${NUM_ARRAYS//,/ }
#do
#  if [ $num != '00' ]
#  then awk NR\>1 $OUTFILE.part$num >> ${OUTFILE}.merged
#  fi
#  done
