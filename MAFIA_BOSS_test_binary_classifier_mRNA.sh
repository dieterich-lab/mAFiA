#!/usr/bin/env bash

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
REF=${HOME}/Data/genomes/GRCh38_96.fa
MOD_FILE=${HOME}/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXTRACTION_LAYER=convlayers.conv21

### new HEK293 ######################################################################################################################
#TEST_DATASET=100_WT_0_IVT
#TEST_DATASET=75_WT_25_IVT
#TEST_DATASET=50_WT_50_IVT
#TEST_DATASET=25_WT_75_IVT
TEST_DATASET=0_WT_100_IVT
#####################################################################################################################################

WORKSPACE=${PRJ_DIR}/${TEST_DATASET}
FAST5_DIR=${WORKSPACE}/fast5
BAM=${WORKSPACE}/filtered_q50.bam
TRAIN_DATASET=WUE_combined
CLASSIFIER_MODEL_DIR=${PRJ_DIR}/MAFIA_classifiers/${TRAIN_DATASET}
OUTDIR=${PRJ_DIR}/results/train_${TRAIN_DATASET}_test_${TEST_DATASET}
OUTFILE=${OUTDIR}/res_train_${TRAIN_DATASET}_test_${TEST_DATASET}.tsv

mkdir -p ${OUTDIR}

NUM_ARRAYS=""
for f in ${MOD_FILE}.part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}

sbatch --array=${NUM_ARRAYS} \
--export=ALL,\
WORKSPACE=${WORKSPACE},\
FAST5_DIR=${FAST5_DIR},\
REF=${REF},\
BAM=${BAM},\
MOD_FILE=${MOD_FILE},\
BACKBONE_MODEL=${BACKBONE_MODEL},\
EXTRACTION_LAYER=${EXTRACTION_LAYER},\
CLASSIFIER_MODEL_DIR=${CLASSIFIER_MODEL_DIR},\
OUTFILE=${OUTFILE} \
${HOME}/git/MAFIA/MAFIA_SWARM_test_binary_classifier_mRNA.sh

### concat output ###
cp ${OUTFILE}.part00 ${OUTFILE}.merged
for num in ${NUM_ARRAYS//,/ }
do
  if [ $num != '00' ]
  then awk NR\>1 $OUTFILE.part$num >> ${OUTFILE}.merged
  fi
  done

#rm ${OUTFILE}.part*