#!/usr/bin/env bash

#for TRAIN_DATASET in ISA WUE ISA-WUE
for TRAIN_DATASET in ISA-WUE
do
#  for TEST_DATASET in ISA_A ISA_m6A WUE_A WUE_m6A ISA-WUE_A ISA-WUE_m6A
#  for TEST_DATASET in ISA_run4_M4M5 ISA_run4_M4M5star ISA_run4_M4starM5 ISA_run4_M4starM5star
for TEST_DATASET in WUE_batch3_AB WUE_batch3_BA WUE_batch3_ABBA
  do
    echo "Train - ${TRAIN_DATASET}; Test - ${TEST_DATASET}"
    sbatch --export=ALL,TRAIN_DATASET=${TRAIN_DATASET},TEST_DATASET=${TEST_DATASET} ${HOME}/git/MAFIA/workflows/test_MAFIA_oligo.sh
  done
done
