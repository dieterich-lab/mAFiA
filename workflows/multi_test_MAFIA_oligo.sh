#!/usr/bin/env bash

for TRAIN_DATASET in ISA WUE ISA-WUE
do
  for TEST_DATASET in ISA_A ISA_m6A WUE_A WUE_m6A ISA-WUE_A ISA-WUE_m6A
  do
    echo "Train - ${TRAIN_DATASET}; Test - ${TEST_DATASET}"
    sbatch --export=ALL,TRAIN_DATASET=${TRAIN_DATASET},TEST_DATASET=${TEST_DATASET} ${HOME}/git/MAFIA/workflows/test_MAFIA_oligo.sh
  done
done
