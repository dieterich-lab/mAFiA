#!/usr/bin/env bash

for TRAIN_DATASET in ISA WUE ISA-WUE
do
  for TEST_DATASET in ISA WUE ISA-WUE
  do
    echo "Train - ${TRAIN_DATASET}; Test - ${TEST_DATASET}"
    sbatch --export=ALL,TRAIN_DATASET=${TRAIN_DATASET},TEST_DATASET=${TEST_DATASET} ${HOME}/git/MAFIA/workflows/test_MAFIA_oligo.sh
  done
done
