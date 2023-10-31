#!/usr/bin/env bash

#DB=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/oligo/db_ISA_run2.tsv
#DB=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/oligo/db_ISA_run4.tsv
#DB=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/oligo/db_WUE_batch3.tsv
#DB=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/db_ISA_mixes1-4.tsv
DB=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/oligo/db_ISA_mixes1-4_m6A.tsv

{
  read -r
  while read -r o r m h l
  do
    echo "ORIG=$o, RUN=$r, MOD=$m, HOMOPOLYMER=$h, LOC=$l"
    bash ${HOME}/git/mAFiA_dev/workflows/prepare_single_oligo_dataset.sh -o "$o" -r "$r" -m "$m" -h "$h" -l "$l" &
  done
} < $DB
