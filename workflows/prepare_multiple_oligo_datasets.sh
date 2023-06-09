#!/usr/bin/env bash

DB=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/oligo/db_ISA_run2.tsv

{
  read -r
  while read -r o r m h l
  do
    echo "ORIG=$o, RUN=$r, MOD=$m, HOMOPOLYMER=$h, LOC=$l"
    bash ${HOME}/git/MAFIA/workflows/prepare_single_oligo_dataset.sh -o "$o" -r "$r" -m "$m" -h "$h" -l "$l" &
  done
} < $DB