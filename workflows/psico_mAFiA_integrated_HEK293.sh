#!/usr/bin/env bash

ds="WT_P2"

jid=$(sbatch --parsable --export=ALL,ds="${ds}" --array=1-23 "${HOME}"/git/mAFiA_dev/workflows/mAFiA_process_reads_HEK293.sh)

sbatch --dependency=afterok:"${jid}" --export=ALL,ds=${ds} --array=1-23 "${HOME}"/git/mAFiA_dev/workflows/mAFiA_pileup_HEK293.sh
