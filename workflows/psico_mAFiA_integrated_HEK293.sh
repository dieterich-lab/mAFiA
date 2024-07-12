#!/usr/bin/env bash

ds="HEK293_TRUB1_OE/rep1"

jid=$(sbatch --parsable --export=ALL,ds="${ds}" --array=1-23 "${HOME}"/git/mAFiA_dev/workflows/mAFiA_process_reads_HEK293.sh)

sbatch --dependency=afterok:"${jid}" --export=ALL,ds=${ds} --array=1-23 "${HOME}"/git/mAFiA_dev/workflows/mAFiA_pileup_HEK293.sh
