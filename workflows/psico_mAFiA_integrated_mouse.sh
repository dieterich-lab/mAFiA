#!/usr/bin/env bash

ds="CM/IVT"

jid=$(sbatch --parsable --export=ALL,ds="${ds}" --array=1-19,23,24 "${HOME}"/git/mAFiA_dev/workflows/mAFiA_process_reads_mouse.sh)

sbatch --dependency=afterok:"${jid}" --export=ALL,ds=${ds} --array=1-19,23,24 "${HOME}"/git/mAFiA_dev/workflows/mAFiA_pileup_mouse.sh
