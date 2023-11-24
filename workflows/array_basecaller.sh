#!/usr/bin/env bash
#
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --verbose
#SBATCH --job-name=array_basecaller
#SBATCH --output=/home/achan/slurm/array_basecaller_%A_%02a.out

set -e -u

#source ${HOME}/git/renata/virtualenv/bin/activate
source ${HOME}/git/mAFiA/mafia-venv/bin/activate

printf -v PART '%03d' "${SLURM_ARRAY_TASK_ID}"
LIST_FILENAMES=${WORKSPACE}/${FILENAME_PREFIX}${PART}
OUTDIR=${WORKSPACE}/part${PART}
echo "Running basecaller on ${LIST_FILENAMES}"

#python3 ${HOME}/git/renata/basecall_viterbi.py \
python3 ${HOME}/git/mAFiA_dev/RODAN/basecall.py \
--list_filenames ${LIST_FILENAMES} \
--model ${MODEL} \
--batchsize 4096 \
--outdir

rm ${LIST_FILENAMES}
