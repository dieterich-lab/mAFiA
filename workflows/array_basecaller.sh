#!/usr/bin/env bash
#
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8GB
#SBATCH --verbose
#SBATCH --job-name=array_basecaller
#SBATCH --output=/home/achan/slurm/array_basecaller_%A_%02a.out

set -e -u

source ${HOME}/git/renata/virtualenv/bin/activate

printf -v PART '%03d' "${SLURM_ARRAY_TASK_ID}"
LIST_FILENAMES=${WORKSPACE}/${FILENAME_PREFIX}${PART}
OUTPUT=${WORKSPACE}/part${PART}.fasta
echo "Running basecaller on ${LIST_FILENAMES}"

python3 ${HOME}/git/renata/basecall_viterbi.py \
--list_filenames ${LIST_FILENAMES} \
--arch ${ARCH} \
--model ${MODEL} \
--batchsize 2048 \
--decoder viterbi \
> ${OUTPUT}

rm ${LIST_FILENAMES}