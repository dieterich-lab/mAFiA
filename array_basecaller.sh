#!/usr/bin/env bash
#
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2GB
#SBATCH --verbose
#SBATCH --job-name=array_basecaller
#SBATCH --output=/home/achan/slurm/array_basecaller_%02a.out

set -e -u

source ${HOME}/git/renata/virtualenv/bin/activate

printf -v PART '%02d' "${SLURM_ARRAY_TASK_ID}"
LIST_FILENAMES=${WORKSPACE}/${FILENAME_PREFIX}${PART}
OUTPUT=${FASTA}${PART}
echo "Running basecaller on ${LIST_FILENAMES}"

python3 ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${FAST5_DIR} \
--list_filenames ${LIST_FILENAMES} \
--arch ${ARCH} \
--model ${MODEL} \
--decoder viterbi \
> ${OUTPUT}

rm ${LIST_FILENAMES}
