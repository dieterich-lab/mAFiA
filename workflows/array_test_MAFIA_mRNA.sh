#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=array_test_MAFIA_mRNA
#SBATCH --output=/home/achan/slurm/array_test_MAFIA_mRNA_%A_%02a.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

set -e -f

printf -v PART '%02d' "${SLURM_ARRAY_TASK_ID}"

python3 -u ${HOME}/git/MAFIA/mRNA_test_binary_classifier.py \
--test_bam_file ${BAM} \
--test_fast5_dir ${FAST5_DIR} \
--ref_file ${REF} \
--mod_file ${MOD_FILE}.part${PART} \
--max_num_reads 1000 \
--min_coverage 50 \
--backbone_model_path ${BACKBONE_MODEL} \
--classifier_model_dir ${CLASSIFIER_MODEL_DIR} \
--outfile ${OUTFILE}.part${PART} \
--output_mod_probs
