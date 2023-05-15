#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=MAFIA_SWARM_test_binary_classifier_mRNA
#SBATCH --output=/home/achan/slurm/MAFIA_SWARM_test_binary_classifier_mRNA_%A_%02a.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

set -e -u -f

printf -v PART '%02d' "${SLURM_ARRAY_TASK_ID}"

python3 ${HOME}/git/MAFIA/mRNA_test_binary_classifier.py \
--test_bam_file ${BAM} \
--test_fast5_dir ${FAST5_DIR} \
--ref_file ${REF} \
--mod_file ${MOD_FILE}.part${PART} \
--max_num_reads 1000 \
--min_coverage 50 \
--backbone_model_path ${BACKBONE_MODEL} \
--extraction_layer ${EXTRACTION_LAYER} \
--classifier logistic_regression \
--classifier_model_dir ${CLASSIFIER_MODEL_DIR} \
--outfile ${OUTFILE}.part${PART} \
--output_mod_probs
