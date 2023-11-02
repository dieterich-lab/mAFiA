#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=test_mAFiA_mouse
#SBATCH --output=/home/achan/slurm/test_mAFiA_mouse_%A.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

set -e -f

DATASET=40-26

FAST5_DIR=/prj/Dewenter_TAC_Backs_lab/raw_data/Nanopore_dRNA/Cologne/${DATASET}/*/fast5_pass
BAM=/prj/Dewenter_TAC_Backs_lab/achan/${DATASET}/genome_filtered_q50.bam
OUTDIR=/prj/Dewenter_TAC_Backs_lab/achan/${DATASET}/mAFiA_predictions

REF=/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa
MOD_FILE=/prj/Dewenter_TAC_Backs_lab/achan/bed_files/m6Anet_sites.bed
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
CLASSIFIER_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/DRACH

python3 -u ${HOME}/git/mAFiA_dev/mAFiA/test_mAFiA.py \
--bam_file ${BAM} \
--fast5_dir ${FAST5_DIR} \
--ref_file ${REF} \
--mod_file ${MOD_FILE} \
--max_num_reads 1000 \
--min_coverage 10 \
--backbone_model_path ${BACKBONE_MODEL} \
--classifier_model_dir ${CLASSIFIER_DIR} \
--out_dir ${OUTDIR} \
--output_mod_probs
