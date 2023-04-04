#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --mem=60GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_test_binary_classifier_tRNA_%A.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

WORKSPACE=/beegfs/prj/tRNA_Berlin/newBatchDec2022_Spombe

#####################################################################################################################################
#DATASET=AEP1_T1
DATASET=AEP1_T2
#DATASET=AEP1_T3
#DATASET=AEP1_Q_T1
#DATASET=AEP1_Q_T2
#DATASET=AEP1_Q_T3

FAST5_DIR=${WORKSPACE}/agg_fast5_pass/${DATASET}

#####################################################################################################################################

REF=${WORKSPACE}/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta
BAM=${WORKSPACE}/achan/mapping/${DATASET}.bam
MOD_FILE=${WORKSPACE}/S_pombe_Q_only.tsv
BACKBON_MODEL=${HOME}/pytorch_models/tRNA_IVT/tRNA_IVT-epoch29.torch
EXTRACTION_LAYER=convlayers.conv21
CLASSIFIER_MODEL_DIR=${WORKSPACE}/achan/rev_filtered_classifier_models
MOD_PROB_THRESH=0.0
OUTFILE=${WORKSPACE}/results/res_${DATASET}_modProbThresh${MOD_PROB_THRESH}.tsv

set -e -u

cd ${HOME}/git/MAFIA

python3 ${HOME}/git/MAFIA/tRNA_test_binary_classifier.py \
--test_bam_file ${BAM} \
--test_fast5_dir ${FAST5_DIR} \
--ref_file ${REF} \
--mod_file ${MOD_FILE} \
--max_num_reads -1 \
--min_coverage 10 \
--backbone_model_path ${BACKBON_MODEL} \
--extraction_layer ${EXTRACTION_LAYER} \
--feature_width 0 \
--classifier logistic_regression \
--classifier_model_dir ${CLASSIFIER_MODEL_DIR} \
--mod_prob_thres ${MOD_PROB_THRESH} \
--outfile ${OUTFILE}
