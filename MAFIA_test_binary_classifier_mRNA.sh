#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --nodelist=gpu-g3-1
#SBATCH --mem=120GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_test_binary_classifier_mRNA_%A.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian

#####################################################################################################################################
DATASET=HEK293A_WT
FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all

#DATASET=HEK293_IVT
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293_IVT_2/fast5_pass

#####################################################################################################################################
#DATASET=HEK293T-WT-0-rep2
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-0-rep2/fast5_pass

#DATASET=HEK293T-WT-25-rep1
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-25-rep1/fast5_pass

#DATASET=HEK293T-WT-50-rep2
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-50-rep2/fast5_pass

#DATASET=HEK293T-WT-50-rep3
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-50-rep3/fast5

#DATASET=HEK293T-WT-75-rep4
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-75-rep4/fast5

#DATASET=HEK293T-WT-100-rep1
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-100-rep1/fast5_pass

#####################################################################################################################################
#DATASET=HEK293T-Mettl3-KO-rep3
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-Mettl3-KO-rep3/fast5

#DATASET=HEK293T-WT-rep3
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-rep3/fast5

#####################################################################################################################################

REF=${HOME}/Data/genomes/GRCh38_96.fa
BAM=${WORKSPACE}/mapping/${DATASET}.bam.sorted
MOD_FILE=${HOME}/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXTRACTION_LAYER=convlayers.conv21
TRAINING_DATA=random_ligation_A_m6A_MaxAbs
CLASSIFIER_MODEL_DIR=${WORKSPACE}/MAFIA_classifiers/${TRAINING_DATA}
MOD_PROB_THRESH=0.5
#OUTFILE=${WORKSPACE}/results/res_${DATASET}_${TRAINING_DATA}_modProbThresh${MOD_PROB_THRESH}.tsv
OUTFILE=${WORKSPACE}/results/res_${DATASET}_${TRAINING_DATA}_modProbPerRead.tsv

set -e -u

cd ${HOME}/git/MAFIA

python3 ${HOME}/git/MAFIA/mRNA_test_binary_classifier.py \
--test_bam_file ${BAM} \
--test_fast5_dir ${FAST5_DIR} \
--ref_file ${REF} \
--mod_file ${MOD_FILE} \
--max_num_reads 1000 \
--min_coverage 50 \
--backbone_model_path ${BACKBONE_MODEL} \
--extraction_layer ${EXTRACTION_LAYER} \
--classifier logistic_regression \
--classifier_model_dir ${CLASSIFIER_MODEL_DIR} \
--mod_prob_thres ${MOD_PROB_THRESH} \
--outfile ${OUTFILE} \
--output_mod_probs
