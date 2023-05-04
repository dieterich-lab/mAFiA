#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --mem=120GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_train_binary_classifier_oligo_%A.out

eval "$(conda shell.bash hook)"
#eval "$(/home/achan/miniconda3/condabin/conda shell.bash hook)"
conda activate MAFIA
set -e -u -f
cd ${HOME}/git/MAFIA

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian

############################################################################################
### feature settings #######################################################################
############################################################################################
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXT_LAYER=convlayers.conv21
CLASSIFIER=logistic_regression
SCALER=MaxAbs

############################################################################################
### source data ############################################################################
############################################################################################
#### Claudia first 3 splint ###
#DATASET=20230221_WUE_splint_lig
#UNM_BAM=${WORKSPACE}/WUE_splint_lig_A_RTA/mapped.bam
#UNM_FAST5="/prj/TRR319_RMaP/Project_BaseCalling/Isabel/${DATASET}/WUE_splint_lig_A_RTA/*/fast5_*"
#MOD_BAM=${WORKSPACE}/WUE_splint_lig_m6A_RTA/mapped.bam
#MOD_FAST5="/prj/TRR319_RMaP/Project_BaseCalling/Isabel/${DATASET}/WUE_splint_lig_m6A_RTA/*/fast5_*"
#REF=${WORKSPACE}/reference/splint_variations_max_blocks_7.fasta
#OUTDIR=${WORKSPACE}/MAFIA_classifiers/${DATASET}_${CLASSIFIER}_${SCALER}

#### Claudia batch 2 ###
#DATASET=20230502_WUE_splint_batch2
#UNM_BAM=${WORKSPACE}/WUE_splint_batch2_A_RTA/mapped.bam
#UNM_FAST5="/prj/TRR319_RMaP/Project_BaseCalling/Isabel/${DATASET}/WUE_splint_batch2_A_RTA/*/fast5_*"
#MOD_BAM=${WORKSPACE}/WUE_splint_batch2_m6A_RTA/mapped.bam
#MOD_FAST5="/prj/TRR319_RMaP/Project_BaseCalling/Isabel/${DATASET}/WUE_splint_batch2_m6A_RTA/*/fast5_*"
#REF=${WORKSPACE}/reference/WUE_batch2_max_blocks_7.fasta
#OUTDIR=${WORKSPACE}/MAFIA_classifiers/${DATASET}_${CLASSIFIER}_${SCALER}

### Isabel random 6 ###
#UNM_BAM=${WORKSPACE}/mapping/RL_RG1-6_A_RTA.bam
#UNM_FAST5=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230418_Random_Ligation_A_m6A/RL_RG1-6_A_RTA/20230418_1325_X1_AOL616_885f620d/fast5
#MOD_BAM=${WORKSPACE}/mapping/RL_RG7-12_m6A_RTA.bam
#MOD_FAST5=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230418_Random_Ligation_A_m6A/RL_RG7-12_m6A_RTA/20230418_1325_X2_AOC149_8138c168/fast5
#REF=${WORKSPACE}/reference/top6_random_permutation_max_blocks_5.fasta
#OUTDIR=${WORKSPACE}/MAFIA_classifiers/random_ligation_A_m6A_${SCALER}_enforceMotif

### Isabel RL Mix 1 & 3 ###
UNM_BAM=${WORKSPACE}/RL_Mix1_A_RTA/mapped.bam
UNM_FAST5=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix1_A_RTA/*
MOD_BAM=${WORKSPACE}/RL_Mix3_m6A_RTA/mapped.bam
MOD_FAST5=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix3_m6A_RTA/*
REF=${WORKSPACE}/reference/RL_Mix1_Mix3_blocks7.fasta
OUTDIR=${WORKSPACE}/MAFIA_classifiers/RL_run2_Mix1_Mix3_${SCALER}

python3 oligo_train_binary_classifier.py \
--unm_bam_file ${UNM_BAM} \
--unm_fast5_dir ${UNM_FAST5} \
--mod_bam_file ${MOD_BAM} \
--mod_fast5_dir ${MOD_FAST5} \
--ref_file ${REF} \
--backbone_model_path ${BACKBONE_MODEL} \
--extraction_layer ${EXT_LAYER} \
--scaler ${SCALER} \
--min_coverage 10 \
--classifier ${CLASSIFIER} \
--classifier_model_dir ${OUTDIR}
