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

ROOTDIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian

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
#### Claudia batch 1 ###
#BATCH=WUE_batch1
#UNM_BAM=${ROOTDIR}/${BATCH}_A/spomlette_q80.bam
#UNM_FAST5="/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_A_RTA/*/fast5_*"
#MOD_BAM=${ROOTDIR}/${BATCH}_m6A/spomlette_q70.bam
#MOD_FAST5="/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_m6A_RTA/*/fast5_*"
#REF=${ROOTDIR}/reference/${BATCH}_combined.fa
#OUTDIR=${ROOTDIR}/MAFIA_classifiers/${BATCH}_spomlette_${CLASSIFIER}_${SCALER}

#### Claudia batch 2 ###
TRAIN_DATASET=WUE_batch2
UNM_BAM=${ROOTDIR}/${TRAIN_DATASET}_A/spomelette_q70.bam
UNM_FAST5=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_batch2_A/fast5
MOD_BAM=${ROOTDIR}/${TRAIN_DATASET}_m6A/spomelette_q70.bam
MOD_FAST5=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_batch2_m6A/fast5
REF=${ROOTDIR}/reference/${TRAIN_DATASET}_combined.fa
OUTDIR=${ROOTDIR}/MAFIA_classifiers/${TRAIN_DATASET}_spomlette_${CLASSIFIER}_${SCALER}

### Isabel RL top 6 ###
#TRAIN_DATASET=RL_top6
#UNM_BAM=${ROOTDIR}/${TRAIN_DATASET}_A/spomlette_q70.bam
#UNM_FAST5=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230418_Random_Ligation_A_m6A/RL_RG1-6_A_RTA/20230418_1325_X1_AOL616_885f620d/fast5
#MOD_BAM=${ROOTDIR}/${TRAIN_DATASET}_m6A/spomlette_q70.bam
#MOD_FAST5=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230418_Random_Ligation_A_m6A/RL_RG7-12_m6A_RTA/20230418_1325_X2_AOC149_8138c168/fast5
#REF=${ROOTDIR}/reference/${TRAIN_DATASET}_combined.fa
#OUTDIR=${ROOTDIR}/MAFIA_classifiers/${TRAIN_DATASET}_spomlette_${CLASSIFIER}_${SCALER}

### Isabel RL Mix 1 & 3 ###
#UNM_BAM=${ROOTDIR}/RL_Mix1_A_RTA/mapped.bam
#UNM_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix1_A_RTA/*/fast5_*"
#MOD_BAM=${ROOTDIR}/RL_Mix3_m6A_RTA/mapped.bam
#MOD_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix3_m6A_RTA/*/fast5_*"
#REF=${ROOTDIR}/reference/RL_Mix1_Mix3_blocks8.fasta
#OUTDIR=${ROOTDIR}/MAFIA_classifiers/RL_run2_Mix1_Mix3_${SCALER}

### Isabel RL Mix 2 & 4 ###
#UNM_BAM=${ROOTDIR}/RL_Mix2_A_RTA/mapped.bam
#UNM_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix2_A_RTA/*/fast5_*"
#MOD_BAM=${ROOTDIR}/RL_Mix4_m6A_RTA/mapped.bam
#MOD_FAST5="/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/RL_Mix4_m6A_RTA/*/fast5_*"
#REF=${ROOTDIR}/reference/RL_Mix2_Mix4_blocks8.fasta
#OUTDIR=${ROOTDIR}/MAFIA_classifiers/RL_run2_Mix2_Mix4_${SCALER}

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
