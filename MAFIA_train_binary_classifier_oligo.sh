#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu-g3-1
#SBATCH --mem=192GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_train_binary_classifier_oligo.out

eval "$(conda shell.bash hook)"
#eval "$(/home/achan/miniconda3/condabin/conda shell.bash hook)"
conda activate MAFIA

set -e -u

cd ${HOME}/git/MAFIA
WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian

python3 train_binary_classifier.py \
--unm_bam_file ${WORKSPACE}/mapping/A_RTA_sorted_filtered.bam \
--unm_fast5_dir /prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_A_RTA/20230221_1328_X1_ANS648_701f60ca \
--mod_bam_file ${WORKSPACE}/mapping/m6A_RTA_sorted_filtered.bam \
--mod_fast5_dir /prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_m6A_RTA/20230221_1328_X2_ANS491_f891b4b9 \
--ref_file ${WORKSPACE}/reference/splint_variations_max_blocks_7.fasta \
--max_num_reads -1 \
--min_coverage 0 \
--backbone_model_path ${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch \
--extraction_layer convlayers.conv21 \
--feature_width 0 \
--classifier logistic_regression \
--classifier_model_dir ${WORKSPACE}/MAFIA_classifiers/A_m6A_NoNorm_PassFail
