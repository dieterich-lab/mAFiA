#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:turing:1
#SBATCH --mem=250GB
#SBATCH --verbose
#SBATCH --job-name=test_Gmorah_Dataset1
#SBATCH --output=/home/achan/slurm/test_Gmorah_Dataset1.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

prj=/prj/TRR319_RMaP/Project_BaseCalling
challenge=Challenge_1
test_ds=Dataset1

set -e -u

python3 ${HOME}/git/Gmorah/mAFiA/test_mAFiA.py \
--ref_file ${prj}/Gm/References/test_reference.fasta \
--backbone_model_path ${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch \
--classifier_model_dir ${prj}/Adrian/Gm/Gmorah_models/combined \
--bam_file ${prj}/Adrian/Gm/Test/${test_ds}/minimap.bam \
--fast5_dir ${prj}/Gm/Test_Datasets/${challenge}/${test_ds} \
--out_dir ${prj}/Adrian/Gm/Test/${test_ds}/Gmorah \
--mod_file ${prj}/Gm/References/Gm_test_sites.bed
