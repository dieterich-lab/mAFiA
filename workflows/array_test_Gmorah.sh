#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:turing:1
#SBATCH --mem=120GB
#SBATCH --verbose
#SBATCH --job-name=test_Gmorah_%A_%a
#SBATCH --output=/home/achan/slurm/test_Gmorah_%A_%a.out

### eg: sbatch --array=0-4 --export=challenge=Challenge_1,test_ds=Dataset1 array_test_Gmorah.sh ###

eval "$(/home/achan/miniconda3/condabin/conda shell.bash hook)"
conda activate MAFIA

prj=/prj/TRR319_RMaP/Project_BaseCalling
#challenge=Challenge_1
#test_ds=Dataset1
#batch=0
batch=${SLURM_ARRAY_TASK_ID}

echo "Running Gmorah on ${challenge}, ${test_ds}, ${batch}"

#set -e -u
#
#python3 ${HOME}/git/Gmorah/mAFiA/test_mAFiA.py \
#--ref_file ${prj}/Gm/References/test_reference.fasta \
#--backbone_model_path ${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch \
#--classifier_model_dir ${prj}/Adrian/Gm/Gmorah_models/combined \
#--bam_file ${prj}/Adrian/Gm/Test/${test_ds}/minimap.bam \
#--fast5_dir ${prj}/Adrian/Gm/Test/${test_ds}/fast5_batch${batch} \
#--out_dir ${prj}/Adrian/Gm/Test/${test_ds}/Gmorah_batch${batch} \
#--mod_file ${prj}/Gm/References/Gm_test_sites.bed
