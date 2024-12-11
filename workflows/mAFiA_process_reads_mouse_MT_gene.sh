#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=psico-mAFiA_mouse_MT
#SBATCH --output=/home/achan/slurm/psico-mAFiA_mouse_MT_%A.out

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
chr="X"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 24 ]]
then
chr="MT"
else
chr=${SLURM_ARRAY_TASK_ID}
fi


#workspace=/prj/Dewenter_TAC_Backs_lab/achan/${ds}
#workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/${ds}
#bam=${workspace}/chr${chr}/sorted.chr${chr}.bam
#fast5_dir=${workspace}/chr${chr}/fast5

chr=MT
ds=HFpEF
cond=HFpEF
bam=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/${ds}/mt-Ti/${cond}_merged_mt-Ti.bam
fast5_dir=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/${ds}/mt-Ti/fast5_${cond}
output=/beegfs/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/${ds}/${cond}_merged_MT-Ti


mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.chr${chr}.bed
backbone=${HOME}/git/mAFiA_dev/models/RODAN_HEK293_IVT.torch
classifiers=${HOME}/git/mAFiA_dev/models/psi-co-mAFiA

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_process_reads_parallel.py \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--mod_file ${mod} \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--num_jobs 4 \
--batchsize 256 \
--out_dir ${output}
