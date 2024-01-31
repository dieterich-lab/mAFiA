#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1,gpu-g2-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=180GB
#SBATCH --verbose
#SBATCH --job-name=PSICO_HEK293
#SBATCH --output=/home/achan/slurm/PSICO_HEK293_%A.out

#ds=Mettl3-KO
ds="100_WT_0_IVT"
#ds=75_WT_25_IVT
#ds=50_WT_50_IVT
#ds=25_WT_75_IVT
#ds=0_WT_100_IVT
#chr=X

#ds=JK_HEK293_DMSO_merged
#ds=JK_HEK293_STM2457_merged

workspace="/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/${ds}"
bam="${workspace}/genome_filtered_q50.bam"
fast5_dir="${workspace}/fast5"

ref="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"
mod="/prj/TRR319_RMaP/Project_BaseCalling/Adrian/psU/site_annotations/BID_seq.bed"

backbone="${HOME}/git/mAFiA/models/RODAN_HEK293_IVT.torch"
classifiers="/prj/TRR319_RMaP/Project_BaseCalling/Adrian/psU/PSICO_classifiers/chosen8"

output="/prj/TRR319_RMaP/Project_BaseCalling/Adrian/psU/PSICO_inference/HEK293/${ds}/chosen8"

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

python3 ${HOME}/git/mAFiA_dev/mAFiA/test_PSICO.py \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--ref_file ${ref} \
--mod_file ${mod} \
--min_coverage 20 \
--max_num_reads 1000 \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--out_dir ${output}
