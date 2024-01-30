#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=mESC_WT_STM_merged_chrX
#SBATCH --output=/home/achan/slurm/mESC_WT_STM_merged_chrX.out

#ds=Mettl3-KO
#ds=100_WT_0_IVT
#ds=75_WT_25_IVT
#ds=50_WT_50_IVT
#ds=25_WT_75_IVT
#ds=0_WT_100_IVT
#chr=X

#ds=JK_HEK293_DMSO_merged
#ds=JK_HEK293_STM2457_merged

#ds=mESC_Mettl3_KO_merged
#ds=mESC_WT_DMSO_merged
ds=mESC_WT_STM_merged

chr=X
#chr=${SLURM_ARRAY_TASK_ID}

#ds=col0
#ds=vir1
#chr=1

#ds=WT
#ds=IME4_KO
#chr=IV

workspace=/prj/TRR319_RMaP/Project_B01/Adrian/${ds}
#workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/${ds}
#workspace=/scratch/achan/Arabidopsis_thaliana/${ds}
#workspace=/scratch/achan/Saccharomyces_cerevisiae/${ds}
#bam=${workspace}/genome_filtered_q50.bam
#fast5_dir=${workspace}/fast5

fast5_dir=${workspace}/chr${chr}/fast5
bam=${workspace}/chr${chr}/sorted.chr${chr}.bam

#ref=/beegfs/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa
#ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/GRCh38_102/GRCh38_102.chr${chr}.fa
#mod=/beegfs/homes/achan/Data/site_annotations/homo_sapiens/GRCh38_102/DRACH.GRCh38_102.chr${chr}.bed

ref=/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa
mod=/beegfs/homes/achan/Data/site_annotations/mus_musculus/GRCm38_102/DRACH.GRCm38_102.chr${chr}.bed

#mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/site_annotations/GLORI_all_ref5mer.bed

#mod=/home/achan/Data/GLORI/bed_files/GLORI_chr${chr}.tsv
#mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/site_annotations/DRACH.GRCh38_102.chr${chr}.bed
#ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/reference/TAIR10_chr_all.fasta
#mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/reference/6motifs.chr1.bed
#ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Saccharomyces_cerevisiae/reference/R64-1-1_96.fa
#mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Saccharomyces_cerevisiae/reference/6motifs.chr${chr}.bed

backbone=${HOME}/git/mAFiA/models/RODAN_HEK293_IVT.torch
classifiers=${HOME}/git/mAFiA/models/DRACH_v1

#backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
#classifiers=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/DRACH_v1
#output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH/${ds}/chr${chr}
#output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Arabidopsis_thaliana/${ds}/chr${chr}
#output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Saccharomyces_cerevisiae/${ds}/chr${chr}
#output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_v1/${ds}

#output=${workspace}/mAFiA/chr${SLURM_ARRAY_TASK_ID}
output=${workspace}/mAFiA/chr${chr}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

test_mAFiA \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--ref_file ${ref} \
--mod_file ${mod} \
--min_coverage 20 \
--max_num_reads 1000 \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--out_dir ${output}
