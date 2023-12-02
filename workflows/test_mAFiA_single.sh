#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=160GB
#SBATCH --verbose
#SBATCH --job-name=DRACH_weighted_IVT_chr1
#SBATCH --output=/home/achan/slurm/DRACH_weighted_IVT_chr1.out

#ds=Mettl3-KO
#ds=100_WT_0_IVT
ds=0_WT_100_IVT
chr=1

#ds=col0
#ds=vir1
#chr=1

#ds=WT
#ds=IME4_KO
#chr=IV

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/HEK293/${ds}
#workspace=/scratch/achan/Arabidopsis_thaliana/${ds}
#workspace=/scratch/achan/Saccharomyces_cerevisiae/${ds}
#bam=${workspace}/chr${chr}/sorted.chr${chr}.bam
bam=${workspace}/filtered_q50.bam

#fast5_dir=${workspace}/chr${chr}/fast5
fast5_dir=${workspace}/fast5

ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/GRCh38_102/GRCh38_102.chr${chr}.fa
#mod=/home/achan/Data/GLORI/bed_files/GLORI_chr${chr}.tsv
mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/site_annotations/DRACH.GRCh38_102.chr${chr}.bed
#ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/reference/TAIR10_chr_all.fasta
#mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Arabidopsis_thaliana/reference/6motifs.chr1.bed
#ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Saccharomyces_cerevisiae/reference/R64-1-1_96.fa
#mod=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/Saccharomyces_cerevisiae/reference/6motifs.chr${chr}.bed

backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
classifiers=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/MAFIA_classifiers/DRACH_weighted_stripped
#output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH/${ds}/chr${chr}
#output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Arabidopsis_thaliana/${ds}/chr${chr}
#output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/6motifs/Saccharomyces_cerevisiae/${ds}/chr${chr}
output=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/DRACH_weighted/${ds}/chr${chr}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

test_mAFiA \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--ref_file ${ref} \
--mod_file ${mod} \
--min_coverage 50 \
--max_num_reads 2500 \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--mod_prob_thresh 0.5 \
--out_dir ${output}
