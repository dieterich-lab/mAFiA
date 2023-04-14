#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --mem=60GB
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --output=/home/achan/slurm/MAFIA_cluster_BID_seq_%A.out

eval "$(conda shell.bash hook)"
conda activate MAFIA

set -e -u

cd ${HOME}/git/MAFIA

python3 ${HOME}/git/MAFIA/mRNA_train_clustering.py \
--unm_bam_file /prj/TRR319_RMaP/Project_BaseCalling/Adrian/mapping/HEK293_IVT.bam.sorted \
--unm_fast5_dir /prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293_IVT_2/fast5_pass \
--mod_bam_file /prj/TRR319_RMaP/Project_BaseCalling/Adrian/mapping/HEK293A_WT.bam.sorted \
--mod_fast5_dir /prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all \
--ref_file ${HOME}/Data/genomes/GRCh38_96.fa \
--mod_file ${HOME}/Data/BID_seq/GSE179798_HEK293T_mRNA_WT_BID-seq.xlsx \
--max_num_reads -1 \
--min_coverage 50 \
--backbone_model_path ${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch \
--extraction_layer convlayers.conv21 \
--feature_width 0 \
--clustering_model_dir ${HOME}/img_out/MAFIA_clustering/BID_seq