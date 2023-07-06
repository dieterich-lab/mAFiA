#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-g4-1
#SBATCH --gres=gpu:turing:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80GB
#SBATCH --verbose
#SBATCH --job-name=m6Anet_CHEUI
#SBATCH --output=/home/achan/slurm/m6Anet_CHEUI_%A.out

DATASET=100_WT_0_IVT
#DATASET=75_WT_25_IVT
#DATASET=50_WT_50_IVT
#DATASET=25_WT_75_IVT
#DATASET=0_WT_100_IVT

echo ${DATASET}

PRJ=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
NANOPOLISH=${PRJ}/nanopolish/${DATASET}
FAST5=${PRJ}/HEK293/${DATASET}/fast5
FASTQ=${PRJ}/nanopolish/${DATASET}/fastq
REF=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.cdna.all.fa
BAM=${NANOPOLISH}/minimap.bam
m6ANET_OUTDIR=${PRJ}/m6Anet/${DATASET}
CHEUI_OUTDIR=${PRJ}/CHEUI/${DATASET}
GIT_CHEUI=/home/achan/git/CHEUI

######################################################################################
### copy fastq #######################################################################
######################################################################################
#mkdir -p ${FASTQ} && cd ${FASTQ}
#
#for f in ${FAST5}/*.fast5
#do
#  f5_path=$(realpath $f)
#  fastq_path=${f5_path//fast5/fastq}.gz
#  ln -s ${fastq_path}
#
#  parent_dir=$(dirname $(dirname $f5_path))
#  seq_summary_name=$(basename $(echo ${parent_dir}/sequencing_summary*.txt))
#  if [ ! -f "${NANOPOLISH}/${seq_summary_name}" ]; then
#    cp ${parent_dir}/sequencing_summary*.txt ${NANOPOLISH}
#  fi
#done
#
#cd ${NANOPOLISH}
#cat ${FASTQ}/*.fastq.gz > all.fastq.gz
#awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' sequencing_summary_*.txt > sequencing_summary.txt
#rm sequencing_summary_*.txt

######################################################################################
### nanopolish #######################################################################
######################################################################################
#module load minimap2
#module load vbz_compression nanopolish
#
#### index reads ###
#echo "Indexing reads"
#nanopolish index -s sequencing_summary.txt -d ${FAST5} all.fastq.gz
#
#### align with minimap ###
#minimap2 -ax map-ont --secondary=no -t 8 ${REF} all.fastq.gz | samtools view -F 2324 -b - | samtools sort -o ${BAM}
#samtools index ${BAM}
#
#### eventalign ###
#echo "Eventalign..."
#nanopolish eventalign \
#--reads all.fastq.gz \
#--bam ${BAM} \
#--genome ${REF} \
#--scale-events \
#--signal-index \
#--samples \
#--summary summary.txt \
#--threads 48 \
#> eventalign.txt

######################################################################################
### m6Anet ###########################################################################
######################################################################################
#module load cuda m6anet
#
#mkdir -p ${m6ANET_OUTDIR}
#cd ${m6ANET_OUTDIR}
#
#rev ${NANOPOLISH}/eventalign.txt | cut -f2- | rev > eventalign.txt
#echo "m6Anet..."
#m6anet dataprep --eventalign eventalign.txt --out_dir ${m6ANET_OUTDIR} --n_processes 4
#m6anet inference --input_dir ${m6ANET_OUTDIR} --out_dir ${m6ANET_OUTDIR} --n_processes 4 --num_iterations 1000

######################################################################################
### CHEUI ############################################################################
######################################################################################
eval "$(conda shell.bash hook)"
conda activate CHEUI

mkdir -p ${CHEUI_OUTDIR}
cd ${CHEUI_OUTDIR}

echo "CHEUI preprocessing..."
python3 ${GIT_CHEUI}/scripts/CHEUI_preprocess_m6A.py \
-i ${NANOPOLISH}/eventalign.txt \
-m ${GIT_CHEUI}/kmer_models/model_kmer.csv \
-o ${CHEUI_OUTDIR}/out_A_signals+IDs.p \
-n 15

echo "CHEUI predict..."
python3 ${GIT_CHEUI}/scripts/CHEUI_predict_model1.py \
-i ${CHEUI_OUTDIR}/out_A_signals+IDs.p/eventalign_signals+IDS.p \
-m ${GIT_CHEUI}/CHEUI_trained_models/CHEUI_m6A_model1.h5 \
-o ${CHEUI_OUTDIR}/read_level_m6A_predictions.txt \
-l ${DATASET}