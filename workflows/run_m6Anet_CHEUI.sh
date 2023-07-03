#!/usr/bin/env bash

module load minimap2
module load vbz_compression nanopolish

#DATASET=WUE_batch3_AB_BA
DATASET=WUE_batch3_ABBA

PRJ=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
NANOPOLISH=${PRJ}/nanopolish/${DATASET}
FAST5=${PRJ}/oligo/${DATASET}/fast5
FASTQ=${PRJ}/nanopolish/${DATASET}/fastq
LIGATION_REF=${PRJ}/oligo/${DATASET}/ligation_ref.fasta
BAM=${NANOPOLISH}/minimap.bam
m6ANET_OUTDIR=${PRJ}/m6Anet/${DATASET}
CHEUI_OUTDIR=${PRJ}/CHEUI/${DATASET}
GIT_CHEUI=/home/achan/git/CHEUI

######################################################################################
### nanopolish #######################################################################
######################################################################################
mkdir -p ${NANOPOLISH}
cd ${NANOPOLISH}

### index reads ###
cat ${FASTQ}/*.fastq.gz > all.fastq.gz
nanopolish index -s sequencing_summary.txt -d ${FAST5} all.fastq.gz

### copy reference ###
sed -e '1,2d' ${LIGATION_REF} > ref.fasta

### align with minimap ###
minimap2 -ax map-ont --secondary=no -t 8 ref.fasta all.fastq.gz | samtools view -F 2324 -b - | samtools sort -o ${BAM}
samtools index ${BAM}

### eventalign ###
nanopolish eventalign \
--reads all.fastq.gz \
--bam ${BAM} \
--genome ref.fasta \
--scale-events \
--signal-index \
--samples \
--summary summary.txt \
--threads 48 \
> eventalign.txt

######################################################################################
### m6Anet ###########################################################################
######################################################################################
module load cuda m6anet

mkdir -p ${m6ANET_OUTDIR}
cd ${m6ANET_OUTDIR}

rev ${NANOPOLISH}/eventalign.txt | cut -f2- | rev > eventalign.txt

m6anet dataprep --eventalign eventalign.txt --out_dir ${m6ANET_OUTDIR} --n_processes 4

srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
m6anet inference --input_dir ${m6ANET_OUTDIR} --out_dir ${m6ANET_OUTDIR} --n_processes 4 --num_iterations 1000

######################################################################################
### CHEUI ############################################################################
######################################################################################
conda activate CHEUI

mkdir -p ${CHEUI_OUTDIR}
cd ${CHEUI_OUTDIR}

python3 ${GIT_CHEUI}/scripts/CHEUI_preprocess_m6A.py \
-i ${NANOPOLISH}/eventalign.txt \
-m ${GIT_CHEUI}/kmer_models/model_kmer.csv \
-o ${CHEUI_OUTDIR}/out_A_signals+IDs.p \
-n 15

srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
python3 ${GIT_CHEUI}/scripts/CHEUI_predict_model1.py \
-i ${CHEUI_OUTDIR}/out_A_signals+IDs.p/eventalign_signals+IDS.p \
-m ${GIT_CHEUI}/CHEUI_trained_models/CHEUI_m6A_model1.h5 \
-o ${CHEUI_OUTDIR}/read_level_m6A_predictions.txt \
-l ${DATASET}