#!/usr/bin/env bash

sterm -N gpu-g3-1 -c 4 -m 64G

module load minimap2
module load vbz_compression nanopolish
module load cuda m6anet

DATASET=WUE_batch3_AB_BA
#DATASET=WUE_batch3_ABBA

PRJ=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
NANOPOLISH=${PRJ}/nanopolish/${DATASET}
FAST5=${PRJ}/oligo/${DATASET}/fast5
BASECALL=${PRJ}/oligo/${DATASET}/renata.fasta
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

### filter basecall ###
sed '$!N;/>.*\n$/d;P;D' ${BASECALL} > basecall.fasta
nanopolish index -s sequencing_summary.txt -d ${FAST5} basecall.fasta

### copy reference ###
sed -e '1,2d' ${LIGATION_REF} > ref.fasta

### align ###
minimap2 -ax map-ont --secondary=no -t 8 ref.fasta basecall.fasta | samtools sort -o ${BAM}
samtools index ${BAM}

### eventalign ###
nanopolish eventalign \
--reads basecall.fasta \
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
mkdir -p ${m6ANET_OUTDIR}
cd ${m6ANET_OUTDIR}

rev ${NANOPOLISH}/eventalign.txt | cut -f2- | rev > eventalign.txt

m6anet dataprep --eventalign eventalign.txt --out_dir ${m6ANET_OUTDIR} --n_processes 4
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

python3 ${GIT_CHEUI}/scripts/CHEUI_predict_model1.py \
-i ${CHEUI_OUTDIR}/out_A_signals+IDs.p/eventalign_signals+IDS.p \
-m ${GIT_CHEUI}/CHEUI_trained_models/CHEUI_m6A_model1.h5 \
-o ${CHEUI_OUTDIR}/read_level_m6A_predictions.txt \
-l ${DATASET}