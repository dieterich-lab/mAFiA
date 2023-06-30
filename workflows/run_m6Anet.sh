#!/usr/bin/env bash

module load minimap2
module load vbz_compression nanopolish
module load cuda m6anet

#DATASET=WUE_batch3_AB_BA
DATASET=WUE_batch3_ABBA

PRJ=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
NANOPOLISH=${PRJ}/nanopolish/${DATASET}
FAST5=${PRJ}/oligo/${DATASET}/fast5
BASECALL=${PRJ}/oligo/${DATASET}/renata.fasta
LIGATION_REF=${PRJ}/oligo/${DATASET}/ligation_ref.fasta
BAM=${NANOPOLISH}/minimap.bam
m6ANET=${PRJ}/m6Anet/${DATASET}

mkdir -p ${NANOPOLISH}
cd ${NANOPOLISH}

### filter basecall ###
sed '$!N;/>.*\n$/d;P;D' ${BASECALL} > basecall.fasta
nanopolish index -d ${FAST5} basecall.fasta

### copy reference ###
sed -e '1,2d' ${LIGATION_REF} > ref.fasta

### align ###
minimap2 -ax map-ont -t 8 ref.fasta basecall.fasta | samtools sort -o ${BAM}
samtools index ${BAM}

### eventalign ###
nanopolish eventalign \
--reads basecall.fasta \
--bam ${BAM} \
--genome ref.fasta \
--scale-events \
--signal-index \
--summary summary.txt \
--threads 50 \
> eventalign.txt

### m6Anet ###
mkdir -p ${m6ANET}
cd ${m6ANET}

m6anet dataprep --eventalign ${NANOPOLISH}/eventalign.txt --out_dir ${m6ANET} --n_processes 4
m6anet inference --input_dir ${m6ANET} --out_dir ${m6ANET} --n_processes 4 --num_iterations 1000
