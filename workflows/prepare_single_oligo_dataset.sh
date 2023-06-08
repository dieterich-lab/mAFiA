#!/usr/bin/env bash

while getopts o:r:m:h:l: flag
do
  case "${flag}" in
    o) ORIG=${OPTARG};;
    r) RUN=${OPTARG};;
    m) MOD=${OPTARG};;
    h) HOMOPOLYMER=${OPTARG};;
    l) LOC=${OPTARG};;
  esac
done

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/oligo
ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
FILTER_SCORE=70
REF=${PRJ_DIR}/reference/${ORIG}_oligo_ref_${RUN}.fasta
DATASET=${ORIG}_${RUN}_${MOD}
WORKSPACE=${PRJ_DIR}/${DATASET}
FAST5_DIR=${WORKSPACE}/fast5
FASTA=${WORKSPACE}/renata.fasta
SAM=${WORKSPACE}/spomelette.sam
LIGATION_REF=${WORKSPACE}/ligation_ref.fasta

### softlink fast5 files ###
mkdir -p ${FAST5_DIR} && cd "$_"
echo "Creating softlinks from ${LOC}"
for f in ${LOC}/*/fast5_*/*.fast5
do
  ln -s $f
done
cd ${WORKSPACE}

#### basecall with Rodan IVT ###
source ${HOME}/git/renata/virtualenv/bin/activate

echo "Basecalling ${FAST5_DIR}"
srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
python3 -u ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${FAST5_DIR} \
--arch ${ARCH} \
--model ${MODEL} \
--batchsize 2048 \
--decoder viterbi \
> ${FASTA}

### align with spomelette ###
echo "Basecalling finished. Now aligning ${FASTA} to ${REF}"
python3 -u ${HOME}/git/MAFIA/spanish_omelette_alignment.py \
--ref_file ${REF} \
--query_file ${FASTA} \
--recon_ref_file ${LIGATION_REF} \
--sam_file ${SAM} \
--homopolymer ${HOMOPOLYMER} \
--write_cs

### filter by quality ###
echo "Filtering and converting ${SAM}"
FILTERED_SAM=${SAM/.sam/_q${FILTER_SCORE}.sam}
samtools view -h -q${FILTER_SCORE} ${SAM} > ${FILTERED_SAM}

### check read num and accuracy ###
echo "Quality control"
samtools flagstats ${FILTERED_SAM} > ${WORKSPACE}/qc_q${FILTER_SCORE}.txt
${HOME}/git/renata/accuracy.py ${FILTERED_SAM} ${LIGATION_REF} >> ${WORKSPACE}/qc_q${FILTER_SCORE}.txt


#### Convert to BAM ###
BAM=${FILTERED_SAM//.sam/.bam}
echo "Converting to ${BAM}"
samtools view -bST ${LIGATION_REF} ${FILTERED_SAM} | samtools sort - > ${BAM}
samtools index ${BAM}

echo "${DATASET} finished"