#!/usr/bin/env bash
shopt -s globstar

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

PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/m6A
ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
FILTER_SCORE=70
#REF=${PRJ_DIR}/reference/${ORIG}_oligo_ref_${RUN}.fasta
REF=${PRJ_DIR}/reference/ISA_m6A_oligos_allDRACH.fasta
DATASET=${ORIG}_${RUN}_${MOD}
WORKSPACE=${PRJ_DIR}/oligo/${DATASET}
READS=${WORKSPACE}/fast5
FASTA=${WORKSPACE}/renata.fasta
SAM=${WORKSPACE}/spomelette.sam
LIGATION_REF=${WORKSPACE}/ligation_ref.fasta

### softlink fast5 files ###
#mkdir -p ${READS} && cd "$_"
#echo "Creating softlinks from ${LOC}"
#for f in ${LOC}/**/*.fast5
#do
#  ln -s $f
#done
#cd ${WORKSPACE}

### check links ###
#for my_link in ${READS}/*.fast5
#do
#  if [ -L ${my_link} ]
#  then
#     if [ -e ${my_link} ]
#     then
#        echo "${my_link} Good link"
#     else
#        echo "${my_link} Broken link"
#     fi
#  elif [ -e ${my_link} ]
#  then
#     echo "${my_link} Not a link"
#  else
#     echo "${my_link} Missing"
#  fi
#done

### convert pod5 back to fast5 ###
#for f in ${LOC}/**/*.pod5; do pod5 convert to_fast5 $f --output ${READS}; done

#### basecall with Rodan IVT ###
#source ${HOME}/git/mAFiA/mafia-venv/bin/activate
#
#echo "Basecalling ${READS}"
#srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
#python3 -u ${HOME}/git/mAFiA/RODAN/basecall.py \
#--fast5dir ${READS} \
#--model ${MODEL} \
#--batchsize 4096 \
#--outdir ${WORKSPACE}

### align with spomelette ###
echo "Basecalling finished. Now aligning ${FASTA} to ${REF}"
python3 -u ${HOME}/git/mAFiA_dev/misc/spanish_omelette_alignment.py \
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


### Convert to BAM ###
BAM=${FILTERED_SAM//.sam/.bam}
echo "Converting to ${BAM}"
samtools view -bST ${LIGATION_REF} ${FILTERED_SAM} | samtools sort - > ${BAM}
samtools index ${BAM}

### align with minimap ###
#echo "Aligning with minimap"
#minimap2 -ax map-ont --secondary=no -t 8 ${LIGATION_REF} ${FASTA} | samtools view -F 2324 -b - | samtools sort -o renata.minimap.bam
#samtools index renata.minimap.bam

echo "${DATASET} finished"

#### merge ligation ref ###
#for RUN in ISA_mix1 ISA_mix2 ISA_mix3 ISA_mix4
#do
#  MERGE_REF_INPUTS=""
#  for MOD in A m6A
#  do
#    MERGE_REF_INPUTS+="${PRJ_DIR}/oligo/${RUN}_${MOD}/ligation_ref.fasta "
#  done
#  mkdir ${PRJ_DIR}/oligo/ligation_ref
#  MERGE_REF_OUTPUT="${PRJ_DIR}/oligo/ligation_ref/ligation_ref_${RUN}_A_m6A.fasta"
#  awk '/^>/{p=seen[$0]++}!p' ${MERGE_REF_INPUTS} > ${MERGE_REF_OUTPUT}
#done
