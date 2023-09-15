#!/usr/bin/env bash

dataset=A1-4
workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm
fast5dir=/prj/TRR319_RMaP/Project_BaseCalling/Gm/Train_Datasets/${dataset}_Ligation/Fast5_Files
fasta=${workspace}/${dataset}/rodan.fasta

oligo_ref=${workspace}/oligo_reference/oligo_${dataset}.fasta
sam=${workspace}/${dataset}/spomelette.sam
ligation_ref=${workspace}/${dataset}/ligation_ref.fasta

arch=${HOME}/git/renata/rnaarch
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

### softlink to fast5 pass files ###
fast5pass=${workspace}/${dataset}/fast5_pass
mkdir -p ${fast5pass}
cd ${fast5pass}
for f in ${fast5dir}/*pass*.fast5; do ln -s $f; done

#### basecall with Rodan IVT ###
source ${HOME}/git/renata/virtualenv/bin/activate

echo "Basecalling ${fast5pass}"
srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
python3 -u ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${fast5pass} \
--arch ${arch} \
--model ${backbone} \
--batchsize 2048 \
--decoder viterbi \
> ${fasta}

### align with spomelette ###
deactivate
conda activate MAFIA

echo "Basecalling finished. Now aligning ${fasta} to ${oligo_ref}"
python3 -u ${HOME}/git/Gmorah/oligo/spanish_omelette_alignment.py \
--ref_file ${oligo_ref} \
--query_file ${fasta} \
--recon_ref_file ${ligation_ref} \
--sam_file ${sam} \
--homopolymer 0 \
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

### align with minimap ###
echo "Aligning with minimap"
minimap2 -ax map-ont --secondary=no -t 8 ${LIGATION_REF} ${FASTA} | samtools view -F 2324 -b - | samtools sort -o renata.minimap.bam
samtools index renata.minimap.bam

echo "${DATASET} finished"
