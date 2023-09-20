#!/usr/bin/env bash

test_ds=Dataset1

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Test/${test_ds}
fast5=/prj/TRR319_RMaP/Project_BaseCalling/Gm/Test_Datasets/Challenge_1/${test_ds}
fasta=${workspace}/rodan.fasta

ref=/prj/TRR319_RMaP/Project_BaseCalling/Gm/References/test_reference.fasta
sam=${workspace}/minimap.sam

arch=${HOME}/git/renata/rnaarch
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

#### basecall with Rodan IVT ###
source ${HOME}/git/renata/virtualenv/bin/activate

echo "Basecalling ${fast5}"
srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
python3 -u ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${fast5} \
--arch ${arch} \
--model ${backbone} \
--batchsize 4096 \
--decoder viterbi \
> ${fasta}

### align with minimap ###
deactivate
module purge
module load minimap2

echo "Basecalling finished. Now aligning ${fasta} to ${ref}"
minimap2 --secondary=no -ax map-ont -k14 -t 36 --cs ${ref} ${fasta} > ${sam}

### filter by quality ###
echo "Filtering and converting ${sam}"
sam_filtered=${sam/.sam/_q${filter_score}.sam}
samtools view -h -q${filter_score} ${sam} > ${sam_filtered}

### check read num and accuracy ###
echo "Quality control"
samtools flagstats ${sam_filtered} > ${workspace}/qc_q${filter_score}.txt
${HOME}/git/renata/accuracy.py ${sam_filtered} ${ref} >> ${workspace}/qc_q${filter_score}.txt

cat ${workspace}/qc_q${filter_score}.txt

#### Convert to BAM ###
bam=${sam_filtered//.sam/.bam}
echo "Converting to ${bam}"
samtools view -bST ${ligation_ref} ${sam_filtered} | samtools sort - > ${bam}
samtools index ${bam}

### align with minimap ###
#echo "Aligning with minimap"
#minimap2 -ax map-ont --secondary=no -t 8 ${LIGATION_REF} ${FASTA} | samtools view -F 2324 -b - | samtools sort -o renata.minimap.bam
#samtools index renata.minimap.bam

echo "${dataset} finished"
