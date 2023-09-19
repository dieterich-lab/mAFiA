#!/usr/bin/env bash

#dataset=A1-4
#dataset=A5-8
#dataset=A9-12

dataset=B1-4
#dataset=B5-8

#dataset=C1-4
#dataset=C5-8
#dataset=C9-12

#dataset=D1-4
#dataset=D5-8
#dataset=D9-12

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm
fast5dir=/prj/TRR319_RMaP/Project_BaseCalling/Gm/Train_Datasets/${dataset}_Ligation/Fast5_Files
fasta=${workspace}/${dataset}/rodan.fasta

oligo_ref=${workspace}/oligo_reference/oligo_${dataset}.fasta
sam=${workspace}/${dataset}/spomelette.sam
ligation_ref=${workspace}/${dataset}/ligation_ref.fasta
filter_score=70

arch=${HOME}/git/renata/rnaarch
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

### softlink to fast5 pass files ###
#fast5pass=${workspace}/${dataset}/fast5_pass
#mkdir -p ${fast5pass}
#cd ${fast5pass}
#for f in ${fast5dir}/*pass*.fast5; do ln -s $f; done

#### basecall with Rodan IVT ###
#source ${HOME}/git/renata/virtualenv/bin/activate
#
#echo "Basecalling ${fast5pass}"
#srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
#python3 -u ${HOME}/git/renata/basecall_viterbi.py \
#--fast5dir ${fast5pass} \
#--arch ${arch} \
#--model ${backbone} \
#--batchsize 4096 \
#--decoder viterbi \
#> ${fasta}

### align with spomelette ###
deactivate
eval "$(/home/achan/miniconda3/condabin/conda shell.bash hook)"
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
echo "Filtering and converting ${sam}"
sam_filtered=${sam/.sam/_q${filter_score}.sam}
samtools view -h -q${filter_score} ${sam} > ${sam_filtered}

### check read num and accuracy ###
echo "Quality control"
samtools flagstats ${sam_filtered} > ${workspace}/${dataset}/qc_q${filter_score}.txt
${HOME}/git/renata/accuracy.py ${sam_filtered} ${ligation_ref} >> ${workspace}/${dataset}/qc_q${filter_score}.txt


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
