#!/usr/bin/env bash
arch=${HOME}/git/renata/rnaarch
backbone=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

challenge=Challenge_1
#test_ds=Dataset1
#test_ds=Dataset2
#test_ds=Dataset3
test_ds=Dataset4

#challenge=Challenge_2
#test_ds=Dataset5
#test_ds=Dataset6

workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/Test/${test_ds}
mkdir -p ${workspace}
cd ${workspace}

fast5=/prj/TRR319_RMaP/Project_BaseCalling/Gm/Test_Datasets/${challenge}/${test_ds}
fasta=${workspace}/rodan.fasta

ref=/prj/TRR319_RMaP/Project_BaseCalling/Gm/References/test_reference.fasta
sam=${workspace}/minimap.sam

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
module purge
module load minimap2

echo "Basecalling finished. Now aligning ${fasta} to ${ref}"
minimap2 --secondary=no -ax map-ont -k14 -t 36 --cs ${ref} ${fasta} > ${sam}

### check read num and accuracy ###
echo "Quality control"
samtools flagstats ${sam} > ${workspace}/qc.txt
${HOME}/git/renata/accuracy.py ${sam} ${ref} >> ${workspace}/qc.txt

cat ${workspace}/qc.txt

#### Convert to BAM ###
bam=${sam//.sam/.bam}
echo "Converting to ${bam}"
samtools view -bST ${ref} ${sam} | samtools view -F 2324 -b - | samtools sort - > ${bam}
samtools index ${bam}

echo "${test_ds} finished"
