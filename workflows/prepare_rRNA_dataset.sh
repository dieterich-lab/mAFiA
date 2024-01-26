model=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

dataset=HeLa_rRNA
workspace=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/${dataset}
mkdir -p ${workspace}
cd ${workspace}

fast5dir=${workspace}/fast5
mkdir -p ${fast5dir}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

########################################################################################################################
### basecall small number of reads #####################################################################################
########################################################################################################################
srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=4 --mem-per-cpu=64GB \
python3 ${HOME}/git/mAFiA_dev/RODAN/basecall.py \
--fast5dir ${fast5dir} \
--model ${model} \
--batchsize 4096 \
--outdir ${workspace}

########################################################################################################################
#### align to genome ###################################################################################################
########################################################################################################################
ref=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/Gm/reference/rRNA_18S_28S.fasta
sam=${workspace}/mapped.sam
bam=${workspace}/filtered_q50.bam

module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${ref} ${workspace}/rodan.fasta > ${sam}

### check stats and accuracy ###
samtools flagstats ${sam} > qc.txt
${HOME}/git/renata/accuracy.py ${sam} ${ref} >> qc.txt

#### Convert to BAM and index ###
samtools view -bST ${ref} -q50 ${sam} | samtools sort - > ${bam}
samtools index ${bam}