ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
FILTER_SCORE=80

########################################################################################################################
### Wuerzburg batch 1 ##################################################################################################
########################################################################################################################
#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_oligos_batch1.fasta

#TRAIN_DATASET=WUE_batch1_A
##FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_A_RTA/*

#TRAIN_DATASET=WUE_batch1_m6A
##FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_m6A_RTA/*

#FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}/fast5

########################################################################################################################
### Wuerzburg batch 2 ##################################################################################################
########################################################################################################################
#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_oligos_batch2.fasta

#TRAIN_DATASET=WUE_batch2_A
##FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230502_WUE_splint_batch2/WUE_splint_batch2_A_RTA/*

#TRAIN_DATASET=WUE_batch2_m6A
##FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230502_WUE_splint_batch2/WUE_splint_batch2_m6A_RTA/*
##FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/${RUN}/WUE_splint_m6A_batch2_RTA_*/*

#FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}/fast5

########################################################################################################################
### Isabel all 6 #######################################################################################################
########################################################################################################################
#REF=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/ISA_run1.fasta
#
#TRAIN_DATASET=ISA_run1_A
#TRAIN_DATASET=ISA_run1_m6A
#
#FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}/fast5
#
#HOMOPOLYMER=1

########################################################################################################################
### Isabel 3+3 #########################################################################################################
########################################################################################################################
#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/RL_Mix1_Mix3_blocks8.fasta
#TRAIN_DATASET=RL_Mix1_A
#TRAIN_DATASET=RL_Mix3_m6A

#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/RL_Mix2_Mix4_blocks8.fasta
#TRAIN_DATASET=RL_Mix2_A
#TRAIN_DATASET=RL_Mix4_m6A

#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/${TRAIN_DATASET}/*

########################################################################################################################
### Isabel hetero ######################################################################################################
########################################################################################################################
REF=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/ISA_M4_M5.fasta
HOMOPOLYMER=0

#TRAIN_DATASET=RL_M4_M5
#TRAIN_DATASET=RL_M4_M5star
#TRAIN_DATASET=RL_M4star_M5
TRAIN_DATASET=RL_M4star_M5star

FAST5_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}/fast5

########################################################################################################################

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}
mkdir -p ${WORKSPACE}
cd ${WORKSPACE}

FASTA=${WORKSPACE}/renata.fasta
SAM=${WORKSPACE}/spomelette.sam

#### basecall with Rodan IVT ###
#deactivate
source ${HOME}/git/renata/virtualenv/bin/activate

srun --partition=gpu --gres=gpu:turing:1 --cpus-per-task=8 --mem-per-cpu=8GB \
python3 ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${FAST5_DIR} \
--arch ${ARCH} \
--model ${MODEL} \
--batchsize 2048 \
--decoder viterbi \
> ${FASTA} &

#### align with minimap ###
#module purge
#module load minimap2
#minimap2 --secondary=no -ax map-ont -t 36 --cs ${REF} ${FASTA} > ${SAM}

### align with spomlette ###
python3 ${HOME}/git/MAFIA/spanish_omelette_alignment.py \
--ref_file ${REF} \
--query_file ${FASTA} \
--recon_ref_file ${WORKSPACE}/ref_recon.fa \
--sam_file ${SAM} \
--homopolymer ${HOMOPOLYMER} \
--write_cs

### filter by quality ###
FILTERED_SAM=${SAM/.sam/_q${FILTER_SCORE}.sam}
samtools view -h -q${FILTER_SCORE} ${SAM} > ${FILTERED_SAM}

### check read num and accuracy ###
samtools flagstats ${FILTERED_SAM}
${HOME}/git/renata/accuracy.py ${FILTERED_SAM} ${WORKSPACE}/ref_recon.fa


#### Convert to BAM ###
BAM=${FILTERED_SAM//.sam/.bam}
samtools view -bST ${REF} ${FILTERED_SAM} | samtools sort - > ${BAM}
samtools index ${BAM}

### filter reads by max indel ###
#deactivate
#conda activate MAFIA
#
#python3 ${HOME}/git/MAFIA/filter_bam_file_by_max_indel_len.py \
#--infile ${BAM} \
#--outfile ${BAM}.filtered \
#--indel_thresh 10
#
#samtools index ${BAM}.filtered

### combine recon references ###
#BATCH=1
#BATCH=2
#awk '/^>/{p=seen[$0]++}!p' /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_batch${BATCH}_A/ref_recon.fa /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_batch${BATCH}_m6A/ref_recon.fa > /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch${BATCH}_ref_recon.fa
#awk '/^>/{p=seen[$0]++}!p' /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/ISA_run1_A/ref_recon.fa /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/ISA_run1_m6A/ref_recon.fa > /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/ISA_run1_ref_recon.fa