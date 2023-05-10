ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

########################################################################################################################
### Wuerzburg batch 1 ##################################################################################################
########################################################################################################################
#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch1_w_splint.fasta

#TRAIN_DATASET=WUE_batch1_A
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_A_RTA/*
#FILTER_SCORE=80

#TRAIN_DATASET=WUE_batch1_m6A
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_m6A_RTA/*
#FILTER_SCORE=70

########################################################################################################################
### Wuerzburg batch 2 ##################################################################################################
########################################################################################################################
REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch2_w_splint.fasta

#TRAIN_DATASET=WUE_batch2_A
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230502_WUE_splint_batch2/WUE_splint_batch2_A_RTA/*
#FILTER_SCORE=80

TRAIN_DATASET=WUE_batch2_m6A
FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230502_WUE_splint_batch2/WUE_splint_batch2_m6A_RTA/*
FILTER_SCORE=70

########################################################################################################################
### Isabel all 6 #######################################################################################################
#REF=${WORKSPACE}/reference/top6_random_permutation_max_blocks_5.fasta

#TRAIN_DATASET=RL_RG1-6_A
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230418_Random_Ligation_A_m6A/RL_RG1-6_A_RTA/20230418_1325_X1_AOL616_885f620d/fast5

#TRAIN_DATASET=RL_RG7-12_m6A
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230418_Random_Ligation_A_m6A/RL_RG7-12_m6A_RTA/20230418_1325_X2_AOC149_8138c168/fast5

### Isabel 3+3 #########################################################################################################
#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/RL_Mix1_Mix3_blocks8.fasta
#TRAIN_DATASET=RL_Mix1_A
#TRAIN_DATASET=RL_Mix3_m6A

#REF=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/RL_Mix2_Mix4_blocks8.fasta
#TRAIN_DATASET=RL_Mix2_A
#TRAIN_DATASET=RL_Mix4_m6A

#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230503_RL_run2/${TRAIN_DATASET}/*

########################################################################################################################
WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TRAIN_DATASET}
mkdir -p ${WORKSPACE}
cd ${WORKSPACE}

FASTA=${WORKSPACE}/basecalled.fasta
SAM=${WORKSPACE}/spomlette.sam

#### basecall with Rodan IVT ###
deactivate
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
BATCH=2
awk '/^>/{p=seen[$0]++}!p' /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_batch${BATCH}_A/ref_recon.fa /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/WUE_batch${BATCH}_m6A/ref_recon.fa > /beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/reference/WUE_batch${BATCH}_combined.fa