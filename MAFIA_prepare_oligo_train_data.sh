sterm -c 2 -m 64GB -N gpu-g3-1

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian
REF=${WORKSPACE}/reference/splint_variations_max_blocks_7.fasta
ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

#####################################################################################################################################
#DATASET=A_RTA
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_A_RTA/20230221_1328_X1_ANS648_701f60ca

DATASET=m6A_RTA
FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230221_WUE_splint_lig/WUE_splint_lig_m6A_RTA/20230221_1328_X2_ANS491_f891b4b9

#####################################################################################################################################

mkdir -p ${WORKSPACE}
#FAST5_DIR=${WORKSPACE}/fast5_all/${DATASET}
#mkdir -p ${FAST5_DIR}
#cp ${RAW_FAST5_DIR}/fast5_pass/*.fast5 ${FAST5_DIR}
#cp ${RAW_FAST5_DIR}/fast5_fail/*.fast5 ${FAST5_DIR}

FASTA=${WORKSPACE}/basecall/${DATASET}.fasta
SAM=${WORKSPACE}/mapping/${DATASET}.sam
BAM=${SAM//.sam/.bam}

#### basecall with Rodan IVT ###
source ${HOME}/git/renata/virtualenv/bin/activate

python3 ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${FAST5_DIR} \
--arch ${ARCH} \
--model ${MODEL} \
> ${FASTA}

#### align and check accuracy ###
module purge
module load minimap2
minimap2 --secondary=no -ax map-ont -t 36 --cs ${REF} ${FASTA} > ${SAM}

samtools flagstats ${SAM}

${HOME}/git/renata/accuracy.py ${SAM} ${REF}


#### Convert to BAM ###
samtools view -bST ${REF} ${SAM} -o ${BAM}
samtools sort ${BAM} > ${BAM//.bam/_sorted.bam}
samtools index ${BAM//.bam/_sorted.bam}

### filter reads by max indel ###
deactivate
conda activate MAFIA

python3 ${HOME}/git/MAFIA/filter_bam_file_by_max_indel_len.py \
--infile ${BAM//.bam/_sorted.bam} \
--outfile ${BAM//.bam/_sorted_filtered.bam} \
--indel_thresh 10

samtools index ${BAM//.bam/_sorted_filtered.bam}
