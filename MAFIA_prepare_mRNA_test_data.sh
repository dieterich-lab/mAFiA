sterm -c 2 -m 64GB -N gpu-g3-1

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian
REF=${HOME}/Data/genomes/GRCh38_96.fa

ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

#####################################################################################################################################
#DATASET=HEK293A_WT
#FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all

DATASET=HEK293_IVT
FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293_IVT_2/20211201_1116_X5_FAR06706_305e3998/fast5_all
#####################################################################################################################################

mkdir -p ${WORKSPACE}
FASTA=${WORKSPACE}/basecall/${DATASET}.fasta
SAM=${WORKSPACE}/mapping/${DATASET}.sam
BAM=${SAM//.sam/.bam}


### basecall large number of reads ###
source ${HOME}/git/renata/virtualenv/bin/activate

FILENAME_PREFIX=${DATASET}_fast5_paths_part
ls -1 ${FAST5_DIR}/*.fast5 > ${WORKSPACE}/${DATASET}_fast5_paths_all
split -l5 -d ${WORKSPACE}/${DATASET}_fast5_paths_all ${WORKSPACE}/${FILENAME_PREFIX}
NUM_STRAND_FILES=`ls ${WORKSPACE}/${DATASET}_fast5_paths_part* | wc -l`
sbatch --array=0-$((${NUM_STRAND_FILES}-1)) --export=ALL,WORKSPACE=${WORKSPACE},FAST5_DIR=${FAST5_DIR},FILENAME_PREFIX=${FILENAME_PREFIX},FASTA=${FASTA},ARCH=${ARCH},MODEL=${MODEL} array_basecaller.sh

cat ${FASTA}* > ${FASTA}_merged

#### align and check accuracy ###
module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF} ${FASTA}_merged > ${SAM}

${HOME}/git/renata/accuracy.py ${SAM} ${REF}

samtools flagstats ${SAM}

#### Convert to BAM ###
samtools view -bST ${REF} ${SAM} -o ${BAM}
samtools sort ${BAM} > ${BAM//.bam/_sorted.bam}
samtools index ${BAM//.bam/_sorted.bam}
