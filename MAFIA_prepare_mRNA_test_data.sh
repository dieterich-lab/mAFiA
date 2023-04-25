WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian
REF=${HOME}/Data/genomes/GRCh38_96.fa
ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

### old HEK293 ######################################################################################################################
#DATASET=HEK293A_WT
#FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all

#DATASET=HEK293_IVT
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293_IVT_2/fast5_pass

### new HEK293 ######################################################################################################################
DATASET=100_WT_0_IVT_RTA
#DATASET=0_WT_100_IVT_RTA
FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230419_HEK293_WT_IVT_Mix/${DATASET}/*/fast5_*

#####################################################################################################################################
#DATASET=HEK293T-WT-0-rep2
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-0-rep2/fast5_pass

#DATASET=HEK293T-WT-25-rep1
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-25-rep1/fast5_pass

#DATASET=HEK293T-WT-50-rep2
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-50-rep2/fast5_pass

#DATASET=HEK293T-WT-50-rep3
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-50-rep3/fast5

#DATASET=HEK293T-WT-75-rep4
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-75-rep4/fast5

#DATASET=HEK293T-WT-100-rep1
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-100-rep1/fast5_pass

#####################################################################################################################################
#DATASET=HEK293T-Mettl3-KO-rep3
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-Mettl3-KO-rep3/fast5

#DATASET=HEK293T-WT-rep3
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293T-WT-Mettl3-Mix/HEK293T-WT-rep3/fast5

#####################################################################################################################################

mkdir -p ${WORKSPACE}
FASTA=${WORKSPACE}/basecall/${DATASET}.fasta
SAM=${WORKSPACE}/mapping/${DATASET}.sam
BAM=${SAM//.sam/.bam}

source ${HOME}/git/renata/virtualenv/bin/activate

### basecall large number of reads ###
FILENAME_PREFIX=${DATASET}_fast5_paths_part
ls -1 ${FAST5_DIR}/*.fast5 > ${WORKSPACE}/${DATASET}_fast5_paths_all
split -l5 -d ${WORKSPACE}/${DATASET}_fast5_paths_all ${WORKSPACE}/${FILENAME_PREFIX}

#NUM_STRAND_FILES=`ls ${WORKSPACE}/${DATASET}_fast5_paths_part* | wc -l`
#sbatch --array=0-$((${NUM_STRAND_FILES}-1)) --export=ALL,WORKSPACE=${WORKSPACE},FAST5_DIR=${FAST5_DIR},FILENAME_PREFIX=${FILENAME_PREFIX},FASTA=${FASTA},ARCH=${ARCH},MODEL=${MODEL} ${HOME}/git/MAFIA/array_basecaller.sh

NUM_ARRAYS=""
for f in ${WORKSPACE}/${DATASET}_fast5_paths_part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}
sbatch --array=${NUM_ARRAYS} --export=ALL,WORKSPACE=${WORKSPACE},FAST5_DIR=${FAST5_DIR},FILENAME_PREFIX=${FILENAME_PREFIX},FASTA=${FASTA},ARCH=${ARCH},MODEL=${MODEL} ${HOME}/git/MAFIA/array_basecaller.sh

#for f in ${FASTA}+([0-9]); do echo $f; grep '>' $f | wc -l; done

cat ${FASTA}+([0-9]) > ${FASTA}_merged

#### align to genome ###
module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF} ${FASTA}_merged > ${SAM}

### check stats and accuracy ###
samtools flagstats ${SAM}
${HOME}/git/renata/accuracy.py ${SAM} ${REF}

#### Convert to BAM and index ###
samtools view -bST ${REF} ${SAM} | samtools sort - > ${BAM}
samtools index ${BAM}

### clean up ###
rm ${FASTA}+([0-9])
rm ${WORKSPACE}/${DATASET}_fast5_paths_all