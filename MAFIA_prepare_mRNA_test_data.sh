REF=${HOME}/Data/genomes/GRCh38_96.fa
ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

### old HEK293 ######################################################################################################################
#DATASET=HEK293A_WT
#FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all

#DATASET=HEK293_IVT
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293_IVT_2/fast5_pass

### new HEK293 ######################################################################################################################
#DATASET=100_WT_0_IVT_RTA
DATASET=0_WT_100_IVT_RTA
FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230419_HEK293_WT_IVT_Mix/${DATASET}/*/fast5_*

#####################################################################################################################################

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${DATASET}
mkdir -p ${WORKSPACE}

SAM=${WORKSPACE}/mapped.sam
BAM=${WORKSPACE}/filtered_q50.bam

source ${HOME}/git/renata/virtualenv/bin/activate

### basecall large number of reads ###
FILENAME_PREFIX=fast5_paths_part
ls -1 ${FAST5_DIR}/*.fast5 > ${WORKSPACE}/fast5_paths_all
split -l10 -d ${WORKSPACE}/fast5_paths_all ${WORKSPACE}/${FILENAME_PREFIX}

NUM_ARRAYS=""
for f in ${WORKSPACE}/fast5_paths_part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}
sbatch --array=${NUM_ARRAYS} --export=ALL,WORKSPACE=${WORKSPACE},FILENAME_PREFIX=${FILENAME_PREFIX},ARCH=${ARCH},MODEL=${MODEL} ${HOME}/git/MAFIA/array_basecaller.sh

#for f in ${WORKSPACE}/part*.fasta; do echo $f; grep '>' $f | wc -l; done

cat ${WORKSPACE}/part*.fasta > ${WORKSPACE}/merged.fasta

#### align to genome ###
module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF} ${WORKSPACE}/merged.fasta > ${SAM}

### check stats and accuracy ###
samtools flagstats ${SAM}
${HOME}/git/renata/accuracy.py ${SAM} ${REF}

#### Convert to BAM and index ###
samtools view -bST ${REF} -q50 ${SAM} | samtools sort - > ${BAM}
samtools index ${BAM}

### clean up ###
rm ${WORKSPACE}/part*.fasta
rm ${WORKSPACE}/fast5_paths_all