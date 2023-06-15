REF=${HOME}/Data/genomes/GRCh38_96.fa
ARCH=${HOME}/git/renata/rnaarch
MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch

### new HEK293 ######################################################################################################################
#DATASET=0_WT_100_IVT
#DATASET=25_WT_75_IVT
#DATASET=50_WT_50_IVT
#DATASET=75_WT_25_IVT
#DATASET=100_WT_0_IVT

DATASET=P2_WT
#####################################################################################################################################

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/${DATASET}
mkdir -p ${WORKSPACE}
cd ${WORKSPACE}

FAST5_DIR=${WORKSPACE}/fast5
mkdir -p ${FAST5_DIR}

#if test -f basecall_merged.fasta
#then
#  cp basecall_merged.fasta basecall_merged.fasta.backup
#  mkdir _old
#  mv basecall_merged.fasta.backup filtered_q50.bam filtered_q50.bam.bai mapped.sam _old
#fi

SAM=${WORKSPACE}/mapped.sam
BAM=${WORKSPACE}/filtered_q50.bam

source ${HOME}/git/renata/virtualenv/bin/activate

### basecall large number of reads ###
FILENAME_PREFIX=fast5_paths_part
ls -1 ${FAST5_DIR}/*.fast5 > ${WORKSPACE}/fast5_paths_all
split -a3 -l10 -d ${WORKSPACE}/fast5_paths_all ${WORKSPACE}/${FILENAME_PREFIX}

NUM_ARRAYS=""
for f in ${WORKSPACE}/fast5_paths_part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}
sbatch --array=${NUM_ARRAYS} --export=ALL,WORKSPACE=${WORKSPACE},FILENAME_PREFIX=${FILENAME_PREFIX},ARCH=${ARCH},MODEL=${MODEL} ${HOME}/git/MAFIA/workflows/array_basecaller.sh

#for f in ${WORKSPACE}/part*.fasta; do echo $f; grep '>' $f | wc -l; done

if test -f basecall_merged.fasta
then
  cat ${WORKSPACE}/part*.fasta >> ${WORKSPACE}/basecall_merged.fasta
else
  cat ${WORKSPACE}/part*.fasta > ${WORKSPACE}/basecall_merged.fasta
fi

#### align to genome ###
module purge
module load minimap2
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${REF} ${WORKSPACE}/basecall_merged.fasta > ${SAM}

### check stats and accuracy ###
samtools flagstats ${SAM}
${HOME}/git/renata/accuracy.py ${SAM} ${REF}

#### Convert to BAM and index ###
samtools view -bST ${REF} -q50 ${SAM} | samtools sort - > ${BAM}
samtools index ${BAM}

### clean up ###
rm ${WORKSPACE}/part*.fasta
rm ${WORKSPACE}/fast5_paths_all