#!/usr/bin/env bash
shopt -s globstar

ds="NTERA/rep2"

MODEL="${HOME}/git/mAFiA_dev/models/RODAN_HEK293_IVT.torch"

WORKSPACE="/prj/TRR319_RMaP_BaseCalling/Adrian/${ds}"
REF_GENOME="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"

#WORKSPACE="/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/Diet/D2_M3KO_WD"
#REF_GENOME="/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa"

cd ${WORKSPACE}

########################################################################################################################
### basecall large number of reads #####################################################################################
########################################################################################################################
source ${HOME}/git/mAFiA/mafia-venv/bin/activate

fast5_dir=${WORKSPACE}/fast5
FILENAME_PREFIX=fast5_paths_part
ls -1 ${fast5_dir}/*.fast5 > ${WORKSPACE}/fast5_paths_all
split -a3 -l10 -d ${WORKSPACE}/fast5_paths_all ${WORKSPACE}/${FILENAME_PREFIX}

echo "Basecalling ${fast5_dir}..."

num_arrays=""
for f in ${WORKSPACE}/fast5_paths_part*; do ff=${f##*part}; num_arrays+="${ff},"; done
num_arrays=${num_arrays%,*}
jid=$(sbatch --parsable --array=${num_arrays} --export=ALL,WORKSPACE=${WORKSPACE},FILENAME_PREFIX=${FILENAME_PREFIX},MODEL=${MODEL} ${HOME}/git/mAFiA_dev/workflows/array_basecaller.sh)

sbatch --dependency=afterok:"${jid}" --export=ALL,WORKSPACE=${WORKSPACE},REF_GENOME=${REF_GENOME} ${HOME}/git/mAFiA_dev/workflows/prepare_mRNA_dataset_part2.sh
