#!/usr/bin/env bash
shopt -s globstar

MODEL="${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch"
#REF_GENOME="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"
REF_GENOME="/biodb/genomes/mus_musculus/GRCm38_102/GRCm38_102.fa"
#WORKSPACE="/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/Diet/D2_M3KO_WD"
WORKSPACE="/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/${ds}"

### new HEK293 #########################################################################################################
#DATASET=0_WT_100_IVT
#DATASET=25_WT_75_IVT
#DATASET=50_WT_50_IVT
#DATASET=75_WT_25_IVT
#DATASET=100_WT_0_IVT
#DATASET=P2_WT
#DATASET=Mettl3-KO

#WORKSPACE=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/${DATASET}

########################################################################################################################
#DATASET=A1_WT_CD
#DATASET=A2_WT_CD
#WORKSPACE=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/Federica_Accornero/${DATASET}

########################################################################################################################

#DATASET=JK_HEK293_DMSO_1_2_RTA
#DATASET=JK_HEK293_DMSO_3_4_RTA
#DATASET=JK_HEK293_STM2457_5_6_RTA
#DATASET=JK_HEK293_STM2457_7_8_RTA
#DATASET=JK_HEK293_DMSO_merged

#DATASET=mESC_WT_DMSO_merged
#DATASET=mESC_Mettl3_KO_merged
#DATASET=mESC_WT_STM_merged

#WORKSPACE=/prj/TRR319_RMaP/Project_B01/Adrian/${DATASET}

########################################################################################################################
#DATASET=40-26
#DATASET=40-27
#DATASET=40-28
#DATASET=40-29
#DATASET=40-30
#DATASET=40-31
#DATASET=40-32
#DATASET=40-33
#DATASET=40-34

#WORKSPACE=/prj/Dewenter_TAC_Backs_lab/achan/${DATASET}
#DATA_DIR=/prj/Dewenter_TAC_Backs_lab/raw_data/Nanopore_dRNA/Cologne
########################################################################################################################

#WORKSPACE=/prj/TRR319_RMaP_BaseCalling/Adrian/HeLa_SRR28796313

########################################################################################################################
#ds=HEK_siCtrl_input_rep1
#ds=HEK_siMETTL3_input_rep1
#ds=HEK_siTRUB1_input_rep1

#ds=HEK_siCtrl_input_rep2
#ds=HEK_siMETTL3_input_rep2
#ds=HEK_siTRUB1_input_rep2

#WORKSPACE=/prj/TRR319_RMaP_BaseCalling/Adrian/NanoSPA/${ds}
########################################################################################################################

#mkdir -p ${WORKSPACE}
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
