### old HEK293 ######################################################################################################################
#TEST_DATASET=HEK293A_WT
#FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all

#TEST_DATASET=HEK293_IVT
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293_IVT_2/fast5_pass

### new HEK293 ######################################################################################################################
#TEST_DATASET=0_WT_100_IVT_RTA
#TEST_DATASET=25_WT_75_IVT_RTA
#TEST_DATASET=100_WT_0_IVT_RTA
#FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230419_HEK293_WT_IVT_Mix/${TEST_DATASET}/*/fast5_*

#TEST_DATASET=50_WT_50_IVT_RTA
TEST_DATASET=75_WT_25_IVT_RTA
FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Isabel/20230426_HEK293_WT_IVT/${TEST_DATASET}/*/fast5_*

#####################################################################################################################################

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${TEST_DATASET}
REF=${HOME}/Data/genomes/GRCh38_96.fa
BAM=${WORKSPACE}/filtered_q50.bam
MOD_FILE=${HOME}/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXTRACTION_LAYER=convlayers.conv21
TRAIN_DATASET=20230221_WUE_splint_lig
CLASSIFIER=logistic_regression
SCALER=MaxAbs
CLASSIFIER_MODEL_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/MAFIA_classifiers/${TRAIN_DATASET}_${CLASSIFIER}_${SCALER}
OUTFILE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_${TEST_DATASET}_${TRAIN_DATASET}_modProbPerRead.tsv

NUM_ARRAYS=""
for f in ${MOD_FILE}.part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}

sbatch --array=${NUM_ARRAYS} \
--export=ALL,\
WORKSPACE=${WORKSPACE},\
FAST5_DIR=${FAST5_DIR},\
REF=${REF},\
BAM=${BAM},\
MOD_FILE=${MOD_FILE},\
BACKBONE_MODEL=${BACKBONE_MODEL},\
EXTRACTION_LAYER=${EXTRACTION_LAYER},\
CLASSIFIER_MODEL_DIR=${CLASSIFIER_MODEL_DIR},\
OUTFILE=${OUTFILE} \
${HOME}/git/MAFIA/MAFIA_SWARM_test_binary_classifier_mRNA.sh

### concat output ###
cp ${OUTFILE}.part00 ${OUTFILE}.merged
for num in ${NUM_ARRAYS//,/ }
do
  if [ $num != '00' ]
  then awk NR\>1 $OUTFILE.part$num >> ${OUTFILE}.merged
  fi
  done

rm ${OUTFILE}.part*