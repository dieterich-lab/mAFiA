#####################################################################################################################################
#DATASET=HEK293A_WT
#FAST5_DIR=/beegfs/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/fast5_all

DATASET=HEK293_IVT
FAST5_DIR=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/fast5/HEK293_IVT_2/fast5_pass

#####################################################################################################################################

WORKSPACE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/${DATASET}
REF=${HOME}/Data/genomes/GRCh38_96.fa
BAM=${WORKSPACE}/genome_mapped_q50.bam
MOD_FILE=${HOME}/Data/GLORI/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv
BACKBONE_MODEL=${HOME}/pytorch_models/HEK293_IVT_2_q50_10M/HEK293_IVT_2_q50_10M-epoch29.torch
EXTRACTION_LAYER=convlayers.conv21
TRAINING_DATA=random_ligation_A_m6A
CLASSIFIER_MODEL_DIR=${WORKSPACE}/MAFIA_classifiers/${TRAINING_DATA}
MOD_PROB_THRESH=0.5
OUTFILE=/beegfs/prj/TRR319_RMaP/Project_BaseCalling/Adrian/results/res_${DATASET}_${TRAINING_DATA}_modProbPerRead.tsv

NUM_ARRAYS=""
for f in ${MOD_FILE}.part*; do ff=${f##*part}; NUM_ARRAYS+="${ff},"; done
NUM_ARRAYS=${NUM_ARRAYS%,*}
sbatch --array=${NUM_ARRAYS} --export=ALL,WORKSPACE=${WORKSPACE},BAM=${BAM},MOD_FILE=${MOD_FILE},FAST5_DIR=${FAST5_DIR},REF=${REF},BACKBONE_MODEL=${BACKBONE_MODEL},EXTRACTION_LAYER=${EXTRACTION_LAYER},REF=${REF},CLASSIFIER_MODEL_DIR=${CLASSIFIER_MODEL_DIR},MOD_PROB_THRESH=${MOD_PROB_THRESH},OUTFILE=${OUTFILE} ${HOME}/git/MAFIA/MAFIA_SWARM_test_binary_classifier_mRNA.sh
