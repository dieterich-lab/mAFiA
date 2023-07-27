PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian/oligo
cd ${PRJ_DIR}

#############################################################################
### WUE #####################################################################
#############################################################################
#INPUT_RUNS="WUE_batch1 WUE_batch2"
#INPUT_MODS="A m6A"
#OUTPUT_RUN="WUE"
#############################################################################
### ISA #####################################################################
#############################################################################
#INPUT_RUNS="ISA_run1 ISA_run2_1 ISA_run2_2 ISA_run3_1 ISA_run3_2"
#INPUT_MODS="A m6A"
#OUTPUT_RUN="ISA"
#############################################################################
### ISA-WUE #################################################################
#############################################################################
INPUT_RUNS="ISA WUE"
INPUT_MODS="A m6A"
OUTPUT_RUN="ISA-WUE"
#############################################################################
#############################################################################
echo "Merging runs: ${INPUT_RUNS// /, }"

for MOD in $INPUT_MODS
do
  OUTDIR=${PRJ_DIR}/${OUTPUT_RUN}_${MOD}

  ### fast5 ###
  echo "Merging fast5 files"
  mkdir -p ${OUTDIR}/fast5
  cd ${OUTDIR}/fast5
  for RUN in ${INPUT_RUNS}
  do
    for f in ${PRJ_DIR}/${RUN}_${MOD}/fast5/*.fast5
    do
      ln -s ${f}
    done
  done

  ### BAM ###
  echo "Merging BAMs"
  cd ${OUTDIR}
  BAM_INPUTS=""
  BAM_OUTPUT=${OUTDIR}/spomelette_q70.bam
  for RUN in ${INPUT_RUNS}
  do
    BAM_INPUTS+="${PRJ_DIR}/${RUN}_${MOD}/spomelette_q70.bam "
  done
  samtools merge - ${BAM_INPUTS} | samtools sort - > ${BAM_OUTPUT}
  samtools index ${BAM_OUTPUT}

  ### reference ###
  echo "Merging ligation references"
  REF_INPUTS=""
  REF_OUTPUT=${OUTDIR}/ligation_ref.fasta
  for RUN in ${INPUT_RUNS}
  do
    REF_INPUTS+="${PRJ_DIR}/${RUN}_${MOD}/ligation_ref.fasta "
  done
  awk '/^>/{p=seen[$0]++}!p' ${REF_INPUTS} > ${REF_OUTPUT}
done

###############################################################################
### merge A and m6A ###########################################################
###############################################################################
MERGE_DIR=${PRJ_DIR}/${OUTPUT_RUN}_${INPUT_MODS/ /_}
MERGE_BAM_INPUTS=""
MERGE_BAM_OUTPUT=${MERGE_DIR}/spomelette_q70.bam
MERGE_REF_INPUTS=""
MERGE_REF_OUTPUT=${MERGE_DIR}/ligation_ref.fasta
mkdir -p ${MERGE_DIR}/fast5
cd ${MERGE_DIR}/fast5

echo "Merging A and m6A"
for MOD in ${INPUT_MODS}
do
  for f in ${PRJ_DIR}/${OUTPUT_RUN}_${MOD}/fast5/*.fast5
  do
    ln -s ${f}
  done

  MERGE_BAM_INPUTS+="${PRJ_DIR}/${OUTPUT_RUN}_${MOD}/spomelette_q70.bam "

  MERGE_REF_INPUTS+="${PRJ_DIR}/${OUTPUT_RUN}_${MOD}/ligation_ref.fasta "
done

cd ${MERGE_DIR}
samtools merge - ${MERGE_BAM_INPUTS} | samtools sort - > ${MERGE_BAM_OUTPUT}
samtools index ${MERGE_BAM_OUTPUT}

awk '/^>/{p=seen[$0]++}!p' ${MERGE_REF_INPUTS} > ${MERGE_REF_OUTPUT}

echo "Done"