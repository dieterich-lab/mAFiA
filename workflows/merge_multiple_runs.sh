PRJ_DIR=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
cd ${PRJ_DIR}

#ORIG=WUE
#ALL_RUNS="batch1 batch2"
#ALL_MODS="A m6A"

ORIG=ISA
ALL_RUNS="run3_1 run3_2"
ALL_MODS="A m6A"

for MOD in $ALL_MODS
do
  OUTDIR=${PRJ_DIR}/${ORIG}_${ALL_RUNS/ /_}_${MOD}

  ### fast5 ###
  mkdir -p ${OUTDIR}/fast5
  cd ${OUTDIR}/fast5
  for RUN in ${ALL_RUNS}
  do
    for f in ${PRJ_DIR}/${ORIG}_${RUN}_${MOD}/fast5/*.fast5
    do
      ln -s ${f}
    done
  done

  ### BAM ###
  cd ${OUTDIR}
  BAM_INPUTS=""
  BAM_OUTPUT=${OUTDIR}/spomelette_q70.bam
  for RUN in ${ALL_RUNS}
  do
    BAM_INPUTS+="${PRJ_DIR}/${ORIG}_${RUN}_${MOD}/spomelette_q70.bam "
  done
  samtools merge ${BAM_OUTPUT} ${BAM_INPUTS}
  samtools index ${BAM_OUTPUT}

  ### reference ###
  REF_INPUTS=""
  REF_OUTPUT=${OUTDIR}/ligation_ref.fasta
  for RUN in ${ALL_RUNS}
  do
    REF_INPUTS+="${PRJ_DIR}/${ORIG}_${RUN}_${MOD}/ligation_ref.fasta "
  done
  awk '/^>/{p=seen[$0]++}!p' ${REF_INPUTS} > ${REF_OUTPUT}
done

###############################################################################
### merge A and m6A ###########################################################
###############################################################################
#MERGE_DIR=${PRJ_DIR}/${ORIG}_${ALL_RUNS/ /_}_${ALL_MODS/ /_}
#MERGE_BAM_INPUTS=""
#MERGE_BAM_OUTPUT=${MERGE_DIR}/spomelette_q70.bam
#MERGE_REF_INPUTS=""
#MERGE_REF_OUTPUT=${MERGE_DIR}/ligation_ref.fasta
#mkdir -p ${MERGE_DIR}/fast5
#cd ${MERGE_DIR}/fast5
#
#for MOD in ${ALL_MODS}
#do
#  ### fast5 ###
#  for f in ${PRJ_DIR}/${ORIG}_${ALL_RUNS/ /_}_${MOD}/fast5/*.fast5
#  do
#    ln -s ${f}
#  done
#
#  MERGE_BAM_INPUTS+="${PRJ_DIR}/${ORIG}_${ALL_RUNS/ /_}_${MOD}/spomelette_q70.bam "
#
#  MERGE_REF_INPUTS+="${PRJ_DIR}/${ORIG}_${ALL_RUNS/ /_}_${MOD}/ligation_ref.fasta "
#done
#
#cd ${MERGE_DIR}
#samtools merge - ${MERGE_BAM_INPUTS} | samtools sort - > ${MERGE_BAM_OUTPUT}
#samtools index ${MERGE_BAM_OUTPUT}
#
#awk '/^>/{p=seen[$0]++}!p' ${MERGE_REF_INPUTS} > ${MERGE_REF_OUTPUT}