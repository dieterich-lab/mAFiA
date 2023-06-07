WORKSPACE=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
cd ${WORKSPACE}

ORIG=WUE
RUNS="batch1 batch2"
MODS="A m6A"

for MOD in $MODS
do
  OUTDIR=${WORKSPACE}/${ORIG}_${RUNS/ /_}_${MOD}

  ### fast5 ###
  mkdir -p ${OUTDIR}/fast5
  cd ${OUTDIR}/fast5
  for RUN in ${RUNS}
  do
    for f in ${WORKSPACE}/${ORIG}_${RUN}_${MOD}/fast5/*.fast5
    do
      ln -s ${f}
    done
  done

  ### BAM ###
  cd ${OUTDIR}
  BAM_INPUTS=""
  BAM_OUTPUT=${OUTDIR}/spomelette_q70.bam
  for RUN in ${RUNS}
  do
    BAM_INPUTS+="${WORKSPACE}/${ORIG}_${RUN}_${MOD}/spomelette_q70.bam "
  done
  samtools merge ${BAM_OUTPUT} ${BAM_INPUTS}
  samtools index ${BAM_OUTPUT}

  ### reference ###
  REF_INPUTS=""
  REF_OUTPUT=${OUTDIR}/ligation_ref.fasta
  for RUN in ${RUNS}
  do
    REF_INPUTS+="${WORKSPACE}/${ORIG}_${RUN}_${MOD}/ligation_ref.fasta "
  done
  awk '/^>/{p=seen[$0]++}!p' ${REF_INPUTS} > ${REF_OUTPUT}
done