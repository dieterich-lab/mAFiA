sterm -c 4 -m 120GB -N gpu-g3-1

BACKBONE_MODEL=/home/achan/pytorch_models/tRNA_IVT/tRNA_IVT-epoch29.torch
ARCH=/home/achan/git/renata/rnaarch

WORKSPACE=/beegfs/prj/tRNA_Berlin/newBatchDec2022_Spombe

#DATASET=tRNA_IVT
DATASET=tRNA_Q
READS=${WORKSPACE}/agg_fast5_pass/${DATASET}
FASTA=${WORKSPACE}/achan/basecall/${DATASET}.fasta
RAND_REF=${WORKSPACE}/rand_ref/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta.rand
REF=${WORKSPACE}/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta
SAM=${WORKSPACE}/achan/mapping/parasail/${DATASET}.sam
CSV=${SAM//.sam/.csv}
BAM=${SAM//.sam/.bam}

mkdir -p `dirname ${FASTA}`
mkdir -p `dirname ${SAM}`
source ${HOME}/git/renata/virtualenv/bin/activate

### basecall ###
python3 ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${READS} \
--arch ${ARCH} \
--model ${BACKBONE_MODEL} \
> ${FASTA}

### align and check accuracy ###
#module load minimap2
#minimap2 --secondary=no -ax map-ont -k 8 -t 36 --cs ${REF} ${FASTA} > ${SAM}
#${HOME}/git/renata/accuracy.py ${SAM} ${REF}

### align with parasail ###
module load parasail/2.5
parasail_aligner -x -a sw -f ${REF} -q ${FASTA} -g ${CSV}
parasail_aligner -x -a sw_trace -f ${REF} -q ${FASTA} -O SAMH -g ${SAM}

### random alignment ###
for i in {0..99}; do parasail_aligner -x -a sw -f ${RAND_REF}$i -q ${FASTA} -g ${CSV}.rand$i; done

### filter ###
python3 /home/achan/git/MAFIA/filter_parasail_mapping.py \
--fasta_file ${FASTA} \
--ref_file ${REF} \
--csv_file ${CSV} \
--sam_file ${SAM}

#### Convert to BAM and index ###
samtools view -bST ${REF} ${SAM}.filtered | samtools sort - > ${BAM}.filtered.sorted
samtools index ${BAM}.filtered.sorted

### clean up ###
rm ${CSV}.rand*
