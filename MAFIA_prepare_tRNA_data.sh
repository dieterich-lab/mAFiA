BACKBONE_MODEL=/home/achan/pytorch_models/tRNA_IVT/tRNA_IVT-epoch29.torch
ARCH=/home/achan/git/renata/rnaarch

WORKSPACE=/beegfs/prj/tRNA_Berlin/newBatchDec2022_Spombe

#DATASET=tRNA_IVT
#DATASET=tRNA_Q
#DATASET=AEP1_T1
#DATASET=AEP1_T2
DATASET=AEP1_T3
#DATASET=AEP1_Q_T1
#DATASET=AEP1_Q_T2
#DATASET=AEP1_Q_T3

READS=${WORKSPACE}/agg_fast5_pass/${DATASET}
FASTA=${WORKSPACE}/achan/basecall/${DATASET}.fasta
BAM=${WORKSPACE}/achan/mapping/${DATASET}.bam
REF=${WORKSPACE}/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta

mkdir -p `dirname ${FASTA}`
mkdir -p `dirname ${BAM}`
source ${HOME}/git/renata/virtualenv/bin/activate

### basecall ###
srun --partition=gpu --exclude=gpu-g4-1 \
python3 ${HOME}/git/renata/basecall_viterbi.py \
--fast5dir ${READS} \
--arch ${ARCH} \
--model ${BACKBONE_MODEL} \
> ${FASTA}

####################################
### Christoph's mapping workflow ###
####################################
export SCRATCH="/scratch/achan/"`basename ${FASTA} | cut -d'.' -f1`
export SCRIPTS="/prj/tRNA_Berlin/newBatchDec2022_Spombe/achan/alignment_workflow/scripts_fasta"

mkdir ${SCRATCH}
grep -Ev 'Using|Finish' ${FASTA} | grep -E '>|A|C|G|T' > ${SCRATCH}/`basename $FASTA`
export CP_FASTA=${SCRATCH}/`basename $FASTA`

#split input and compute alignments (forward and reverse)
${SCRIPTS}/workflow_tRNA_fasta.sh ${CP_FASTA} ${REF}

#merge BAM files and keep only best alignments with minimal score of 20
#Wait for all slurm jobss to finish first (from previous step)
${SCRIPTS}/workflow_tRNA_fasta_step2.sh ${CP_FASTA} ${REF}

#determine score cutoff and filter alignments
srun -p small ${SCRIPTS}/tRNA_alignment_score_cutoff.sh ${CP_FASTA}

cp ${SCRATCH}/final_output_pass_filter.bam ${BAM}
cp ${SCRATCH}/final_output_pass_filter.bam.bai ${BAM}.bai