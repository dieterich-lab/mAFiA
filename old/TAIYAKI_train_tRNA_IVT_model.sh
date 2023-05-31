sterm -c 4 -m 120GB -N gpu-g3-1

BASE_DIR=/prj/tRNA_Berlin/newBatchDec2022_Spombe
REF=/prj/tRNA_Berlin/newBatchDec2022_Spombe/S.pombe_mature_tRNAs_adapters_nuclear_and_mt.fasta
READS=/prj/tRNA_Berlin/newBatchDec2022_Spombe/agg_fast5_pass/4tRNAs_IVT
BAM=/prj/tRNA_Berlin/newBatchDec2022_Spombe/agg_pass_bams/4tRNAs_IVT.bam

export PYTHONPATH="/home/achan/git/taiyaki"
MODEL=/home/achan/git/taiyaki/models/r941_rna_minion.checkpoint
TAIYAKI=${BASE_DIR}/taiyaki

source ${HOME}/.venv/taiyaki/bin/activate
mkdir -p ${TAIYAKI}

### get read references ###
python3 ${HOME}/git/taiyaki/bin/get_refs_from_sam.py ${REF} ${BAM} \
--min_coverage 0.8 \
--reverse \
> ${TAIYAKI}/read_references.fasta

### per-read params ###
python3 ${HOME}/git/taiyaki/bin/generate_per_read_params.py --jobs 36 ${READS} > ${TAIYAKI}/read_params.tsv

### prepare mapped reads ###
python3 ${HOME}/git/taiyaki/bin/prepare_mapped_reads.py \
${READS} \
${TAIYAKI}/read_params.tsv \
${TAIYAKI}/mapped_reads.hdf5 \
${MODEL} \
${TAIYAKI}/read_references.fasta \
--jobs 36

### convert hdf5 ###
${HOME}/git/taiyaki/misc/merge_mappedsignalfiles.py \
--input ${TAIYAKI}/mapped_reads.hdf5 None \
${TAIYAKI}/merged.hdf5

### generate rodan chunks ###
deactivate
source ${HOME}/git/renata/virtualenv/bin/activate

${HOME}/git/renata/gendata.py \
--infile ${TAIYAKI}/merged.hdf5 \
--outdir ${TAIYAKI} \
--seqlen 1024 \
--maxchunks 35000 \
--maxvalchunks 4000

### train ###
${HOME}/git/renata/train.py \
--config /home/achan/pytorch_models/tRNA_IVT/tRNA_IVT.config \
--name tRNA_IVT \
--savedir /home/achan/pytorch_models/tRNA_IVT \
--workers 4
