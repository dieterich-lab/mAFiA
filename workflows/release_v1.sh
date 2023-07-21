conda activate MAFIA

workspace="/prj/TRR319_RMaP/Project_BaseCalling/Adrian/release_v1"
backbone="${workspace}/models/backbone.torch"
classifier="${workspace}/models/mAFiA"
fast5dir="${workspace}/data/fast5_chrX"
basecall="${workspace}/output/basecall.fasta"
mkdir -p "${workspace}/output"

########################################################################################################################
### basecall ###########################################################################################################
########################################################################################################################
python3 ${HOME}/git/MAFIA/rodan_viterbi.py \
--fast5dir ${fast5dir} \
--model ${backbone} \
--batchsize 2048 \
--decoder viterbi \
--arch ${HOME}/git/MAFIA/rnaarch \
> ${basecall} &

########################################################################################################################
#### align #############################################################################################################
########################################################################################################################
ref="${workspace}/data/GRCh38_96.X.fa"
bam="${workspace}/output/aligned_filtered_q50.bam"

module purge
module load minimap2

minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${ref} ${basecall} \
| samtools view -bST ${ref} -q50 - \
| samtools sort - > ${bam}

samtools index ${bam}

########################################################################################################################
#### mAFiA #############################################################################################################
########################################################################################################################
mod="${workspace}/data/GLORI_chrX.bed"

python3 ${HOME}/git/MAFIA/test_mAFiA.py \
--test_bam_file ${bam} \
--test_fast5_dir ${fast5dir} \
--ref_file ${ref} \
--mod_file ${mod} \
--min_coverage 50 \
--max_num_reads 1000 \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifier} \
--mod_prob_thresh 0.5 \
--out_dir "${workspace}/output" &