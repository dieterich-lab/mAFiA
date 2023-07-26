conda activate MAFIA

workspace="/prj/TRR319_RMaP/Project_BaseCalling/Adrian/release_v1"
models="${workspace}/models"
data="${workspace}/data"
output="/prj/TRR319_RMaP/Project_BaseCalling/Adrian/output_101"

backbone="${models}/backbone.torch"
classifier="${models}/mAFiA"
fast5dir="${data}/fast5_chrX"
ref="${data}/GRCh38_96.X.fa"
mod="${data}/GLORI_chrX.bed"

mkdir -p "${output}"
basecall="${output}/rodan.fasta"
bam="${output}/minimap.q50.bam"

########################################################################################################################
### basecall ###########################################################################################################
########################################################################################################################
srun --job-name=basecall --partition=gpu --exclude=gpu-g4-1,gpu-g2-1,gpu-g1-1,gpu-g1-2 -c 16 --mem 256GB \
python3 ${HOME}/git/MAFIA/basecall_with_feature_extraction.py \
--fast5dir ${fast5dir} \
--model ${backbone} \
--batchsize 4096 \
--outdir ${output} &

########################################################################################################################
#### align #############################################################################################################
########################################################################################################################
module purge
module load minimap2

minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${ref} ${basecall} \
| samtools view -bST ${ref} -q50 - \
| samtools sort - > ${bam}

samtools index ${bam}

########################################################################################################################
#### mAFiA #############################################################################################################
########################################################################################################################
srun --job-name=mafia --partition=gpu --exclude=gpu-g4-1,gpu-g2-1,gpu-g1-1,gpu-g1-2 -c 4 --mem 128GB \
python3 ${HOME}/git/MAFIA/test_mAFiA.py \
--in_bam_file ${bam} \
--in_feat_file ${output}/features.h5 \
--ref_file ${ref} \
--mod_file ${mod} \
--min_coverage 50 \
--max_num_reads 1000 \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifier} \
--mod_prob_thresh 0.5 \
--out_dir ${output} &