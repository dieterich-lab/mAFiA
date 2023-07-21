# m6AFiA - (Another) m6A Finding Algorithm

Brief demo on chr X:

Preliminary:
1. Download models and data from https://data.dieterichlab.org/index.php/apps/files?dir=/Dieterich%20lab%20group%20folder/mAFiA_releave_v1
2. The folder "models" contains:
- backbone.torch: RODAN-based neural network for basecalling and feature extraction
- backbone.config: cfg file for backbone
- mAFiA: pickled models for mAFiA
3. The folder "data" contains:
  - fast5_chrX: dRNA-Seq raw data from HEK293 WT mRNA, chr X only
  - GRCh38_96.X: genome reference
  - GLORI_chrX.bed: GLORI mod-sites on chr X, bed format
4. Assume that data and model are unzipped to ${workspace}
  backbone="${workspace}/models/backbone.torch"
  classifier="${workspace}/models/mAFiA"
  fast5dir="${workspace}/data/fast5_chrX"
  ref="${workspace}/data/GRCh38_96.X.fa"
  mod="${workspace}/data/GLORI_chrX.bed"

  mkdir -p "${workspace}/output"
  basecall="${workspace}/output/basecall.fasta"
  bam="${workspace}/output/aligned_filtered_q50.bam"

5. Activate virtual environment
   conda activate MAFIA

Basecalling:
python3 ${HOME}/git/MAFIA/rodan_viterbi.py \
--fast5dir ${fast5dir} \
--model ${backbone} \
--batchsize 2048 \
--decoder viterbi \
--arch ${HOME}/git/MAFIA/rnaarch \
> ${basecall} &

Should take about 20 mins on a GPU computer

Alignment:
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${ref} ${basecall} \
| samtools view -bST ${ref} -q50 - \
| samtools sort - > ${bam}

samtools index ${bam}

m6A Detection:
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

Should take less than 1 hour. We are currently working on combined the feature extraction step with basecalling. The run-time should then be significantly compressed.
