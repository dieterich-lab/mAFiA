#!/usr/bin/env bash

#DATASET=100_WT_0_IVT
DATASET=75_WT_25_IVT
#DATASET=50_WT_50_IVT
#DATASET=25_WT_75_IVT
#DATASET=0_WT_100_IVT

echo ${DATASET}

PRJ=/prj/TRR319_RMaP/Project_BaseCalling/Adrian
SCRATCH=/scratch/achan
NANOPOLISH=${SCRATCH}/nanopolish/${DATASET}
FAST5=${PRJ}/HEK293/${DATASET}/fast5
FASTQ=${SCRATCH}/nanopolish/${DATASET}/fastq
REF=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.cdna.all.fa
BAM=${NANOPOLISH}/minimap.bam
m6ANET_OUTDIR=${SCRATCH}/m6Anet/${DATASET}
CHEUI_OUTDIR=${SCRATCH}/CHEUI/${DATASET}
GIT_CHEUI=/home/achan/git/CHEUI

######################################################################################
### copy fastq #######################################################################
######################################################################################
#mkdir -p ${FASTQ} && cd ${FASTQ}
#
#for f in ${FAST5}/*.fast5
#do
#  f5_path=$(realpath $f)
#  fastq_path=${f5_path//fast5/fastq}.gz
#  ln -s ${fastq_path}
#
#  parent_dir=$(dirname $(dirname $f5_path))
#  seq_summary_name=$(basename $(echo ${parent_dir}/sequencing_summary*.txt))
#  if [ ! -f "${NANOPOLISH}/${seq_summary_name}" ]; then
#    cp ${parent_dir}/sequencing_summary*.txt ${NANOPOLISH}
#  fi
#done
#
#cd ${NANOPOLISH}
#cat ${FASTQ}/*.fastq.gz > all.fastq.gz
#awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' sequencing_summary_*.txt > sequencing_summary.txt
#rm sequencing_summary_*.txt

######################################################################################
### nanopolish #######################################################################
######################################################################################
#module load minimap2
#module load vbz_compression nanopolish
#
#### index reads ###
#echo "Indexing reads"
#nanopolish index -s sequencing_summary.txt -d ${FAST5} all.fastq.gz
#
#### align with minimap ###
#minimap2 -ax map-ont --secondary=no -t 8 ${REF} all.fastq.gz | samtools view -F 2324 -b - | samtools sort -o ${BAM}
#samtools index ${BAM}
#
#### eventalign ###
#echo "Eventalign..."
#nanopolish eventalign \
#--reads all.fastq.gz \
#--bam ${BAM} \
#--genome ${REF} \
#--scale-events \
#--signal-index \
#--samples \
#--summary summary.txt \
#--threads 48 \
#> eventalign.txt

######################################################################################
### m6Anet ###########################################################################
######################################################################################
#module load cuda m6anet
#
#mkdir -p ${m6ANET_OUTDIR}
#cd ${m6ANET_OUTDIR}
#
#cut -d "    " -f 1-15 ${NANOPOLISH}/eventalign.txt > eventalign.txt
#echo "m6Anet dataprep...":
#m6anet dataprep --eventalign eventalign.txt --out_dir ${m6ANET_OUTDIR} --n_processes 4
#echo "m6Anet inference...":
#m6anet inference --input_dir ${m6ANET_OUTDIR} --out_dir ${m6ANET_OUTDIR} --n_processes 4 --num_iterations 1000

######################################################################################
### CHEUI ############################################################################
######################################################################################
eval "$(conda shell.bash hook)"
conda activate CHEUI

#mkdir -p ${CHEUI_OUTDIR}
#cd ${CHEUI_OUTDIR}
echo "Splitting nanopolish file..."
srun -p general \
python3 /home/achan/git/MAFIA/workflows/split_large_nanopolish_File.py \
--infile ${NANOPOLISH}/eventalign.txt \
--outfile_prefix ${CHEUI_OUTDIR}/eventalign_part &

num_parts=$(ls ${CHEUI_OUTDIR}/eventalign_part*.txt | wc -l)

cd ${GIT_CHEUI}/scripts/preprocessing_CPP
for i in {0..3}
do
  PART=$(printf %02d $i)
  echo "CHEUI preprocessing part${PART}..."
  srun -p general -c 4 --mem 128GB \
  ./CHEUI \
  -i ${CHEUI_OUTDIR}/eventalign_part${PART}.txt \
  -m ${GIT_CHEUI}/kmer_models/model_kmer.csv \
  -o ${CHEUI_OUTDIR}/prep_m6A_part${PART} \
  -n 16 \
  --m6A &
done

#wait

for i in {0..$((num_parts-1))}
do
  PART=$(printf %02d $i)
  echo "CHEUI predict model 1 on part${PART}..."
  srun --partition=gpu --cpus-per-task=8 \
  python3 ${GIT_CHEUI}/scripts/CHEUI_predict_model1.py \
  -m ${GIT_CHEUI}/CHEUI_trained_models/CHEUI_m6A_model1.h5 \
  -i ${CHEUI_OUTDIR}/prep_m6A/eventalign_part${PART}_signals+IDS.p \
  -o ${CHEUI_OUTDIR}/read_level_m6A_predictions_part${PART}.txt \
  -l ${DATASET} &
done

cat ${CHEUI_OUTDIR}/read_level_m6A_predictions_part00.txt > ${CHEUI_OUTDIR}/read_level_m6A_predictions_combined.txt
for i in {1..3}
do
  PART=$(printf %02d $i)
  cat ${CHEUI_OUTDIR}/read_level_m6A_predictions_part${PART}.txt >> ${CHEUI_OUTDIR}/read_level_m6A_predictions_combined.txt
done

echo "Sort reads..."
srun -p general -c 8 --mem 64GB \
sort -k1  --parallel=16 ${CHEUI_OUTDIR}/read_level_m6A_predictions_combined.txt > ${CHEUI_OUTDIR}/read_level_m6A_predictions_sorted.txt

echo "CHEUI predict model 2..."
srun --partition=gpu --cpus-per-task=16 \
python3 ${GIT_CHEUI}/scripts/CHEUI_predict_model2.py \
-m  ${GIT_CHEUI}/CHEUI_trained_models/CHEUI_m6A_model2.h5 \
-i ${CHEUI_OUTDIR}/read_level_m6A_predictions_sorted.txt \
-o ${CHEUI_OUTDIR}/site_level_m6A_predictions.txt &

tail -n +2 ${CHEUI_OUTDIR}/site_level_m6A_predictions.txt | split -a 2 -d -l 100000 - ${CHEUI_OUTDIR}/site_level_m6A_predictions.txt.part
for file in ${CHEUI_OUTDIR}/site_level_m6A_predictions.txt.part*
do
    head -n 1 ${CHEUI_OUTDIR}/site_level_m6A_predictions.txt > with_header_tmp
    cat "$file" >> with_header_tmp
    mv -f with_header_tmp "$file"
done

for file in ${CHEUI_OUTDIR}/site_level_m6A_predictions.txt.part*
do
  srun python3 ${HOME}/git/MAFIA/misc/convert_CHEUI_results_to_GLORI_gDNA_coords.py ${file} &
done

cat site_level_m6A_predictions.txt.part00.glori > site_level_m6A_predictions.txt.glori
for i in {1..19}
do
  part=$(printf %02d $i)
  tail -n+2 site_level_m6A_predictions.txt.part${part}.glori >> site_level_m6A_predictions.txt.glori
done

#cat site_level_m6A_predictions.txt.part00.glori.filtered > site_level_m6A_predictions.txt.glori.filtered
#for i in {1..13}
#do
#  part=$(printf %02d $i)
#  tail -n+2 site_level_m6A_predictions.txt.part${part}.glori.filtered >> site_level_m6A_predictions.txt.glori.filtered
#done
