#!/bin/sh

annot_bed="/home/adrian/Data/genomes/mus_musculus/mm10/mm10_annot.sorted.bed"
region_sizes="/home/adrian/Data/genomes/mus_musculus/mm10/region_sizes.txt"
metaPlotR="/home/adrian/git/metaPlotR"

########################################################################################################################
ds="HFpEF"
conditions="ctrl HFpEF"

#ds="Diet"
#conditions="WT_CD WT_WD"

#ds="TAC"
#conditions="SHAM TAC"

#res_dir="/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR"
########################################################################################################################

#res_dir="/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/polyA"
res_dir="/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR"

#polyAs="below_50 above_150"
thresholds="0.0 50.0"

for this_cond in ${conditions}
do
  file_prefix="${ds}_${this_cond}"
  for thresh_modRatio in ${thresholds}
  do
    for this_mod in "m6A" "psi"
    do
      in_bed="${res_dir}/${file_prefix}_${this_mod}_modRatio${thresh_modRatio}.bed"
      out_bed="${res_dir}/${file_prefix}_${this_mod}_modRatio${thresh_modRatio}.dist.measures.txt"

      sort -k1,1 -k2,2n ${in_bed} > ${in_bed}.sorted
      perl ${metaPlotR}/annotate_bed_file.pl --bed ${in_bed}.sorted --bed2 ${annot_bed} > ${in_bed}.sorted.annot
      perl ${metaPlotR}/rel_and_abs_dist_calc.pl --bed ${in_bed}.sorted.annot --regions ${region_sizes} > ${out_bed}

      rm ${in_bed}.sorted ${in_bed}.sorted.annot
    done
  done
done