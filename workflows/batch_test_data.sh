#!/usr/bin/env bash

prj=/prj/TRR319_RMaP/Project_BaseCalling

#challenge=1
#for test_ds in {1..4}
challenge=2
for test_ds in {5..6}
do
  source_dir=${prj}/Gm/Test_Datasets/Challenge_${challenge}/Dataset${test_ds}
  dest_dir=${prj}/Adrian/Gm/Test/Dataset${test_ds}

  cd ${dest_dir}

  for i in {0..4}
  do
    mkdir -p ${dest_dir}/fast5_batch${i}
    for ((ii=$(($i*5));ii<=$(($i*5+4));ii++))
    do
      ln -s ${source_dir}/dataset${test_ds}_${ii}.fast5 fast5_batch${i}/dataset${test_ds}_${ii}.fast5
    done
  done
done