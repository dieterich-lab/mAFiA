#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup_bambu
#SBATCH --output=/home/achan/slurm/pileup_bambu_chr%a.out

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
chr="X"
else
chr=${SLURM_ARRAY_TASK_ID}
fi

workspace="/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT/chr${chr}"
out_dir=${workspace}/transcript

echo "Running bambu..."
Rscript /home/achan/git/mAFiA_dev/workflows/run_bambu_on_mAFiA_bam.R ${workspace} ${out_dir}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/homo_sapiens/GRCh38_102/m6A.psi.GRCh38_102.chr${chr}.bed

for ids in ${out_dir}/*.ids
do
  bam=${ids//.ids/.bam}
  output=${out_dir}

  temp=$(basename $ids)
  tx_name=${temp//.ids/}
  echo ${tx_name}

  samtools view -@ 36 -N ${ids} -o ${bam} ${workspace}/mAFiA.reads.bam
  samtools index ${bam}

  python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
    --bam_file ${bam} \
    --mod_file ${mod} \
    --min_coverage 10 \
    --out_dir ${output} \
    --out_filename ${tx_name}.bed \
    --num_jobs 36
done
