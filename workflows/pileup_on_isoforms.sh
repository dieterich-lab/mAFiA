#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup_Fhl1
#SBATCH --output=/home/achan/slurm/pileup_Fhl1_%A.out

workspace="/beegfs/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/${ds}/bambu"
out_dir=${workspace}

gene="Fhl1"
chr="X"
isoforms="ENSMUST00000023854 ENSMUST00000114772 unclassified.new_exon"

cd ${workspace}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.chr${chr}.bed

for isoform in ${isoforms}
do
  bam="${workspace}/${gene}.${isoform}.bam"
  samtools index ${bam}
  python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
    --bam_file ${bam} \
    --mod_file ${mod} \
    --min_coverage 1 \
    --out_dir ${out_dir} \
    --out_filename "${gene}.${isoform}.bed" \
    --num_jobs 36
done

echo "Finished"
