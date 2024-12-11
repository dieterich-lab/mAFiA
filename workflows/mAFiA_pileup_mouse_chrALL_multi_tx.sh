#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=90GB
#SBATCH --verbose
#SBATCH --job-name=pileup_mouse_chrALL_multi_tx
#SBATCH --output=/home/achan/slurm/pileup_mouse_chrALL_multi_tx_%A.out

ds="TAC"
cond="TAC"

workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/DTU/${ds}
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.bed
output=${workspace}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate

for tx in `cat ${workspace}/tx_${ds}.txt`
do
	bam=${workspace}/${cond}_${tx}.bam
	echo ${bam}

	if [[ `samtools view -c ${bam}` > 0 ]]
	then
		python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
		--bam_file ${bam} \
		--mod_file ${mod} \
		--min_coverage 1 \
		--out_dir ${output} \
		--out_filename ${bam//.bam/.bed} \
		--num_jobs 36
	fi
done
