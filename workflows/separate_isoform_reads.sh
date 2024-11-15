#!/bin/bash

ref_all=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/linearized.GRCm38.cdna.all.fa
annot=/biodb/genomes/mus_musculus/GRCm38_102/GRCm38.102.gtf

backbone=${HOME}/git/mAFiA_dev/models/RODAN_HEK293_IVT.torch
classifiers=${HOME}/git/mAFiA_dev/models/psi-co-mAFiA

#ds=TAC
#conditions="SHAM TAC"
#bambu="/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/DTU/DTU_TAC_Dewenter_ENS102.tsv"

ds=HFpEF
conditions="ctrl HFpEF"
bambu="/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/DTU/DTU_HFpEF_Laura_ENS102.tsv"

gene=Mfn2

gene_gtf=`grep $gene /biodb/genomes/mus_musculus/GRCm38_102/genes.gtf`
chrom=`echo $gene_gtf | cut -d' ' -f1`
chromStart=`echo $gene_gtf | cut -d' ' -f4`
chromEnd=`echo $gene_gtf | cut -d' ' -f5`
filter="${chrom}:${chromStart}-${chromEnd}"
transcripts=`grep ${gene} $bambu | cut -f5 | tr '\n' ' '`

workspace=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/${ds}/isoforms/${gene}
ref=${workspace}/isoforms.fa
mkdir ${workspace}
cd ${workspace}

source ${HOME}/git/mAFiA/mafia-venv/bin/activate
mod=/prj/TRR319_RMaP_BaseCalling/Adrian/site_annotations/mus_musculus/GRCm38_102/m6A.psi.GRCm38_102.chr${chrom}.bed

for tx in ${transcripts}
do
	grep ${tx} -A+1 ${ref_all} >> ${ref}
	grep ${tx} ${annot} > ${tx}.gtf
done

for cond in ${conditions}
do
	in_bam=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/${ds}/${cond}_merged.sorted.bam

	samtools view ${in_bam} ${filter} -o ${cond}.bam
	samtools fasta ${cond}.bam > ${cond}.fasta

	minimap2 -ax map-ont ${ref} ${cond}.fasta > ${cond}_isoforms.sam

	for tx in ${transcripts}
	do
		echo $tx
		samtools view -q1 -F260 ${cond}_isoforms.sam | grep ${tx} | cut -f1 > read_ids_${tx}_${cond}.txt
		samtools view -N read_ids_${tx}_${cond}.txt ${cond}.bam -o ${tx}_${cond}.bam
		samtools index ${tx}_${cond}.bam
		samtools view -c ${tx}_${cond}.bam
	done
	
	fast5_dir=/prj/TRR319_RMaP_BaseCalling/Adrian/mouse_heart/${ds}/fast5_${cond}
	for tx in ${transcripts}
	do
		srun -p gpu -c 4 --mem 90GB --gres gpu:turing:1 \
		python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_process_reads_parallel.py \
		--bam_file ${workspace}/${tx}_${cond}.bam \
		--fast5_dir ${fast5_dir} \
		--mod_file ${mod} \
		--backbone_model_path ${backbone} \
		--classifier_model_dir ${classifiers} \
		--num_jobs 4 \
		--batchsize 256 \
		--out_dir ${workspace} \
		--out_filename ${tx}_${cond}.mAFiA.bam
		
		srun -p general -c 40 --mem 120GB \
		python3 ${HOME}/git/mAFiA_dev/mAFiA/mAFiA_pileup.py \
		--bam_file ${tx}_${cond}.mAFiA.bam \
		--mod_file ${mod} \
		--min_coverage 1 \
		--out_dir ${workspace} \
		--out_filename ${tx}_${cond}.mAFiA.bed \
		--num_jobs 36
		
		samtools view -h ${tx}_${cond}.mAFiA.bam | sed -e 's/N+21891/N+a/g' | sed -e 's/N-21891/N-a/g' | samtools view -o relabelled.${tx}_${cond}.mAFiA.bam
		samtools index relabelled.${tx}_${cond}.mAFiA.bam
	done
done
