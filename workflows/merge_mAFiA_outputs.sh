#!/usr/bin/env bash
shopt -s globstar

cd ${WORKSPACE}

########################################################################################################################
### filter and merge bam ###############################################################################################
########################################################################################################################

echo "Merging mAFiA outputs at ${WORKSPACE}..."

### filter ###
#for chr in {1..19} X
#do
#  mv chr${chr}/mAFiA.reads.bam chr${chr}/_mAFiA.reads.bam
#  samtools index chr${chr}/_mAFiA.reads.bam
#  samtools view -h chr${chr}/_mAFiA.reads.bam ${chr} | samtools sort - > chr${chr}/mAFiA.reads.bam
#done

## merge bam ###
samtools merge -o chrALL.mAFiA.reads.bam chr*/mAFiA.reads.bam
samtools index chrALL.mAFiA.reads.bam

### merge bed ###
cp chr1/mAFiA.sites.bed chrALL.mAFiA.sites.bed
for chr in {2..22} X MT
do
  tail -n+2 chr${chr}/mAFiA.sites.bed >> chrALL.mAFiA.sites.bed
done
