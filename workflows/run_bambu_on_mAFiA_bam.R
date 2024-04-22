#!/usr/bin/env Rscript

library(bambu)
args = commandArgs(trailingOnly=TRUE)
workspace <- args[1]
out.dir <- args[2]
# workspace <- "/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA/HEK293/100_WT_0_IVT/chrX"
# out.dir <- file.path(workspace, "bambu")

bam <- file.path(workspace, "mAFiA.reads.bam")
dir.create(out.dir)
gtf <- "/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf"
annot <- prepareAnnotations(gtf)
fa <- "/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"

se <- bambu(reads = bam, annotations = annot, genome = fa, trackReads = TRUE, discovery=FALSE, ncore=36)

# mask <- !sapply(metadata(se)$readToTranscriptMaps[[1]][["equalMatches"]], is.null)
mask <- sapply(metadata(se)$readToTranscriptMaps[[1]][["equalMatches"]], length)==1
readIds <- metadata(se)$readToTranscriptMaps[[1]][["readId"]][mask]
bambu_indices <- unlist(metadata(se)$readToTranscriptMaps[[1]][["equalMatches"]][mask])
txNames <- rowData(se)[["TXNAME"]][bambu_indices]
geneIds <- rowData(se)[["GENEID"]][bambu_indices]

df_out <- data.frame(readIds=readIds, txNames=txNames, geneIds=geneIds)
df_out <- df_out[order(df_out$txNames),]
write.table(df_out, file.path(out.dir, "readIds_txNames.tsv"), sep="\t", row.names=FALSE)

min_coverage <- 20
for (this_tx in unique(df_out$txNames)) {
    sub_df <- df_out[df_out$txNames==this_tx, ]
    if (dim(sub_df)[1]>=min_coverage) {
        write(sub_df[["readIds"]], file.path(out.dir, sprintf("%s_%s.ids", sub_df[["geneIds"]][1], this_tx)))
    }
}
