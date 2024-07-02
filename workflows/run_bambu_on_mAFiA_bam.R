#!/usr/bin/env Rscript

library(bambu)
# args = commandArgs(trailingOnly=TRUE)
# workspace <- args[1]
# out.dir <- args[2]
workspace <- "/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/TAC/SHAM_day56"
out.dir <- file.path(workspace, "bambu")

bam <- file.path(workspace, "Fhl1.mAFiA.reads.bam")
dir.create(out.dir)
# fa <- "/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"
# gtf <- "/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf"
fa <- '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/GRCm38_102_X.fa'
gtf <- '/home/adrian/Data/genomes/mus_musculus/GRCm38_102/Fhl1.GRCm38.102.gtf'
annot <- prepareAnnotations(gtf)

se <- bambu(reads = bam, annotations = annot, genome = fa, trackReads = TRUE, discovery=TRUE, NDR=1, ncore=4)
tx.map <- metadata(se)$readToTranscriptMaps[[1]]

write.table(apply(data.frame(tx.map), 2, as.character), file.path(out.dir, "Fhl1.tx_map.tsv"), sep="\t", row.names = FALSE)
write.table(rowData(se)[["TXNAME"]], file.path(out.dir, "Fhl1.tx_id.tsv"), sep="\t", col.names = FALSE, row.names = FALSE)

# mask <- !sapply(metadata(se)$readToTranscriptMaps[[1]][["equalMatches"]], is.null)
# mask.equalMatches <- sapply(tx.map[["equalMatches"]], length)==1
# readIds <- tx.map[["readId"]][mask]
# bambu_indices.equalMatches <- unlist(tx.map[["equalMatches"]][mask])
# txNames.equalMatches <- rowData(se)[["TXNAME"]][bambu_indices.equalMatches]
# geneIds <- rowData(se)[["GENEID"]][bambu_indices.equalMatches]
# df_out <- data.frame(readIds=readIds, txNames=txNames, geneIds=geneIds)
# df_out <- df_out[order(df_out$txNames),]
# 
# write.table(df_out, file.path(out.dir, "readIds_txNames.tsv"), sep="\t", row.names=FALSE)
# 
# for (this_tx in unique(df_out$txNames)) {
#     sub_df <- df_out[df_out$txNames==this_tx, ]
#     write(sub_df[["readIds"]], file.path(out.dir, sprintf("%s_%s.ids", sub_df[["geneIds"]][1], this_tx)))
# }
