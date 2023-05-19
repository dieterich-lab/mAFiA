library(ensembldb)
# library(EnsDb.Hsapiens.v102)
# edbx <- filter(EnsDb.Hsapiens.v102, filter = ~ seq_name == "X")

library(AnnotationHub)
query(ah, c("EnsDb", "v102", "Homo Sapiens"))
ensdb <- ah[["AH89180"]]

in_csv <- "/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Christoph/m6anet/workflow_tx/inference/75_WT_25_IVT/data.site_proba.csv"
in_df <- read.csv(in_csv)

ids <- in_df[['transcript_id']]
starts <- in_df[['transcript_position']]
widths <- rep(1, length(ids))

for (i in (1:length(ids))) {
  ids[[i]] <- strsplit(ids[[i]], '\\.')[[1]][1]
}

rng_tx <- IRanges(start = starts[1:10], width = widths[1:10], names = ids[1:10])
rng_gnm <- transcriptToGenome(rng_tx, ensdb)