if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("bambu")

library(bambu)

### test ###
# test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
# fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
# gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
# bambuAnnotations <- prepareAnnotations(gtf.file)

### HEK293 ###
# test.bam <- "/prj/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/100_WT_0_IVT/filtered_q50.bam"
# fa.file <- "/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38_96.fa"
# gtf.file <- "/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf"

test.bam <- "/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/100_WT_0_IVT/filtered_q50.bam"
fa.file <- "/home/adrian/Data/genomes/GRCh38_96.fa"
gtf.file <- "/home/adrian/Data/transcriptomes/GRCh38.96.gtf"
outpath <- "/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/HEK293/100_WT_0_IVT/bambu"

### map by bambu ###
bambuAnnotations <- prepareAnnotations(gtf.file)
# se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file)
# se.quantOnly <- bambu(reads = test.bam, annotations = gtf.file, genome = fa.file, discovery = FALSE)
se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, trackReads = TRUE)
print(metadata(se)$readToTranscriptMaps[[1]])

### write output ###
writeBambuOutput(se, path = outpath)

df <- apply(as.data.frame(metadata(se)$readToTranscriptMaps[[1]]), 2, as.character)
write.table(
  df, 
  file = file.path(outpath, "readToTransciptMaps.tsv"),
  sep = "\t",
  row.names = FALSE
  )