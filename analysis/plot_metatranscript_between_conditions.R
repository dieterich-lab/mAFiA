library("scales")
library (ggplot2)


get_metagene_coord <- function(in_path, rescale=FALSE) {
  dist <- read.delim(in_path, header =T)
  
  if (rescale == TRUE) {
    utr5.SF <- median(dist$utr5_size, na.rm = T)/median(dist$cds_size, na.rm = T)
    utr3.SF <- median(dist$utr3_size, na.rm = T)/median(dist$cds_size, na.rm = T)
    utr5.dist <- dist[dist$rel_location < 1, ]
    cds.dist <- dist [dist$rel_location < 2 & dist$rel_location >= 1, ]
    utr3.dist <- dist[dist$rel_location >= 2, ]
    utr5.dist$rel_location <- rescale(utr5.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
    utr3.dist$rel_location <- rescale(utr3.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
    metagene.coord <- c(utr5.dist$rel_location, cds.dist$rel_location, utr3.dist$rel_location)
    return(metagene.coord)
  }
  else {
    return(dist$rel_location)
  }
}

compare_mods <- function(in.metagene.coord, in_cond) {
  transcript.region <- c(in.metagene.coord[[in_cond]][["m6A"]], 
                           in.metagene.coord[[in_cond]][["psi"]])
  mods <- c(rep("m6A", length(in.metagene.coord[[in_cond]][["m6A"]])), 
           rep("psi", length(in.metagene.coord[[in_cond]][["psi"]]))) 
  df <- data.frame(transcript.region, mods)
  
  ggplot(df) + 
    geom_density(aes(x = transcript.region, colour = mods)) + 
    xlim(0, 3) +
    theme_bw() +
    geom_vline(xintercept = 1:2, col = "grey")
}


compare_conditions <- function(in.metagene.coord, in_mod) {
  transcript.region <- c(in.metagene.coord[["ctrl_merged"]][[in_mod]], 
                         in.metagene.coord[["HFpEF_merged"]][[in_mod]])
  conditions <- c(rep("ctrl_merged", length(in.metagene.coord[["ctrl_merged"]][[in_mod]])),
                  rep("HFpEF_merged", length(in.metagene.coord[["HFpEF_merged"]][[in_mod]])))
  df <- data.frame(transcript.region, conditions)
  
  ggplot(df) + 
    geom_histogram(aes(x = transcript.region, colour = conditions)) +
    xlim(0, 3) +
    theme_bw() +
    geom_vline(xintercept = 1:2, col = "grey")
}


res_dir <- "/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR"

ds <- "HFpEF"
all_conditions <- c("ctrl_merged", "HFpEF_merged")
all_mods <- c("m6A", "psi")
# thresholds <- c("modRatio50.0", "modRatio0.0")
# thresholds <- c("modRatio0.0")

### get gene body coverage ###
coverage.file <- "/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/HFpEF/ctrl_merged.geneBodyCoverage.txt"
coverage_df <- t(read.table(coverage.file))
colnames(coverage_df) <- coverage_df[1,]
coverage_df <- coverage_df[-1,]
row.names(coverage_df) <- NULL
coverage_df <- as.data.frame(apply(coverage_df, 2, as.numeric))
coverage_profile <- coverage_df$chrALL.mAFiA.reads / sum(coverage_df$chrALL.mAFiA.reads)

metagene.coord <- list()
bins <- seq(0, 3, 0.05)
for (this_cond in all_conditions) {
  for (this_mod in all_mods) {
    metagene.coord[[this_cond]][[this_mod]] <- get_metagene_coord(file.path(res_dir, sprintf("%s_%s_%s_%s.dist.measures.txt", ds, this_cond, this_mod, "modRatio50.0")))
  }
}

### m6A vs psi ###
compare_mods(metagene.coord, "ctrl_merged")
compare_mods(metagene.coord, "HFpEF_merged")


### ctrl vs HFpEF ###
# m6A
compare_conditions(metagene.coord, "m6A")
compare_conditions(metagene.coord, "psi")
