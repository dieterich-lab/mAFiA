library("scales")
library (ggplot2)


get_metagene_coord <- function(in_path) {
  dist <- read.delim(in_path, header =T)
  
  ### rescale ###
  utr5.SF <- median(dist$utr5_size, na.rm = T)/median(dist$cds_size, na.rm = T)
  utr3.SF <- median(dist$utr3_size, na.rm = T)/median(dist$cds_size, na.rm = T)
  utr5.dist <- dist[m6a.dist$rel_location < 1, ]
  cds.dist <- dist [dist$rel_location < 2 & dist$rel_location >= 1, ]
  utr3.dist <- dist[dist$rel_location >= 2, ]
  utr5.dist$rel_location <- rescale(utr5.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
  utr3.dist$rel_location <- rescale(utr3.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
  metagene.coord <- c(utr5.dist$rel_location, cds.dist$rel_location, utr3.dist$rel_location)
  
  return(metagene.coord)
}


res_dir <- "/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1/mouse_heart/metaPlotR"
ds <- "HFpEF"

cond <- "ctrl_merged"
ctrl.m6a.metagene.coord <- get_metagene_coord(file.path(res_dir, sprintf("%s_%s_m6a.dist.measures.txt", ds, cond)))
ctrl.psi.metagene.coord <- get_metagene_coord(file.path(res_dir, sprintf("%s_%s_psi.dist.measures.txt", ds, cond)))

cond <- "HFpEF_merged"
HFpEF.m6a.metagene.coord <- get_metagene_coord(file.path(res_dir, sprintf("%s_%s_m6a.dist.measures.txt", ds, cond)))
HFpEF.psi.metagene.coord <- get_metagene_coord(file.path(res_dir, sprintf("%s_%s_psi.dist.measures.txt", ds, cond)))

### m6A vs psi ###
# ctrl
ctrl.metagene.cord <- c(ctrl.m6a.metagene.coord, ctrl.psi.metagene.coord)
mod <- c(rep("m6A", length(ctrl.m6a.metagene.coord)), 
         rep("psi", length(ctrl.psi.metagene.coord))) 
df <- data.frame(ctrl.metagene.cord, mod)

ggplot(df) + geom_density(aes(x = ctrl.metagene.cord, colour = mod)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")

# HFpEF
HFpEF.metagene.cord <- c(HFpEF.m6a.metagene.coord, HFpEF.psi.metagene.coord)
mod <- c(rep("m6A", length(HFpEF.m6a.metagene.coord)), 
         rep("psi", length(HFpEF.psi.metagene.coord))) 
df <- data.frame(HFpEF.metagene.cord, mod)

ggplot(df) + geom_density(aes(x = HFpEF.metagene.cord, colour = mod)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")

### ctrl vs HFpEF ###
# m6A
m6a.metagene.cord <- c(ctrl.m6a.metagene.coord, HFpEF.m6a.metagene.coord)
mod <- c(rep("ctrl", length(ctrl.m6a.metagene.coord)), 
         rep("HFpEF", length(HFpEF.m6a.metagene.coord))) 
df <- data.frame(m6a.metagene.cord, mod)

ggplot(df) + geom_density(aes(x = m6a.metagene.cord, colour = mod)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")

# psi
psi.metagene.cord <- c(ctrl.psi.metagene.coord, HFpEF.psi.metagene.coord)
mod <- c(rep("ctrl", length(ctrl.psi.metagene.coord)), 
         rep("HFpEF", length(HFpEF.psi.metagene.coord))) 
df <- data.frame(psi.metagene.cord, mod)

ggplot(df) + geom_density(aes(x = psi.metagene.cord, colour = mod)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")