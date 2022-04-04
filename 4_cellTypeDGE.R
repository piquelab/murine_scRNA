
#############################################################
### cluster marker identification
### 
#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)

sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")

outFolder="./5_harmonyClustersDGE-soupx-doubletfinder-newfilter_chrM_res0.5/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

dim(sc)

table(sc$Library)

table(sc$Location) 


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 60 * 1024 ^ 3)



markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

opfn <- paste0(outFolder,"markers.rds") 
write_rds(markers, opfn)
opfn <- paste0(outFolder,"markers.tsv") 
write_tsv(markers,opfn)

m2 <- markers %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)


top20 <- m2 %>% top_n(n = 20, wt = avg_log2FC)

fname=paste0(outFolder,"ClusterDEG.tsv");
write_tsv(m2,fname)
m2<-read_tsv(fname)
write.csv(m2,file=paste0(outFolder,"ClusterDEG.csv"))

fname=paste0(outFolder,"ClusterDEGtop20.tsv");
write_tsv(top20,fname)
top20<-read_tsv(fname)
write.csv(top20,file=paste0(outFolder,"ClusterDEGtop20.csv"))


