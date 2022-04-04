library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
library(SingleR)

########################################################
# cell type classification using SingleR 
# reference: Cross-tissue organization of the fibroblast lineage
########################################################

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

outFolder="./4_harmony_cellClass_fibroblast/"
system(paste0("mkdir -p ", outFolder))

#######################################################
# seurat R objects

# uploaded in grid: 
Mouse_SS_Fibro<-read_rds("fibroblast_cross-tissue/Mouse_SS_Fibro.RDS")
sc1<-Mouse_SS_Fibro

sc2<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/seuratObj-after-mt-filtering.2021-06-28.rds")


sc <- merge(sc1,list(sc2))


sc <- NormalizeData(sc, verbose=TRUE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE)

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

loc<-names(sc$Location)[which(is.na(sc$Location))]

tissues<-sc1@meta.data$Tissue
names(tissues)<-rownames(sc1@meta.data)
names(tissues)<-paste0(names(tissues),"_1")
sc$Location[which( colnames(sc) %in% loc)]<- tissues[loc ] #sc1$Tissue [ colnames(sc)[which(colnames(sc) %in% loc)]]

sc <- RunHarmony(sc,c("Location"),reduction="pca")



sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)

###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE, resolution=0.5)





he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

query.he <- he[,is.na(sc@meta.data$ClustName)]
ref.he <- he[,!is.na(sc@meta.data$ClustName)]
ref.labels <- sc@meta.data$ClustName[!is.na(sc@meta.data$ClustName)]


pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))

fname=paste0(outFolder,"sc.NormByLibrary.ref.Harmony.singler.refFibroblast.res0.5.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
  rownames_to_column("BARCODES") %>%
  left_join(sc2@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLibrary.refFibroblast.Harmony.singler.res0.5.csv")
write_csv(md,fname)

## save object.
fname=paste0("4_harmony_cellClass_fibroblast/sc.NormByLibFullIntegrated.refFibroblast.Harmony.res0.5.rds")
write_rds(sc,fname)


pred.labels<-read_rds(paste0("4_harmony_cellClass_fibroblast/sc.NormByLibrary.ref.Harmony.singler.refFibroblast.res0.5.rds"))



