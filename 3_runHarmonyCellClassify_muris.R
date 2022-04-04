
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
library(SingleR)

########################################################
# cell type classification using SingleR 
# reference: Tabula Muris Consortium. "Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris."Nature 562, no. 7727 (2018): 367-372.

########################################################

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

outFolder="./muris_newsubset/"
system(paste0("mkdir -p ", outFolder))

# data at /wsu/home/groups/prbgenomics/mouse_pilot/mouse_pilot_analysis/tabula-muris/tabula-muris/00_data_ingest


# FACS 
# uploaded in grid: "tabula-muris/tabula-muris/00_data_ingest/00_facs_raw_data"
# https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells/5715040

# droplet 
# uploaded in grid: "tabula-muris/tabula-muris/00_data_ingest/01_droplet_raw_data" 
# https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_microfluidic_emulsion/5715025


#######################################################
# seurat R objects

# uploaded in grid: 
#"tabula-muris/tabula-muris/00_data_ingest/04_tissue_robj_generated"

seurat_files<-list.files("tabula-muris/tabula-muris/00_data_ingest/04_tissue_robj_generated/")
rm_files<-c("droplet_Marrow_seurat_tiss.Robj" ,"droplet_Tongue_seurat_tiss.Robj","droplet_Tongue_seurat_tiss.Robj","facs_Brain_Non-Myeloid_seurat_tiss.Robj",
            "facs_Skin_seurat_tiss.Robj" ,"facs_Tongue_seurat_tiss.Robj" ,"facs_Diaphragm_seurat_tiss.Robj" ,"facs_Marrow_seurat_tiss.Robj","facs_Brain_Myeloid_seurat_tiss.Robj" ,"facs_Brain_Non-Myeloid_seurat_tiss.Robj"    )

readfolder<-"./tabula-muris/tabula-muris/00_data_ingest/04_tissue_robj_generated/"
seurat_files<-seurat_files[!seurat_files %in% rm_files]


sc_list<-lapply(seurat_files, function(x)
  {
  print(x)
  load(paste0(readfolder,x))
  sc<-tiss
  sc_updated<-UpdateSeuratObject(object = sc)
  sc_updated
  })


sc1 <- merge(sc_list[[1]],sc_list[-1], project="muris")
# sc1@meta.data$cell_ontology_class


sc2<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/seuratObj-after-mt-filtering.2021-06-28.rds")


sc <- merge(sc1,list(sc2))


sc <- NormalizeData(sc, verbose=TRUE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE)

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

lib<-names(sc$Library)[which(is.na(sc$Library))]
sc$Library[which( colnames(sc) %in% lib)]<- sc1$channel [ colnames(sc)[which(colnames(sc) %in% lib)]]

sc <- RunHarmony(sc,c("Library"),reduction="pca")



sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)

###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE, resolution=0.5)





he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

query.he <- he[,is.na(sc@meta.data$cell_ontology_class)]
ref.he <- he[,!is.na(sc@meta.data$cell_ontology_class)]
ref.labels <- sc@meta.data$cell_ontology_class[!is.na(sc@meta.data$cell_ontology_class)]


pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))

fname=paste0(outFolder,"sc.NormByLibrary.ref.Harmony.singler.Muris.res0.5.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
  rownames_to_column("BARCODES") %>%
  left_join(sc2@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLibrary.refMuris.Harmony.singler.res0.5.csv")
write_csv(md,fname)

## save object.
fname=paste0(outFolder,"sc.NormByLibFullIntegrated.refMuris.Harmony.res0.5.rds")
write_rds(sc,fname)


sc <- read_rds("muris/sc.NormByLibrary.ref.Harmony.singler.Muris.res0.5.rds")
