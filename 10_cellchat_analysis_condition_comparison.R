library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)

#######################################################################
# Location- specific CellChat analysis (preterm Ecoli and control)
#######################################################################

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster

clust2Name<-cell.type.annotation$Cell.type
clust2Name<-paste0(c(0:30),"_",clust2Name)
names(clust2Name)<-c(0:30)

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                                    "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                                    "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
names(cluster.Colors)<-clust2Name 


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)



outFolder="./10_CellChat_analysis_default_after_filter_200/"
system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

# CellChat requires two user inputs: 
# one is the gene expression data of cells, 
# and the other is either user assigned cell labels (i.e., label-based mode) 
# or a low-dimensional representation of the single-cell data

####################################
# Load data
####################################

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)


m2 = read_tsv("5_harmonySubTypesDGE/Epithelial/ClusterDEG.tsv")
m2<-m2 %>% filter(p_val_adj<0.1)
subtypes<-unique(m2$Celltype)


# Here we load a scRNA-seq data matrix and its associated cell meta data

sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")
sc2<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/seuratObj-after-mt-filtering.2021-06-28.rds")
sc2 <- NormalizeData(sc2, verbose=TRUE)
#data.input<-sc2@assays$RNA@data
#meta<-sc@meta.data
#data.input<-data.input[,rownames(meta)]

# meta$labels<-clust2Names[meta$seurat_clusters]
# data.input<-sc2@assays$RNA@data
# meta<-sc@meta.data
# data.input<-data.input[,rownames(meta)]
# meta$labels<-clust2Names[meta$seurat_clusters]


locations<-unique(sc2$Location)
conditions<-unique(sc2$Condition)
locations<-c("Cervix"   ,  "Decidua" )
mclapply(locations, function(xlocation)
  {
  
  for (xcondition in conditions)
  {
    sc2_filter<-subset(sc2,Location==xlocation & Condition==xcondition) 
    
    
    metadata<-sc@meta.data
    mapping<-metadata$seurat_clusters
    names(mapping)<-rownames(metadata)
    
    
    sc2_filter$barcode<-colnames(sc2_filter)
    sc2_filter$seurat_clusters<-mapping[sc2_filter$barcode]
    sc2_filter$celltype<-clust2Name[sc2_filter$seurat_clusters]
    cluster_filter<-names(which(table(sc2_filter$seurat_clusters)>200))
    
    

    metadata<-metadata %>% filter (Location==xlocation & Condition==xcondition) 
    metadata<-metadata %>% filter( seurat_clusters %in% cluster_filter )
    sc2_filter<-subset(sc2_filter,barcode %in% rownames(metadata) & seurat_clusters %in% cluster_filter) 
    
    metadata$labels<-clust2Names[metadata$seurat_clusters]
    data<-sc2_filter@assays$RNA@data
    
    cellchat <- createCellChat(object = data, meta = metadata, group.by = "labels")
    cellchat <- addMeta(cellchat, meta = metadata)
    cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    
    
    ####################################
    # Set the ligand-receptor interaction database
    ####################################
    
    CellChatDB <- CellChatDB.mouse
    
    # pdf(paste0(outFolder,"showDatabaseCategory.pdf"))
    # showDatabaseCategory(CellChatDB)
    # dev.off()
    
    # Show the structure of the database
    dplyr::glimpse(CellChatDB$interaction)
    
    
    # use a subset of CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    
    # simply use the default CellChatDB
    
    CellChatDB.use <- CellChatDB 
    
    
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    
    ########################################################################
    # Preprocessing the expression data for cell-cell communication analysis
    ########################################################################
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    # project gene expression data onto PPI network (optional)
    cellchat <- projectData(cellchat, PPI.human)
    
    
    # Part II: Inference of cell-cell communication network
    
    # Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    
    
    ########################################################################
    # Infer the cell-cell communication at a signaling pathway level
    ########################################################################
    
    # NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
    cellchat <- computeCommunProbPathway(cellchat)
    
    ########################################################################
    # Calculate the aggregated cell-cell communication network
    ########################################################################
    
    cellchat <- aggregateNet(cellchat)
    
    write_rds(cellchat,file=paste0(outFolder,"cellchat_",xlocation,"_",xcondition,"_",Sys.Date(),".rds"))
    
  }
  
},mc.cores=6)









