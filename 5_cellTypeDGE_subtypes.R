#############################################################
### sub-type marker identification
### 
#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
library(clusterProfiler)
##library(SingleR)


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################

outFolder="./5_harmonySubTypesDGE/"
system(paste0("mkdir -p ", outFolder))



#anno<-read_rds("3_MergeDemux_Output/anno.rds")

sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")


dim(sc)

table(sc$Library)

table(sc$Location) 

system(paste0("mkdir -p ", outFolder))



cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster

celltypes<-cell.type.annotation$Cell.type
names(celltypes)<-cell.type.annotation$Cluster
   
   
   # clust2Name<-cell.type.annotation$Cell.type
# clust2Name<-paste0(c(0:30),"_",clust2Name)
# names(clust2Name)<-c(0:30)

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
names(cluster.Colors)<-clust2Names #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")




Stromal<-names(clust2Names)[which(celltypes %in% c("Stromal-2 (Decidua)" ,"Stromal-1","Stromal-3"))]
SMC<-names(clust2Names)[which(celltypes %in% c("SMC-1","SMC-2"))]
Macrophage<-names(clust2Names)[which(celltypes %in% c("Macrophage-2 (Progenitor)","Macrophage"))]
Epithelial<-names(clust2Names)[which(celltypes %in% c("Epithelial-1 (Basal)","Epithelial-2 (Squamous)","Epithelial-3 (Squamous)","Epithelial-4 (Glandular)","Epithelial-10 (Proliferative)","Epithelial-5 (Luminal)","Epithelial-6 (Secretory)","Epithelial-7 (Glandular)","Epithelial-8 (Enterocyte)","Epithelial-9 (Secretory)"))]
Fibroblast<-names(clust2Names)[which(celltypes %in% c("Fibroblast-1","Fibroblast-2","Fibroblast-3"))]
NKcell<-names(clust2Names)[which(celltypes %in% c("NK-cell-1","NK-cell-2"))]


subtypes<-c("SMC","Macrophage","Stromal","Epithelial","Fibroblast","NKcell")

subtypes_list<-list(SMC,Macrophage,Stromal,Epithelial,Fibroblast,NKcell)

 sapply(1:length(subtypes),function(x){
   
   print(subtypes[x])
   sc1 <- subset(sc, subset = seurat_clusters %in% subtypes_list[[x]])
   markers <- FindAllMarkers(sc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
   markers$Celltype<-as.character(clust2Names[as.character(markers$cluster)])
   
   
   eg = bitr(markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
   names(eg)[1]="gene"
   head(eg)
   e2g <- eg$gene
   names(e2g) <- eg$ENTREZID
   
   
   m2 <- markers %>% left_join( eg) 
   m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)
   
   system(paste0("mkdir -p ", outFolder,subtypes[x],"/"))
   outFolder<-paste0(outFolder,subtypes[x],"/")
   
   fname=paste0(outFolder,"ClusterDEG.tsv");
   write_tsv(m2,fname)
    print(dim(m2))
   top100 <- m2 %>% top_n(n = 100, wt = avg_log2FC)
   
   fname=paste0(outFolder,"ClusterDEGtop100.tsv");
   write_tsv(top100,fname)
 })
 

 