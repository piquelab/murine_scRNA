
######################################################################
# 1. preprocessing and clustering steps
######################################################################
library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(future)
library(readr)
library(DoubletFinder)
library(SoupX)
library(diem)
library(tidyverse)


outFolder="./4_harmony_cellClass_soupx_doubletfinder_chrM/"
system(paste0("mkdir -p ", outFolder))

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


samples<-read_tsv("Bacteria_inducedPTLproject-Details.txt")


#Load in the required files (barcodes.tsv, genes.tsv, and matrix.mtx), from the raw data files (raw_feature_bc_matrix).

forders <- list.files("../counts_BIPTL_2021-05-01",pattern="^BIPTL_*")


sc_list<-sapply(forders, function(x){
  #x<-forders[i]
  print(x)
  #mouse.data <- Read10X(data.dir = paste0("../counts_BIPTL_2021-05-01/",x,"/outs/filtered_feature_bc_matrix"))
  mouse.data<-load10X(paste0("../counts_BIPTL_2021-05-01/",x,"/outs/"))

  
 
  ##############################################################################
  # SoupX automatic way
  ##############################################################################
  
  sc = autoEstCont(mouse.data)
  out = adjustCounts(sc)
  mouse.data<-out
  
  #################################################################################
  # creating seurat object
  #################################################################################
  
  sc4 <- CreateSeuratObject(counts = mouse.data, project = "Mouse",min.cells = 3, min.features=200)
  sc4@meta.data$Library<-rep(x,nrow(sc4@meta.data))
  sc4@meta.data$Location<-rep(samples %>% filter(Sample==x) %>%select(Tissue)%>% unlist, nrow(sc4@meta.data))
  sc4@meta.data$Condition<-rep(samples %>% filter(Sample==x) %>%select(Condition)%>% unlist, nrow(sc4@meta.data))
  sc4@meta.data$SampleindexID<-rep(samples %>% filter(Sample==x) %>%select(SampleindexID)%>% unlist, nrow(sc4@meta.data))
  
  sc4 <- NormalizeData(sc4, verbose=TRUE)

  sc4 <- FindVariableFeatures(sc4, selection.method = "vst", nfeatures = 3000)

  sc4 <- ScaleData(sc4, verbose = TRUE)
  sc4 <- RunPCA(sc4)
  sc4 <- RunUMAP(sc4, dims = 1:10)

  #################################################################################
  # doublet removal
  #################################################################################
  
  nExp_poi <- round(0.075*nrow(sc4@meta.data))  ## Assuming 7.5% doublet formation rate
 
  #as maxima in BCmvn distributions
  sc_dbf <- doubletFinder_v3(sc4, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)

  DF.name = colnames(sc_dbf@meta.data)[grepl("DF.classification", colnames(sc_dbf@meta.data))]
  sc_dbf = sc_dbf[, sc_dbf@meta.data[, DF.name] == "Singlet"]
  sc_dbf

  } )

sc <- merge(sc_list[[1]],sc_list[-1], project="mouse-pilot")

sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 ) # change to: #>200  #<20000 , #doubletfinder #soupx  #diem

opfn <- paste0(outFolder,"seuratObj-newfilter-merge.",Sys.Date(),".rds") 
write_rds(sc, opfn)

# before MT filtering
sc<-readRDS("4_harmony_cellClass_soupx_doubletfinder/seuratObj-newfilter-merge.2021-06-21.rds")



#################################################################
# Quality control
#################################################################


####################################################
# percentage of reads mapping to mitochondrial
####################################################





cmd <- paste0("cat ","/wsu/home/groups/piquelab/data/refGenome10x/refdata-gex-mm10-2020-A/genes/genes.gtf",
              " | awk '$3~/gene/'",
              " | sed 's/gene_id //;s/;.*gene_type/\t/;s/; gene_name /\t/;s/; level.*//'")
cat(cmd,"\n")



aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5)) %>%
  select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11) 


anno <- tibble(gene_name=rownames(sc),rs=rowSums(sc@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>% filter(!is.na(Chr))

table(is.na(anno$Chr))

table(anno$Chr)

table(anno$Type)

head(anno)

head(aux)


sc <- sc[anno$gene_name,]


sc[["percent.mt"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrM",]$gene_name)

sc[["percent.mt"]] %>% summary()



# dens <- density(sc[["percent.mt"]]$percent.mt)
# # plot density
# plot(dens, frame = FALSE, col = "steelblue", 
#      main = "mitochondrial genes ") 


# ## Filter sc for things matching genotype. chrM or number of RNAs.  
# 
# sc[["percent.Y"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrY",]$gene_name)
# sc[["percent.Y"]] %>% summary()
# 
# ##anno[anno$gene_name=="XIST",]$kbid 
# sc[["percent.XIST"]] <- PercentageFeatureSet(sc,features=anno[anno$gene_name=="XIST",]$kbid )
# sc[["percent.XIST"]] %>% summary()

# table(sc[["percent.XIST"]]>0.01,md$SNG.BEST.GUESS)


## Check 0 or 1-based coordinates. 




#sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
#VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))


sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 10)

opfn <- paste0(outFolder,"seuratObj-after-mt-filtering.",Sys.Date(),".rds") 
write_rds(sc, opfn)

sc<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/seuratObj-after-mt-filtering.2021-06-28.rds")

############################################################

## Clustering
############################################################

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


sc <- NormalizeData(sc, verbose=TRUE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE)

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

#before clustering
opfn <- paste0(outFolder,"seuratObj-before-clustering.",Sys.Date(),".rds") 
write_rds(sc, opfn)


sc <- FindClusters(sc, verbose = TRUE,resolution=0.5)
opfn <- paste0(outFolder,"sc.NormByLibrary.cellclassify_newfilter-res0.2.",Sys.Date(),".rds") 
write_rds(sc, opfn)


sc<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.4.2021-06-28.rds")
fname=paste0("5_harmony_cellClass_soupx-doubletfinder_chrM_plots/UMAP_Harmony-res0.4.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) #+ NoLegend()
dev.off()
