library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)
library(NMF)
library(ggalluvial)


#######################################################################
# Location- specific CellChat analysis -spreadsheet details
#######################################################################



#######################################################################

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
#res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
res <- res %>% filter(!is.na(pvalue) & padj<0.1)

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster


cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
names(cluster.Colors)<-clust2Names #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



tobechanged<-c("5_Epithelial-1 (Basal)","28_Epithelial-9 (Secretory)","7_Epithelial-2 (Squamous)","8_Epithelial-3 (Squamous)","10_Epithelial-4 (Glandular)","13_Epithelial-5 (Luminal)" , "14_Epithelial-6 (Secretory)" ,"23_Epithelial-8 (Enterocyte)",
               "11_Epithelial-10 (Proliferative)" ,"20_Epithelial-7 (Glandular)")






#color_tobechanged<-cluster.Colors[clust2Names[which(clust2Names %in% tobechanged)]]

outFolder="./10_CellChat_locations/"



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
names(cluster.Colors)<-clust2Name #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")




locations<-c("Cervix" ,    "Decidua",    "Myometrium")
#locations<-"Myometrium"
conditions<-c("Control" ,"E. coli")


sapply(locations, function(xlocation){
  pathways_genes_celltypes_db_list<-lapply (locations ,function(xlocation){
    
    subFolder<-xlocation
    
    system(paste0("mkdir -p ", outFolder,subFolder,"/"))
    
    
    load(paste0("./10_CellChat_analysis_default/","cellchat_",xlocation,"_","2021-12-07",".RData"))
    pathways<-cellchat@netP$pathways
    groupSize <- as.numeric(table(cellchat@idents))
    
    
    
    # 
    system(paste0("mkdir -p ", outFolder,subFolder,"/Pathways_circleplots_top25"))

    sapply(pathways,function(pathways.show){
      pdf(paste0(outFolder,subFolder,"/Pathways_circleplots_top25/",pathways.show,".pdf"),width=20,height=20)
      par(mfrow=c(1,1))
      #arrow(length = unit(.02, "inches"),type = "closed",angle = 40)
      netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=2,top=0.25,arrow.width = 15)
      dev.off()
    })

    
  })
})


pathways_genes_celltypes_db_final<-lapply(locations, function(xlocation){
  
subFolder<-xlocation
system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  #pathways_genes_celltypes_db_list<-lapply (conditions ,function(xcondition){
    
pathways_genes_celltypes_db_list<-list()

    
load(paste0("./10_CellChat_analysis_default/","cellchat_",xlocation,"_","2021-12-07",".RData"))
    pathways<-cellchat@netP$pathways
    groupSize <- as.numeric(table(cellchat@idents))
    
      system(paste0("mkdir -p ", outFolder,subFolder,"/"))
    

    pathways_genes_celltypes<-lapply(pathways,function(x){
      
      pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
      LR_genes<-unique(unlist(str_split(pairLR$interaction_name,"_")))
      
      LR_genes<-tolower(LR_genes)
      
      
      celltypes<-res %>% filter(tolower(kbid) %in% LR_genes)%>% dplyr::select(Cell_type)%>% unlist %>% unique()
      
      LR_DEGs <-res %>% filter (tolower(kbid) %in% LR_genes) %>% dplyr::select(kbid)%>% unlist %>% unique()
      LR_DEGs<-paste(LR_DEGs,collapse = ", ")
      
      LR_genes<-paste(LR_genes,collapse = ", ")
      
      celltype_len<-length(celltypes)
      celltypes<-paste(celltypes,collapse = ", ")
      if (celltypes=="") celltypes<-"NA"
      return(c(x,LR_genes,LR_DEGs,celltypes,celltype_len))
    })
    pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes)
    
    pathways_genes_celltypes_db<-as.data.frame(pathways_genes_celltypes_db)
    pathways_genes_celltypes_db$location<-xlocation
    colnames(pathways_genes_celltypes_db)<-c("Pathway","Genes_in_pathway","DEGs_in_pathway","Celltype","celltype_len","Location")
    
    write.csv(pathways_genes_celltypes_db,file=paste0(outFolder,subFolder,"/pathways_genes_celltypes_db.csv"))
    
  
  return(pathways_genes_celltypes_db) 
  
})

pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes_db_final)
pathways_genes_celltypes_db<-pathways_genes_celltypes_db[order(pathways_genes_celltypes_db$celltype_len,decreasing = TRUE),]
write.csv(pathways_genes_celltypes_db,file=paste0(outFolder,"pathways_genes_celltypes_db.csv"))
























                                                            
  

