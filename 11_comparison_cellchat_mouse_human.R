
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)

####################################################################################
# Comparison between mouse and human
####################################################################################

outFolder<-"./10_human_mouse_myometrium/"
system(paste0("mkdir -p ", outFolder))

######################################################
# mouse 
######################################################
pathways_mouse<-read_rds(paste0("./10_CellChat_comparison_conditions_after_filter_200_plots/pathways_genes_celltypes_db.rds"))
pathways_mouse<-pathways_mouse %>% select(Pathway,DEGs_mouse=DEGs_in_pathway,DEGs_len_mouse=LR_DEGs_len,Celltype_mouse=Celltype,Condition_mouse=Condition,Location_mouse=Location)
pathways_mouse<-pathways_mouse %>% filter (Location_mouse=="Myometrium" & Condition_mouse=="Ecoli")
pathways_mouse$DEGs_len_mouse<-as.numeric(pathways_mouse$DEGs_len_mouse)
pathways_mouse<- pathways_mouse %>% arrange(desc(DEGs_len_mouse))
pathways_mouse$Pathway[1:10]


# DEGs
cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster
clust2Name<-cell.type.annotation$Cell.type
clust2Name<-paste0(c(0:30),"_",clust2Name)
names(clust2Name)<-c(0:30)
# cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
#                   "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
#                   "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
# names(cluster.Colors)<-clust2Name #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
res <- read_tsv("/wsu/home/groups/prbgenomics/mouse_pilot/mouse_pilot_analysis/7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
res_mouse<-res
res_mouse<-res_mouse %>% filter(Location=="Myometrium")
res_mouse$Cell_type_main_mouse<-
res_mouse$Cell_type_main<-gsub("[^a-zA-Z]", "", res_mouse$Cell_type)
res_mouse$Cell_type_main[which(res_mouse$Cell_type_main %in% c("EpithelialGlandular","EpithelialSecretory","EpithelialSquamous","EpithelialProliferative","EpithelialLuminal"))]<-"Epithelial"
res_mouse$Cell_type_main[which(res_mouse$Cell_type_main %in% c("MacrophageProgenitor"))]<-"Macrophage"

res_mouse$Cell_type_main[which(res_mouse$Cell_type_main %in% c("StromalDecidua"))]<-"Stromal"




######################################################
### human
######################################################
pathways_human<-read_csv(paste0("/wsu/home/groups/prbgenomics/labor_myo/myometrium_analysis/10_CellChat_analysis_customized/pathways_DEgenes_celltypes_roles_db.csv"))
pathways_human<-pathways_human[,-1]
colnames(pathways_human)<-c("Pathway","DEGs_in_pathway","Celltype","sender","receiver")
pathways_human<-pathways_human %>% select (Pathway,DEGs_human=DEGs_in_pathway,Celltype_human=Celltype)
pathways_human$DEGs_len_human<-sapply(pathways_human$DEGs_human, function(x)return(as.numeric(length(unlist(str_split(x,", "))))))
pathways_human<- pathways_human %>% arrange(-DEGs_len_human)

pathways_human<- pathways_human %>% arrange(desc(DEGs_len_human))
pathways_human$Pathway[1:10]

res <- read_tsv("/wsu/home/groups/prbgenomics/labor_myo/myometrium_analysis/7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-08-30.tsv")
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res<-res %>% filter(padj<0.1)
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]
res_human<-res
res_human$Cell_type_main<-gsub("[^a-zA-Z]", "", res_human$Cell_type)
res_human$Cell_type_main[which(res_human$Cell_type_main=="Smoothmusclecells")]<-"SMC"
res_human$Cell_type_main[which(res_human$Cell_type_main=="CDTcell")]<-"Tcell"



selected_pathways<-intersect(pathways_mouse$Pathway[1:10],pathways_human$Pathway[1:10])


##############################
# pathways and DEG#
##############################


pathway_map<-pathways_human %>% inner_join(pathways_mouse)

pathway_map<-pathway_map %>% arrange(desc(DEGs_len_mouse,DEGs_len_human))


pathway_map$DEGs_overlap<-sapply(1:nrow(pathway_map),function(x){
  hm<-unlist(str_split(pathway_map$DEGs_human[x],", "))
  mm<-unlist(str_split(pathway_map$DEGs_mouse[x],", "))
  overlap<-intersect(toupper(hm),toupper(mm))
  overlap_len<-length(overlap)
  
  })

pathway_map<-pathway_map %>% arrange(desc(DEGs_overlap,DEGs_len_mouse,DEGs_len_human))

pathway_map<-pathway_map %>% select(Pathway,DEGs_overlap,DEGs_len_human,DEGs_len_mouse,DEGs_human,DEGs_mouse,Celltype_human,Celltype_mouse)

write.csv(pathway_map,file="./10_human_mouse_myometrium/pathway_map.csv")





####################################
#### comparison per pathway
####################################

cor_all<-sapply(selected_pathways,function(x){
  
  
  system(paste0("mkdir -p ", outFolder,"/",x))
  
  
  genes_human<-pathways_human %>% filter(Pathway ==x) %>% select(DEGs_human)%>% unlist()
  genes_human<-unlist(str_split(genes_human,", "))
  
  res_human_select<-res_human%>% filter (gene_name %in% genes_human  & padj<0.1)
  res_human_select<-res_human_select %>% select(gene_name,Cell_type_human=Cell_type,log2FoldChange_human=log2FoldChange,padj_human=padj,Cell_type_main)
  
  
  genes_mouse<-pathways_mouse %>% filter(Pathway ==x) %>% select(DEGs_mouse)%>% unlist()
  genes_mouse<-unlist(str_split(genes_mouse,", "))
  
  res_mouse_select<-res_mouse%>% filter (kbid %in% genes_mouse  & padj<0.1)
  res_mouse_select<-res_mouse_select %>% select(gene_name=kbid,Cell_type_mouse=Cell_type,log2FoldChange_mouse=log2FoldChange,padj_mouse=padj,Cell_type_main)
  res_mouse_select$gene_name<-toupper(res_mouse_select$gene_name)
  
  
  resJoin <- res_human_select %>% inner_join(res_mouse_select) 
  resJoin_all <- res_human_select %>% inner_join(res_mouse_select,by="gene_name") 
  
  fname<-paste0(outFolder,"/",x,"/")
  
  cor_result<-cor.test(resJoin$log2FoldChange_human,resJoin$log2FoldChange_mouse,method="spearman",na.rm=TRUE)
  print(x)
  print(cor.test(resJoin$log2FoldChange_human,resJoin$log2FoldChange_mouse,method="spearman",na.rm=TRUE))
  
  res_mouse_select<-res_mouse_select %>% group_by(Cell_type_mouse) %>% ungroup()
  res_human_select<-res_human_select %>% group_by(Cell_type_human) %>% ungroup()
  write.csv(res_mouse_select,file=paste0( fname, x,"_DEGs_mouse.csv"  ))
  write.csv(res_human_select,file=paste0( fname, x,"_DEGs_human.csv"  ))
  write.csv(resJoin,file=paste0( fname, x,"_map_human_mouse.csv"  ))
  return (c(cor_result$estimate,cor_result$p.value))
  
})


#######

pathways_human<-pathways_human %>% arrange(desc(DEGs_len_human))
pathways_mouse<-pathways_mouse %>% arrange(desc(DEGs_len_mouse))
overlap<-intersect(pathways_human$Pathway , pathways_mouse$Pathway)
pathways_human<-pathways_human %>% filter (Pathway %in% overlap)
pathways_mouse<-pathways_mouse %>% filter (Pathway %in% overlap)

colnames(pathways_human)[1]<-"Pathway_human"
colnames(pathways_mouse)[1]<-"Pathway_mouse"
ranked_pathways<-cbind(pathways_human,pathways_mouse)

ranked_pathways<- ranked_pathways %>% select(Pathway_human,DEGs_len_human,Pathway_mouse,DEGs_len_mouse,DEGs_human,DEGs_mouse,Celltype_human,Celltype_mouse)

write.csv(ranked_pathways,file="./10_human_mouse_myometrium/ranked_pathways.csv")

