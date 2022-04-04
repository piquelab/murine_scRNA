library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)
library(NMF)
library(ggalluvial)


#######################################################################
# CellChat plots 
# alluvial plots, pathway-specific circle plots, information flow 
#######################################################################




outFolder="./10_CellChat_comparison_conditions_after_filter_200_plots/"



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

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)



#######################################################################
# alluvial plots
#######################################################################




# pathways to focus 
#informationflow_data<-read_rds(paste0("10_CellChat_comparison_conditions_after_filter_200_plots/informationflow_data.rds"))
selected_pathways_mat<-read.delim("selectedpathways.txt")

# locations
locations<-c("Cervix" ,    "Decidua",    "Myometrium")
conditions<-c("Control" ,"E. coli")


sapply(locations, function(xlocation){
  
subFolder<-xlocation
system(paste0("mkdir -p ", outFolder,subFolder,"/"))

lapply (conditions ,function(xcondition){

    cellchat<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_",xcondition,"_","2021-12-21",".rds"))
    pathways<-cellchat@netP$pathways
    groupSize <- as.numeric(table(cellchat@idents))
    
    
    if (xcondition=="E. coli") xcondition<-"Ecoli" 
    
  
    selected_pathways<-selected_pathways_mat[,xlocation]
    selected_pathways<-selected_pathways[selected_pathways!=""]
    system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))
    
  
  ############################################################################
  # alluvial plots
  ############################################################################
  
  
  
  nPatterns<-2 #
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

  cutoff<-0.5

  ############################################################################
  # outgoing data
  ############################################################################
  outgoing_signaling<-cellchat@netP$pattern$outgoing$pattern$signaling
  outgoing_signaling<-outgoing_signaling %>% filter(Contribution>cutoff)

  colnames(outgoing_signaling)[3]<-"Contribution_outgoing_signaling"
  colnames(outgoing_signaling)[1]<-paste0(colnames(outgoing_signaling)[1],"_outgoing")

  outgoing_cell<-cellchat@netP$pattern$outgoing$pattern$cell
  outgoing_cell<-outgoing_cell %>% filter(Contribution>cutoff)

  colnames(outgoing_cell)[3]<-"Contribution_outgoing_cell"
  colnames(outgoing_cell)[2]<-paste0(colnames(outgoing_cell)[2],"_outgoing")

  outgoing_data<-outgoing_signaling %>% inner_join(outgoing_cell)
  outgoing_data<-outgoing_data %>% filter(Signaling %in% selected_pathways)


  ############################################################################
  # incoming data
  ############################################################################
  incoming_signaling<-cellchat@netP$pattern$incoming$pattern$signaling
  incoming_signaling<-incoming_signaling %>% filter(Contribution>cutoff)
  #incoming_signaling$Contribution[incoming_signaling$Contribution < cutoff] <- 0

  colnames(incoming_signaling)[3]<-"Contribution_incoming_signaling"
  colnames(incoming_signaling)[1]<-paste0(colnames(incoming_signaling)[1],"_incoming")

  incoming_cell<-cellchat@netP$pattern$incoming$pattern$cell
  incoming_cell<-incoming_cell %>% filter(Contribution>cutoff)
  #incoming_cell$Contribution[incoming_cell$Contribution < cutoff] <- 0

  colnames(incoming_cell)[3]<-"Contribution_incoming_cell"
  colnames(incoming_cell)[2]<-paste0(colnames(incoming_cell)[2],"_incoming")

  incoming_data<-incoming_signaling %>% inner_join(incoming_cell)

  incoming_data<-incoming_data %>% filter(Signaling %in% selected_pathways)
  
  outgoing_data<-outgoing_data %>% filter(Signaling %in% selected_pathways)
  
  outgoing_data<-outgoing_data %>% mutate(CellGroup_outgoing=CellGroup)#,Signaling,Contribution-outgoing_cell,Contribution_outgoing_signaling)
  incoming_data<-incoming_data %>% mutate(CellGroup_incoming=CellGroup)#,Signaling,Contribution_incoming_signal,Contribution_incoming_cell)
  signalins_incoming<-unique(incoming_data$Signaling)
  outgoing_incoming<-unique(outgoing_data$Signaling)
  overlap_signalings<-intersect(signalins_incoming,outgoing_incoming)


  ############################################################################
  # outgoing alluvial plot
  ############################################################################

  
  # summary 
  majorcells<-rowSums(table(outgoing_data$CellGroup,outgoing_data$Signaling))
  majorcells<-majorcells[order(majorcells,decreasing = TRUE)]
  majorsignalings<-rowSums(table(outgoing_data$Signaling,outgoing_data$Contribution_outgoing_signaling))
  majorsignalings<-majorsignalings[order(majorsignalings,decreasing = TRUE)]
  
  major_cell_signaling<-list(majorcells,majorsignalings)
  names(major_cell_signaling)<-c("majorcells","majorsignalings")
  
  
  outgoing_data$Signaling<-as.character(outgoing_data$Signaling)
  outgoing_data<-outgoing_data[order(outgoing_data$Signaling,decreasing = FALSE),]
  outgoing_data$Signaling <- factor(outgoing_data$Signaling,levels=unique(outgoing_data$Signaling))

  #system(paste0("mkdir -p ", outFolder,subFolder,"/"))

   fname=paste0(outFolder,subFolder,"/",xcondition,"/",xcondition,"_outgoing_ggalluvial_selectedPathways_cuttoff_",Sys.Date(),"_.pdf");
  
  pdf(fname,width=20,height=10)
  ggplot(data = outgoing_data,
         aes(axis1 = CellGroup_outgoing, axis2 = Signaling, y = Contribution_outgoing_signaling,fill=CellGroup_outgoing)) +
    geom_alluvium(aes(fill=CellGroup_outgoing)) +#aes(fill = Signaling)
    geom_stratum() +
    geom_flow()+
    scale_fill_manual("Cell type",values=cluster.Colors) +
    geom_text(stat = "stratum",aes(label = after_stat(stratum)),size=5) +
    scale_x_discrete(limits = c("Outgoing cell", "Pathway"),
                     expand = c(0.3, 0.1)) +
    theme_bw()
  dev.off()


   ############################################################################
  # incoming alluvial plot
  ############################################################################

  incoming_data$Signaling<-as.character(incoming_data$Signaling)
  incoming_data<-incoming_data[order(incoming_data$Signaling,decreasing = FALSE),]
  incoming_data$Signaling <- factor(incoming_data$Signaling,levels=unique(incoming_data$Signaling))


  #fname=paste0(outFolder,subFolder,"/",xcondition,"/",xcondition,"_incoming_ggalluvial_cuttoff.0.7.pdf");
  
  fname=paste0(outFolder,subFolder,"/",xcondition,"/",xcondition,"_incoming_ggalluvial_selectedPathways_cuttoff_",Sys.Date(),"_.pdf")
  pdf(fname,width=20,height=10)
  ggplot(data = incoming_data,
         aes(axis1 = Signaling, axis2 = CellGroup_incoming, y = Contribution_incoming_signaling)) +
    geom_alluvium() +#aes(fill = Signaling)
    geom_stratum(aes(fill=CellGroup_incoming)) +
    geom_flow(aes(fill=CellGroup_incoming))+
    scale_fill_manual("Cell type",values=cluster.Colors) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),size=5) +
    scale_x_discrete(limits = c("Pathway", "Incoming cell"),
                     expand = c(0.15, 0.05)) +
    theme_bw()
  dev.off()


  }  )

})



########################################################### 
# circle plots per pathways
###########################################################



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


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)

selected_pathways_mat<-read.delim("selectedpathways.txt")
locations<-c("Cervix" ,    "Decidua",    "Myometrium")
conditions<-c("Control","E. coli")
nolabel<-TRUE
sapply(locations, function(xlocation){
  
  
  print(xlocation)
  subFolder<-xlocation
  
  selected_pathways<-selected_pathways_mat[,xlocation]
  selected_pathways<-selected_pathways[selected_pathways!=""]
  
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  
sapply (conditions ,function(xcondition){
    #for (xcondition in conditions){
    print(xcondition)
    
    
    
    cellchat<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_",xcondition,"_","2021-12-21",".rds"))
    pathways<-cellchat@netP$pathways
    
    selected_pathways<-intersect(pathways,selected_pathways)
    groupSize <- as.numeric(table(cellchat@idents))
    
    
    if (xcondition=="E. coli") xcondition<-"Ecoli" else xcondition<-"control"
    
    sapply(selected_pathways,function(pathways.show){
      #pdf(paste0(outFolder,subFolder,"/",xcondition,"/Pathways_circleplots_top50/",pathways.show,".pdf"),width=20,height=20)
      #pdf(paste0(outFolder,subFolder,"/",xcondition,"/Pathways_circleplots_top25/",pathways.show,".pdf"),width=20,height=20)
      
      if (nolabel) 
        {vertexlabelcex=0.000001
        filename<-"/Pathways_circleplots_nolable_top25/"}else {
           vertexlabelcex=2
           filename<-"/Pathways_circleplots_top25/"}
      
      #pdf(paste0(outFolder,subFolder,"/",xcondition,"/Pathways_circleplots_nolable_top25/",pathways.show,".pdf"),width=15,height=15)
      pdf(paste0(outFolder,subFolder,"/",xcondition,filename,pathways.show,".pdf"),width=15,height=15)
      
      system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,filename))
      
      
      par(mfrow=c(1,1))
      netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",color.use=cluster.Colors[rownames(cellchat@net$weight)],vertex.label.cex=vertexlabelcex,top=0.25,arrow.width = 15)
      dev.off()
    })

    
  })
})






########################################################################
# pathways, DEGs , LR pairs 
########################################################################

outFolder="./10_CellChat_comparison_conditions_after_filter_200_plots/"


cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster
clust2Name<-cell.type.annotation$Cell.type
clust2Name<-paste0(c(0:30),"_",clust2Name)
names(clust2Name)<-c(0:30)

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)

selected_pathways_mat<-read.delim("selectedpathways.txt")
locations<-c("Cervix" ,    "Decidua",    "Myometrium")
conditions<-c("Control","E. coli")


pathways_genes_celltypes_db_final<-lapply(locations, function(xlocation){
  
  subFolder<-xlocation
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  #pathways_genes_celltypes_db_list<-lapply (conditions ,function(xcondition){
    
  pathways_genes_celltypes_db_list<-list()
   for (i in 1:length(conditions)){
    xcondition<-conditions[i]
    cellchat<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_",xcondition,"_","2021-12-21",".rds"))
    pathways<-cellchat@netP$pathways
    groupSize <- as.numeric(table(cellchat@idents))
    
    
    if (xcondition=="E. coli") xcondition<-"Ecoli" else xcondition<-"control"
    
    
    system(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))
    print(paste0("mkdir -p ", outFolder,subFolder,"/",xcondition,"/"))

    pathways_genes_celltypes<-lapply(pathways,function(x){
      
      pairLR <- extractEnrichedLR(cellchat, signaling = x, geneLR.return = FALSE)
      LR_genes<-unique(unlist(str_split(pairLR$interaction_name,"_")))
      LR_genes<-tolower(LR_genes)
      
      celltypes<-res %>% filter(tolower(kbid) %in% LR_genes)%>% dplyr::select(Cell_type)%>% unlist %>% unique()
      
      LR_DEGs <-res %>% filter (tolower(kbid) %in% LR_genes) %>% dplyr::select(kbid)%>% unlist %>% unique()
      
      LR_DEGs_len<-length(LR_DEGs)
      LR_DEGs<-paste(LR_DEGs,collapse = ", ")
      
      LR_genes<-paste(LR_genes,collapse = ", ")
      
      celltype_len<-length(celltypes)
      celltypes<-paste(celltypes,collapse = ", ")
      if (celltypes=="") celltypes<-"NA"
      return(c(x,LR_genes,LR_DEGs,celltypes,celltype_len,LR_DEGs_len))
    })
    pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes)
    
    pathways_genes_celltypes_db<-as.data.frame(pathways_genes_celltypes_db)
    
    print(xcondition)
    pathways_genes_celltypes_db$condition<-xcondition
    pathways_genes_celltypes_db$location<-xlocation
    colnames(pathways_genes_celltypes_db)<-c("Pathway","Genes_in_pathway","DEGs_in_pathway","Celltype","celltype_len","LR_DEGs_len","Condition","Location")
    
    pathways_genes_celltypes_db_list[[i]]<-pathways_genes_celltypes_db
    #return(pathways_genes_celltypes_db)
    
  }  #)
  
  pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes_db_list)
  return(pathways_genes_celltypes_db) 
  
})


pathways_genes_celltypes_db<-do.call(rbind,pathways_genes_celltypes_db_final)

pathways_genes_celltypes_db<-pathways_genes_celltypes_db[order(as.numeric(pathways_genes_celltypes_db$celltype_len),decreasing = TRUE),]

write.csv(pathways_genes_celltypes_db,file=paste0(outFolder,"pathways_genes_celltypes_db.csv"))
write_rds(pathways_genes_celltypes_db,file=paste0(outFolder,"pathways_genes_celltypes_db.rds"))



########################################################################
# pathways, information flow
########################################################################


informationflow_list<-lapply(locations,function(xlocation){
  
  
  subFolder<-xlocation
  
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  
  cellchat_Control<-cellchat_Control_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","Control","_","2021-12-21",".rds"))
  cellchat_Ecoli<-cellchat_Ecoli_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","E. coli","_","2021-12-21",".rds"))
  
  
  
  pathways_Control<-cellchat_Control_original@netP$pathways
  controls<-cbind(pathways_Control, rep("control",length(pathways_Control)))
  controls<-as.data.frame(controls)
  colnames(controls)<-c("pathway","control")
  
  
  pathways_Ecoli<-cellchat_Ecoli_original@netP$pathways
  Ecolis<-cbind(pathways_Ecoli, rep("Ecoli",length(pathways_Ecoli)))
  Ecolis<-as.data.frame(Ecolis)
  colnames(Ecolis)<-c("pathway","Ecoli")
  
  allpathways<-Ecolis %>% full_join(controls)
  
  allpathways<-allpathways %>% arrange(desc(!is.na(control),(!is.na(Ecoli))))
  
  write.csv(allpathways,file=paste0(outFolder,subFolder,"/","pathways_",xlocation,".csv"))
  
  shared_pathways<-intersect(pathways_Control,pathways_Ecoli)
  pathways_Control_only<-pathways_Control[!pathways_Control %in% pathways_Ecoli]
  pathways_Ecoli_only<-pathways_Ecoli[! pathways_Ecoli %in% pathways_Control]
  
  
  
  
  # there is an additional population  30_B cell specific to cellchat_Control compared to Ecoli
  # we lift up Ecoli by lifting up the cell groups to the same cell labels as cellchat_Control 
  
  control_cells<-rownames(cellchat_Control@net$prob)
  Ecoli_cells<-rownames(cellchat_Ecoli@net$prob)
  
  Ecoli_cells_only<-Ecoli_cells [!Ecoli_cells %in%control_cells ]
  control_cells_only<-control_cells [! control_cells  %in%Ecoli_cells ]
  
  if (length(Ecoli_cells_only)==0)
  {
    group.new = levels(cellchat_Control@idents)
    cellchat_Ecoli <- liftCellChat(cellchat_Ecoli, group.new)
  }else if (length(control_cells_only)==0){
    group.new = levels(cellchat_Ecoli@idents)
    cellchat_Control <- liftCellChat(cellchat_Control, group.new)} else 
    {
      group.new = union(levels( cellchat_Ecoli@idents), levels(cellchat_Control@idents))
      cellchat_Control <- liftCellChat(cellchat_Control, group.new)
      cellchat_Ecoli <- liftCellChat(cellchat_Ecoli, group.new)
    }
  
  
  # now merge
  object.list <- list(Control = cellchat_Control, Ecoli = cellchat_Ecoli)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  
  
  
  #pdf(paste0(outFolder,subFolder,"/","informationflow_",xlocation,".pdf"),width=10,height=18)
  #gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#333399","#A50021"),font.size=15)
  # gg2 
  # dev.off()
  data<-gg2$data
  colnames(data)<-c("Pathway","contribution","contribution.scaled","Condition","contribution.relative.1","pvalues" )
  data$Location<-xlocation
  #data$Condition[which(data$Condition=="Control")]<-"control"
  data<-data %>% select(Pathway,Location,Condition,contribution,contribution.scaled,pvalues)
  
  return (data)
  })



informationflow_data<-do.call(rbind,informationflow_list)

write.csv(informationflow_data,file=paste0(outFolder,"informationflow_data.csv"))
write_rds(informationflow_data,paste0(outFolder,"informationflow_data.rds"))


informationflow_data<-read_csv("10_CellChat_comparison_conditions_after_filter_200_plots/informationflow_data.csv")

informationflow_data$Condition<-tolower(informationflow_data$Condition)
pathways_genes_celltypes_db$Condition<-tolower(pathways_genes_celltypes_db$Condition)


pathways_genes_celltypes_contribution<-pathways_genes_celltypes_db %>% left_join(informationflow_data)

write.csv(pathways_genes_celltypes_contribution,file=paste0(outFolder,"pathways_genes_celltypes_contribution.csv"))



####################################################################

# selecting top pathways based on information flows
####################################################################
informationflow_data<-read_rds(paste0("10_CellChat_comparison_conditions_after_filter_200_plots/informationflow_data.rds"))

pathways<-sapply(locations, function(xlocation){
  
  subFolder<-xlocation
  system(paste0("mkdir -p ", outFolder,subFolder,"/"))

  #for (xcondition in conditions){
  
  xcondition<-"E. coli"
  
  cellchat<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_",xcondition,"_","2021-12-21",".rds"))
  pathways<-cellchat@netP$pathways
  groupSize <- as.numeric(table(cellchat@idents))
  
  
  if (xcondition=="E. coli") xcondition<-"Ecoli" 
  
  informationflow_data_specific<-informationflow_data %>% filter(Location==xlocation) 
  
  # top 10 %
  #informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" ) %>% arrange(pvalues,-contribution.scaled) 
  
  informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" ) %>% arrange(pvalues) 
  
  #top10<-round(0.1 *nrow(informationflow_data_specific))
  
  top10<-30
  #informationflow_data_specific<-informationflow_data_specific %>% top_n(top10,w=c(pvalues))
  informationflow_data_specific<-informationflow_data_specific %>% head(top10)
  
  
  #informationflow_data_specific<-informationflow_data_specific %>% filter(Condition=="Ecoli" & pvalues<=0.01 ) 
  selected_pathways<-informationflow_data_specific %>% select (Pathway) %>% unlist %>% unique() %>% as.character()
 
return (selected_pathways)
})

