
#######################################################################
# CellChat plots 
# Comprison between E.coli and control - differential plots, vector plots, circle plots
#######################################################################


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(tidyverse)
library(Seurat)


##########################################################################
# Whether the cell-cell communication is enhanced or not
# The interaction between which cell types is significantly changed
# How the major sources and targets change from one condition to another
##########################################################################


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



res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)



cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster

tobechanged<-c("5_Epithelial-1 (Basal)","28_Epithelial-9 (Secretory)","7_Epithelial-2 (Squamous)","8_Epithelial-3 (Squamous)","10_Epithelial-4 (Glandular)","13_Epithelial-5 (Luminal)" , "14_Epithelial-6 (Secretory)" ,"23_Epithelial-8 (Enterocyte)",
              "11_Epithelial-10 (Proliferative)" ,"20_Epithelial-7 (Glandular)")


outFolder="./10_CellChat_comparison_conditions_after_filter_200_plots/"

system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


####################################
# Load data
####################################

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)

names(cluster.Colors)<-clust2Name #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")


locations<-c("Cervix" ,    "Decidua",    "Myometrium")
conditions<-c("Control" ,"E. coli")


##################################################################
# circle plots for paper 
##################################################################


# sapply(locations, function(xlocation){
#   
#   
#   
#   subFolder<-xlocation
#   
# 
#   cellchat_Control<-cellchat_Control_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","Control","_","2021-12-21",".rds"))
#   cellchat_Ecoli<-cellchat_Ecoli_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","E. coli","_","2021-12-21",".rds"))
#   # 
# 
#   pathways_Control<-cellchat_Control_original@netP$pathways
#   controls<-cbind(pathways_Control, rep("control",length(pathways_Control)))
#   controls<-as.data.frame(controls)
#   colnames(controls)<-c("pathway","control")
# 
#   # there is an additional population  30_B cell specific to cellchat_Control compared to Ecoli
#   # we lift up Ecoli by lifting up the cell groups to the same cell labels as cellchat_Control 
#   
#   control_cells<-rownames(cellchat_Control@net$prob)
#   Ecoli_cells<-rownames(cellchat_Ecoli@net$prob)
#   
#   Ecoli_cells_only<-Ecoli_cells [!Ecoli_cells %in%control_cells ]
#   control_cells_only<-control_cells [! control_cells  %in%Ecoli_cells ]
#   
#   if (length(Ecoli_cells_only)==0)
#   {
#     group.new = levels(cellchat_Control@idents)
#     cellchat_Ecoli <- liftCellChat(cellchat_Ecoli, group.new)
#   }else if (length(control_cells_only)==0){
#     group.new = levels(cellchat_Ecoli@idents)
#     cellchat_Control <- liftCellChat(cellchat_Control, group.new)} else 
#     {
#       group.new = union(levels( cellchat_Ecoli@idents), levels(cellchat_Control@idents))
#       cellchat_Control <- liftCellChat(cellchat_Control, group.new)
#       cellchat_Ecoli <- liftCellChat(cellchat_Ecoli, group.new)
#     }
#   
#   
#   # now merge
#   object.list <- list(Control = cellchat_Control, Ecoli = cellchat_Ecoli)
#   cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#   
#   
#     color.use.Ecoli=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)]
#     color.use.Control=cluster.Colors[rownames(cellchat@netP$Control$prob)]
#     coloruses<-list(color.use.Ecoli,color.use.Control)
#   
#   
#     system(paste0("mkdir -p ", outFolder,subFolder,"/nolabel","/"))
#   
#     pdf(paste0(outFolder,subFolder,"/nolabel/","diffInteraction_",xlocation,".pdf"),width=15,height=15)
#     gg2<-netVisual_diffInteraction(cellchat, vertex.weight=6,vertex.label.cex=0.000001, weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)],arrow.width=15)
#     
#     gg2
#     dev.off()
#   
#     pdf(paste0(outFolder,subFolder,"/","diffInteraction_",xlocation,".pdf"),width=15,height=15)
#     gg2<-netVisual_diffInteraction(cellchat, vertex.weight=6, weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)],arrow.width=15)
#     gg2
#     dev.off()
#     
#     
#   
#   
#   
# })


##################################################
#comparison between ecoli and pbs plots
######################################################

sapply(locations, function(xlocation){
  
  
  
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
  
  color.use.Ecoli=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)]
  color.use.Control=cluster.Colors[rownames(cellchat@netP$Control$prob)]
  coloruses<-list(color.use.Control,color.use.Ecoli)
  
  overalcoloruses<-cluster.Colors[unique(c(rownames(cellchat@netP$Ecoli$prob), rownames(cellchat@netP$Control$prob)))]
  
  
#   
  pdf(paste0(outFolder,subFolder,"/","informationflow_",xlocation,".pdf"),width=10,height=18)
  #gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#333399","#A50021"),font.size=15)
  gg2
  dev.off()
#   
#   
#   
#   
  # Compare the total number of interactions and interaction strength
  # pdf(paste0(outFolder,subFolder,"/","total_number_interactions_",xlocation,".pdf"),width=10,height=4)
  # gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use=c("#A50021","#333399"))
  # gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use=c("#A50021","#333399"))
  # gg1 + gg2
  # dev.off()
#   
#   
#   # Compare the number of interactions and interaction strength among different cell populations
#   
#   # The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, 
#  # where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#   # 
#   
  
  
  pdf(paste0(outFolder,subFolder,"/","diffInteraction_",xlocation,".pdf"),width=24,height=24)
  gg2<-netVisual_diffInteraction(cellchat, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)],vertex.label.cex = 1,edge.width.max=20)
  gg2
  dev.off()
  
  pdf(paste0(outFolder,subFolder,"/","diffInteraction_nolabel_",xlocation,".pdf"),width=24,height=24)
  gg2<-netVisual_diffInteraction(cellchat, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,weight.scale = T, measure = "weight",top=0.25,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)],vertex.label.cex = 0.00001,edge.width.max=20)
  gg2
  dev.off()
  
  
  
  weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  pdf(paste0(outFolder,subFolder,"/","interactions_strength_conditions_nolabel_",xlocation,".pdf"),width=26,height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight, arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8, vertex.label.cex=0.000001,top=0.25,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  }
  dev.off()
  
  
  weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  pdf(paste0(outFolder,subFolder,"/","interactions_strength_conditions_",xlocation,".pdf"),width=26,height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight,arrow.width = 2,arrow.size=1, vertex.weight = 15, vertex.size.max=8,top=0.25,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  }
  dev.off()
  
  
  # weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
  # pdf(paste0(outFolder,subFolder,"/","interactions_strength_conditions_",xlocation,".pdf"),width=25,height=20)
  # par(mfrow = c(1,2), xpd=TRUE)
  # for (i in 1:length(object.list)) {
  #   netVisual_circle(object.list[[i]]@net$weight, top=0.25,arrow.width=4,weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(object.list)[i]),color.use = coloruses[[i]])
  # }
  # dev.off()
  
  
#   
#   
#   #We can also show differential number of interactions or interaction strength in a greater details using a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, 
#   # red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.
#   
#   
  pdf(paste0(outFolder,subFolder,"/","diffInteraction_heatmap_",xlocation,".pdf"),width=10,height=10)
  #gg1 <- netVisual_heatmap(cellchat,color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)])
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)])
  #> Do heatmap based on a merged object
  #gg1 + gg2
  gg2
  dev.off()
#   
#   
#   
#  # To better control the node size and edge weights of the inferred networks across different datasets, 
#   # we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.
#   
#   
#   
#   
#   
  color.use.Ecoli=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)]
  color.use.Control=cluster.Colors[rownames(cellchat@netP$Control$prob)]
  coloruses<-list(color.use.Ecoli,color.use.Control)

  
pdf(paste0(outFolder,subFolder,"/","outgoing_incoming_conditions_",xlocation,".pdf"),width=19,height=10)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.listi <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")

  gg[[i]] <- netAnalysis_signalingRole_scatter(object.listi, title = names(object.list)[i], weight.MinMax = weight.MinMax,color.use = coloruses[[i]])
}
 p<-patchwork::wrap_plots(plots = gg)
 plot(p)
 dev.off()

label.size = 4
dot.size = c(2, 6)
dot.alpha = 0.6
font.size.title = 13
font.size = 15
df1<-gg[[1]]$data
df2<-gg[[2]]$data
newdf<-rbind(df1,df2)

newdf$color_custome<-cluster.Colors[newdf$labels]

xlabel = "Outgoing interaction strength"
ylabel = "Incoming interaction strength"

newdf<-newdf %>% arrange(labels)
# newdf$color_custome <- factor(newdf$color_custome,levels=unique(newdf$color_custome))
newdf$labels <- factor(newdf$labels,levels=unique(newdf$labels))
overalcoloruses<-overalcoloruses[newdf$labels]

require(grid)
gg <- ggplot(data = newdf, aes(x, y),show.legend = F) + geom_point(aes(size = Count, colour = labels, fill = labels),show.legend = F)
gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in")) + labs(title = "Ecoli + control", x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, face = "plain")) + theme(axis.line.x = element_line(size = 0.25),                                                                                                                                                                               axis.line.y = element_line(size = 0.25))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)
gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), size = label.size, show.legend = F, segment.size = 0.2, segment.alpha = 0.5)

gg<-gg+geom_line(aes(group = labels,colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_colour_manual(values =overalcoloruses, drop = FALSE) + guides(colour = FALSE)+ guides(fill = FALSE)

pdf(paste0(outFolder,subFolder,"/","outgoing_incoming_conditions_both_",xlocation,".pdf"),width=17,height=12)
gg+ theme(legend.position = "none")
dev.off()


require(grid)
#newdf <- data.frame(x = outgoing.cells, y = incoming.cells, labels = names(incoming.cells), Count = num.link)
#newdf$labels <- factor(newdf$labels, levels = names(incoming.cells))
gg <- ggplot(data = newdf, aes(x, y)) + geom_point(aes(size = Count, colour = labels, fill = labels),show.legend = FALSE)
gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in")) + labs(title = "Ecoli + control", x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, face = "plain")) + theme(axis.line.x = element_line(size = 0.25),axis.line.y = element_line(size = 0.25))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)+ guides(fill=FALSE, color=FALSE)
gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), size = 0.0001, show.legend = F, segment.size = 0.2, segment.alpha = 0.5)

gg<-gg+geom_line(aes(group = labels,colour=labels,size = 1.5),arrow = arrow(type = "closed",length=unit(0.2, "inches")))
#gg <- gg + scale_fill_manual(values = ggplot2::alpha(overalcoloruses, alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
gg <- gg + scale_colour_manual(values =overalcoloruses, drop = FALSE) + guides(colour = FALSE)+ guides(fill = FALSE)
#gg + theme(legend.position = "none")
pdf(paste0(outFolder,subFolder,"/","outgoing_incoming_conditions_nolabel_both_",xlocation,".pdf"),width=17,height=12)
gg+ theme(legend.position = "none")

dev.off()


})







