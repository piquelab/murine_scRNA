
#######################################################################
# CellChat analysis 
# receiving and sending strength
#######################################################################



locations<-c("Cervix" ,    "Decidua",    "Myometrium")
conditions<-c("Control" ,"E. coli")


strength_locations<-lapply(locations, function(xlocation){
  
subFolder<-xlocation
  
system(paste0("mkdir -p ", outFolder,subFolder,"/"))
  


cellchat_Control<-cellchat_Control_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","Control","_","2021-12-21",".rds"))
cellchat_Ecoli<-cellchat_Ecoli_original<-read_rds(paste0("./10_CellChat_analysis_default_after_filter_200/","cellchat_",xlocation,"_","E. coli","_","2021-12-21",".rds"))

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



object=cellchat
weight.scale = T
measure = "weight"
top=0.25
color.use=cluster.Colors[rownames(cellchat@netP$Ecoli$prob)]
arrow.width=4
vertex.label.cex = 1
comparison = c(1, 2)
title.name = NULL
sources.use = NULL
targets.use = NULL 
remove.isolate = FALSE 

weight.scale = FALSE
vertex.weight = 20
vertex.weight.max = NULL
vertex.size.max = 15
vertex.label.cex = 1
vertex.label.color = "black"
edge.weight.max = NULL
edge.width.max = 8
alpha.edge = 0.6
label.edge = FALSE
edge.label.color = "black"
edge.label.cex = 0.8
edge.curved = 0.2
shape = "circle"
layout = in_circle() 
margin = 0.2
arrow.width = 1
arrow.size = 0.2

  #measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed",  weighted = T)
  
  sending_receiving_strength<-strength(g)
  sending_receiving_strength<-sending_receiving_strength[order(sending_receiving_strength,decreasing = TRUE)]
  sending_receiving_strength<-as.data.frame(sending_receiving_strength)
  sending_receiving_strength$cell<-rownames(sending_receiving_strength)
  sending_receiving_strength<-sending_receiving_strength[,c(2,1)]
    
  sending_strength<-strength(g, mode="out")
  sending_strength<-sending_strength[order(sending_strength,decreasing = TRUE)]
  sending_strength<-as.data.frame(sending_strength)
  sending_strength$cell<-rownames(sending_strength)
  sending_strength<-sending_strength[,c(2,1)]
  
  receiving_strength<-strength(g, mode="in")
  receiving_strength<-receiving_strength[order(receiving_strength,decreasing = TRUE)]
  
  receiving_strength<-as.data.frame(receiving_strength)
  receiving_strength$cell<-rownames(receiving_strength)
  receiving_strength<-receiving_strength[,c(2,1)]
  
  strength_location<-cbind(sending_receiving_strength,sending_strength,receiving_strength)
  rownames(strength_location)<-NULL
  
  strength_location$location<-xlocation
  return(strength_location)
  
})

outFolder="./10_CellChat_comparison_conditions_after_filter_200_plots/"
strength_locations_mat<-do.call(rbind,strength_locations)
write.csv(strength_locations_mat,file=paste0(outFolder,"differential_strength_locations_mat.csv"))
