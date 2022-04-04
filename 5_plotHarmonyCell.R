library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(future)
library(harmony)

#############################################################
### UMAP plots  

#############################################################


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


outFolder="./5_harmony_cellClass_soupx-doubletfinder_chrM_res0.5_plots/"
system(paste0("mkdir -p ", outFolder))


###########################################

sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Name<-cell.type.annotation$Cell.type
clust2Name<-paste0(c(0:30),"_",clust2Name)
names(clust2Name)<-c(0:30)

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
names(cluster.Colors)<-clust2Name
md<-sc@meta.data

dim(sc)
table(sc$Library)
table(sc$Location) 

cluster.Colors2<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
names(cluster.Colors2)<-c(0:30) #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



fname=paste0(outFolder,"UMAP_Harmony.png");
png(fname,width=1000,height=1000)
#DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) #+ NoLegend()
DimPlot(sc, cols = cluster.Colors2, label=TRUE , repel=TRUE)
dev.off()


# UMAP with colors

aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Library","Location","Condition")) 
aa$cluster_name <- clust2Name[as.character(aa$seurat_clusters)]

fname=paste0(outFolder,"UMAP_LocationHarmony.png")
png(fname,width=800,height=600)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Location)) +
    geom_point(size=0.1) +
    scale_colour_manual(values=c("Cervix"="#F4B183","Decidua"="#DDABD6","Myometrium"="#8AD2CD") )+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    theme(legend.text=element_text(size=30,face="bold"), axis.text=element_text(size=30,face="bold"), axis.title=element_text(size=20,face="bold"))+
    theme_bw()
p1
##    theme_black()
dev.off()


# UMAP with colors
fname=paste0(outFolder,"UMAP_LocationHarmony.png");
png(fname,width=1200,height=800)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
    geom_point(size=0.1) +
    theme(text = element_text(size=30,face = "bold"))+
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"), 
          legend.text=element_text(size=25,face="bold"))+
    #
    theme_bw()
p1
##    theme_black()
dev.off()



source("theme_black.R")
## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony_black.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
    guides(colour = guide_legend(override.aes = list(size=10),title="Condition")) +
    theme(legend.text=element_text(size=30), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))+
    scale_color_manual(values=c("Control"="#333399","E. coli"="#A50021"))+
    theme_black()+
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 30, face = "bold"),
          legend.title=element_text(size=30,face="bold"), 
          legend.text=element_text(size=30,face="bold"))

p1
dev.off()


fname=paste0(outFolder,"UMAP_LocationHarmony.Cell_annotation.png");
png(fname,width=2000,height=1000)
#p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
    geom_point(size=0.1) +
    #theme(legend.text=element_text(size=30), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))+
    scale_color_manual(values=cluster.Colors) +
    #guides(colour = guide_legend(override.aes = list(size=5),title="Cluster")) +
    facet_grid(Condition ~ Location) +
    
    theme_bw()+
    theme(text = element_text(size=25,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"), 
          legend.text=element_text(size=25,face="bold"))
p1
dev.off()


aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Library","Location","Condition")) 
aa$cluster_name <- clust2Name[as.character(aa$seurat_clusters)]

aa<-aa[order(aa$seurat_clusters,decreasing = TRUE),]
aa$seurat_clusters <- factor(aa$seurat_clusters,levels=unique(aa$seurat_clusters))
aa$cluster_name <- factor(aa$cluster_name,levels=unique(aa$cluster_name))

write.csv(table(aa$cluster_name, aa$Location),file=paste0(outFolder,"cell.count.csv"))

fname=paste0(outFolder,"UMAP_LocationCondition.Barplot.pdf");
pdf(fname,width=15,height=6)
#p2 <- ggplot(aa,aes(x=seurat_clusters,fill=Condition))+ #fill=Condition)) +
p2 <- ggplot(aa,aes(x=cluster_name,fill=Condition))+ #fill=Condition)) +
    geom_bar(position = position_stack(reverse = TRUE)) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    facet_grid(.~Location) + coord_flip() +
    scale_fill_manual(values=c("Control"="#333399","E. coli"="#A50021"))+
    xlab("")+
    theme_bw()
p2
dev.off()



fname=paste0(outFolder,"UMAP_ConditionLocation.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=seurat_clusters,fill=Location))+ #fill=Condition)) +
    geom_bar(position = position_stack(reverse = TRUE)) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    facet_grid(.~Condition) + coord_flip() +
    #scale_fill_manual(values=c("Cervix"="","Decidua"="","Myometrium"=""))+
    xlab("")+
    theme_bw()
p2
dev.off()



fname=paste0(outFolder,"UMAP_LocationCondition.overall.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=seurat_clusters,fill=Condition))+ #fill=Condition)) +
    geom_bar(position = position_stack(reverse = TRUE)) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    coord_flip() +
    scale_fill_manual(values=c("Control"="#333399","E. coli"="#A50021"))+
    xlab("")+
    theme_bw()
p2
dev.off()



fname=paste0(outFolder,"UMAP_LocationConditionHarmony.png");
png(fname,width=1400,height=1000)
#p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
    geom_point(size=0.2) +
    scale_color_manual(values=cluster.Colors) +
    #guides(colour = guide_legend(override.aes = list(size=5),title="Cluster")) +
    ##    facet_wrap(~LocTime) +
    facet_grid(Condition ~ Location) +
    theme_bw()
p1
dev.off()



fname=paste0(outFolder,"UMAP_LibraryLocationConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Library)) +
    geom_point(size=0.1) +
    #scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Library")) +
    ##    facet_wrap(~LocTime) +
    facet_grid(Condition ~ Location) +
    theme_bw()
p1
dev.off()


#######################################################

## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.Cell_annotation.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
    facet_wrap(Condition ~ Location) +
    theme_bw()
p1
##    theme_black()
dev.off()



########################################################
library(reshape)

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)

res_up<-res %>% filter(padj<0.1,log2FoldChange>0) %>% dplyr::count(cname)
colnames(res_up)[2]<-"up-regulated"
res_down<-res %>% filter(padj<0.1,log2FoldChange<0) %>% dplyr::count(cname)
colnames(res_down)[2]<-"down-regulated"
res_down[,2]<--1*res_down[,2]
res_total<-res_down  %>% left_join(res_up)


cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Name<-cell.type.annotation$Cell.type
clust2Name<-paste0(c(0:30),"_",clust2Name)
names(clust2Name)<-c(0:30)

res_total <- res_total %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)

res_total$sum<-res_total$`down-regulated`+res_total$`up-regulated`
rw<-as.character(clust2Name[res_total$Cell_type])

res_total<-as.matrix(res_total)
rownames(res_total)<-rw

res_total[is.na(res_total)]<-0
rw<-clust2Name[res_total[,"Cell_type"]]

res_total2<-as.matrix(res_total[,2])
res_total2<-as.matrix(res_total2)

dat1<- melt(res_total2)
dat1$Var1<-as.character(rw)
colnames(dat1)<-c("cluster","var2","location")
res_total<-res_total[,c(-1,-2,-3,-6)]
res_total<-as.matrix(res_total)
res_total<-apply(res_total, c(1,2),as.numeric)
    
library(reshape2)
library(ggplot2)
dat2 <- melt(res_total)
dat2$Var1<-as.character(rw)
colnames(dat2)<-c("cluster","DE","DEnumber")
dat2$Location<-dat1$location

dat2$clusternumber<-as.numeric(sapply(dat2$cluster,function(x){return(unlist(strsplit(as.character(x),"_"))[1])}))
dat2<-dat2[order(dat2$clusternumber,decreasing = TRUE),]


fname=paste0(outFolder,"barplot-DE-up-down.pdf")
pdf(fname,width=15,height=7)
ggplot(dat2, aes(x=reorder(cluster,-clusternumber), y=DEnumber,fill=DE))+#,color=mycolor)+#,fill=DE) +
    geom_bar(stat='identity') +
    theme_bw()+
    scale_fill_manual("legend",values=c("down-regulated"="#F2B8C6","up-regulated"="#990F02") )+
    facet_grid(.~Location,scales="free") + coord_flip() +
    theme(legend.title=element_blank())+
    ylab("Number of differentially expressed genes")+
    xlab("")
dev.off() 


#######################
# cell counts per location
cell_counts_location<-as(table(aa$cluster_name, aa$Location),"matrix")
cell_counts_location<-as.data.frame(cell_counts_location)
cell_counts_location<-cell_counts_location %>% select(Cervix_cellcount=Cervix, Decidua_cellcount=Decidua, Myometrium_cellcount=Myometrium)

cell_counts_location$cluster.seurat<-unlist(strsplit(rownames(cell_counts_location),"_"))[seq(1, 2* length(rownames(cell_counts_location)),by=2)]
cell_counts_location$cluster.seurat<-as.integer(cell_counts_location$cluster.seurat)
all_annotations_joined<-read.csv(paste0("5_cellTypeAnnotation/all_annotations_joined.csv"))

all_annotations_joined<-all_annotations_joined %>% inner_join(cell_counts_location)
write.csv(all_annotations_joined,file=paste0("5_cellTypeAnnotation/all_annotations_joined.csv"),row.names = FALSE,quote = FALSE)



for (loc in unique(aa$Location))
{
    aa2<-aa %>% filter(loc==Location)
    #table(aa2$cluster_name, aa$Condition)
    write.csv(table(aa2$cluster_name, aa2$Condition),file=paste0(outFolder,loc,".cell.count.csv"))
    
    }

