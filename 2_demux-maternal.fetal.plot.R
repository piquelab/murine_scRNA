######################################
### maternal fetal umap plot ###
######################################
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster
clust2Name<-cell.type.annotation$Cell.type
clust2Name<-paste0(c(0:30),"_",clust2Name)
names(clust2Name)<-c(0:30)


dd <- read_rds("./1_demux_output/1_demux_New.ALL.rds") %>% filter(DIFF.LLK.BEST.NEXT>20 & NUM.SNPS>200 & NUM.READS>150 ) %>%
  select(NEW_BARCODE,SNG.BEST.GUESS,DROPLET.TYPE,NUM.SNPS,NUM.READS,DIFF.LLK.BEST.NEXT,EXP)

sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")

md<-sc@meta.data

md = md %>% rownames_to_column("BARCODE") %>% mutate(NEW_BARCODE=paste0(Library,"_", gsub("-1_.*","",BARCODE))) %>% left_join(dd) 


sc@meta.data <- md
rownames(md)<-md$BARCODE

aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","SNG.BEST.GUESS")) 

aa$seurat_clusters <- md[rownames(aa),"seurat_clusters"]
aa$seurat_clusters <- clust2Name[aa$seurat_clusters]

aa$SNG.BEST.GUESS<-md[rownames(aa),"SNG.BEST.GUESS"]
aa$Location<-md[rownames(aa),"Location"]
write.csv(table(aa$seurat_clusters,aa$SNG.BEST.GUESS),file=paste0(outFolder,"cluster.origin_t120.csv"))


aa<-aa %>% filter(!is.na(SNG.BEST.GUESS))
outFolder <- "./1_demux_output/"

fname=paste0(outFolder,"UMAP_LocationHarmony.Origin_",Sys.Date(),".pdf");
pdf(fname,width=10,height=4)
p2 <- ggplot(aa ,aes(UMAP_1,UMAP_2,color=SNG.BEST.GUESS)) +
  geom_point(size=0.1,alpha = 1) +
  ##    scale_color_manual(values=group.colors) +
  #guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
  scale_color_manual("Origin",values=c("maternal"="#D1D1D1","fetal"="#A61BB5"))+
  facet_wrap(~Location) +
  theme_bw()
p2
##    theme_black()
dev.off()

fname=paste0(outFolder,"UMAP_LocationHarmony.Origin_",Sys.Date(),".alpha0.05.png");
png(fname,width=1000,height=500)
p2 <- ggplot(aa ,aes(UMAP_1,UMAP_2,color=SNG.BEST.GUESS)) +
  geom_point(size=0.1,alpha = 0.05) +
  ##    scale_color_manual(values=group.colors) +
  #guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
  scale_color_manual("Origin",values=c("maternal"="#D1D1D1","fetal"="#A61BB5"))+
  facet_wrap(~Location) +
  theme_bw()+
theme(text = element_text(size=20,face = "bold"),
      plot.title = element_text(size = 20, face = "bold"),
      legend.title=element_text(size=20,face="bold"), 
      legend.text=element_text(size=20,face="bold"))
p2
##    theme_black()
dev.off()


fname=paste0(outFolder,"UMAP.Origin_",Sys.Date(),".pdf");
pdf(fname,width=5,height=4)
p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=SNG.BEST.GUESS)) +
  geom_point(size=0.1,alpha = 1) +
  ##    scale_color_manual(values=group.colors) +
  #guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
  #scale_color_manual(values=c("fetal"="#333399","E. coli"="#maternal"))+
  scale_color_manual("Origin",values=c("maternal"="#D1D1D1","fetal"="#A61BB5"))+
  #facet_wrap(~Location) +
  theme_bw()+
  theme(text = element_text(size=20,face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20,face="bold"))

p2
##    theme_black()
dev.off()


fname=paste0(outFolder,"UMAP.Origin_",Sys.Date(),".alpha0.05");
png(fname,width=900,height=700)
p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=SNG.BEST.GUESS)) +
  geom_point(size=0.1,alpha = 0.05) +
  ##    scale_color_manual(values=group.colors) +
  #guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
  #scale_color_manual(values=c("fetal"="#333399","E. coli"="#maternal"))+
  scale_color_manual("Origin",values=c("maternal"="#D1D1D1","fetal"="#A61BB5"))+
  #facet_wrap(~Location) +
  theme_bw()+
  theme(text = element_text(size=20,face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20,face="bold"))

p2
##    theme_black()
dev.off()