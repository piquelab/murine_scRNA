#############################################################
###  Heatmap plot showing the comparison between the clusters based of log2fc
### 
#############################################################

library(gplots)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)


outFolder="6_celltype_cor_heatmap_plot/"
system(paste0("mkdir -p ", outFolder))

# cell type labels
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


# union of DE genes across cell types
res_de <-read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/SIG.combined.2021-02-17.tsv")
res_de <- res_de %>% separate(cname,c("Location","Cluster"),sep="_",remove=FALSE)
de_gene<-unique(res_de$kbid)
res_de$Cell_type<-clust2Names[res_de$Cluster]
res_de$cname<-paste0(res_de$Cell_type,"_",res_de$Location)  

# calculating correlation based on union of DE genes 
res <-read_tsv("7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res  %>% filter (kbid %in% de_gene)
res <- res %>% separate(cname,c("Location","Cluster"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location)  
clusters<-unique(res_de$cname)


cor_matrix<-matrix(NA,nrow=length(clusters),ncol=length(clusters))
for (i in 1:length(clusters))
{
  for (j in 1: length(clusters))
  {
      resi<-res %>% dplyr::filter(cname == clusters[i]) %>% dplyr::select(kbid,log2FoldChange,lfcSE)
      resj<-res %>% dplyr::filter(cname == clusters[j])%>% dplyr::select(kbid,log2FoldChange,lfcSE)
      colnames(resj)<-c("kbid","log2FoldChange2","lfcSE2")
      res_intersect<-resi%>% inner_join(resj)
      res_intersect<-res_intersect%>% dplyr::select(log2FoldChange,log2FoldChange2)
      if(nrow(res_intersect)>5)
      {
        cr<-cor(res_intersect)[1,2]
        cor_matrix[i,j]<-cr
        cor_matrix[j,i]<-cr
      }
      }
  }
 
rownames(cor_matrix)<-clusters
colnames(cor_matrix)<-clusters
write_rds(cor_matrix,file=paste0(outFolder,"cor_matrix.rds"))



# heatmap plot

fname=paste0(outFolder,"heatmap_celltype_cor.pdf");
pdf(fname,width=30,height=30)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#333399", "white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))

#myBreaks <- c(seq( min(cor_matrix),0,length.out=20), seq(0,max(cor_matrix), length.out=30) )


pheatmap(cor_matrix,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=20)
dev.off()


fname=paste0(outFolder,"heatmap_celltype_cor2.pdf");
pdf(fname,width=30,height=30)

pheatmap(cor_matrix,cluster_rows=TRUE,fontsize=20)
dev.off()



### new plot #### 

fname=paste0(outFolder,"heatmap_celltype_cor3.pdf");
pdf(fname,width=38,height=28)
#pheatmap(cor_matrix,cluster_rows=TRUE,cluster_cols=TRUE,scale="none")
#my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 201)
paletteLength<-30
# my_palette <- colorRampPalette(colors = c("#333399", "white", "#A50021"))(n = paletteLength+1)
#  myBreaks <- c(seq(min(cor_matrix),0, length.out=ceiling(paletteLength/2) + 1), 
#                seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))

my_palette <- colorRampPalette(colors = c("#D8EAF0","white", "#A50021"))(n = paletteLength+1)
myBreaks <- c(seq(min(cor_matrix),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(cor_matrix)/paletteLength, max(cor_matrix), length.out=floor(paletteLength/2)))

pheatmap(cor_matrix,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize =10)
dev.off()


############



clusters_location<-sapply(clusters,function(x){
  y<-unlist(strsplit(x,"_"))
  return(y[length(y)])

})


clusters_number<-sapply(clusters,function(x){
  y<-unlist(strsplit(x,"_"))
  return(as.numeric(y[1]))
  
})

clusters_number<-clusters_number[order(clusters_number,decreasing = FALSE)]

cor_matrix<-cor_matrix[names(clusters_number),]
cor_matrix_cervix<-cor_matrix[,clusters_location[colnames(cor_matrix)]=="Cervix"]



########################################################
# Same cell type across locations
# Decidua-Cervix
# Decidua-Myometrium
# Myometrium-Cervix
########################################################
D<-names(clusters_location)[which(clusters_location =="Decidua")]
M<-names(clusters_location)[which(clusters_location =="Myometrium")]
C<-names(clusters_location)[which(clusters_location =="Cervix")]

#D_C<-names(clusters_location)[which(clusters_location !="Myometrium")]
cor_matrix_D_C<-cor_matrix[D,C]

rw<-unlist(str_split(rownames(cor_matrix_D_C),"_"))[seq(2,3*nrow(cor_matrix_D_C),by=3)]
names(rw)<-rownames(cor_matrix_D_C)
rw<-rw[order(rw)]
cl<-unlist(str_split(colnames(cor_matrix_D_C),"_"))[seq(2,3*ncol(cor_matrix_D_C),by=3)]
names(cl)<-colnames(cor_matrix_D_C)
cl<-cl[order(cl)]
cl<-cl[cl %in% intersect(rw,cl)]
rw<-rw[rw %in% intersect(rw,cl)]
cor_matrix_D_C<-cor_matrix_D_C[names(rw),names(cl)]
cor_matrix_D_M<-cor_matrix[D,M]


rw<-unlist(str_split(rownames(cor_matrix_D_M),"_"))[seq(2,3*nrow(cor_matrix_D_M),by=3)]
names(rw)<-rownames(cor_matrix_D_C)
rw<-rw[order(rw)]
cl<-unlist(str_split(colnames(cor_matrix_D_M),"_"))[seq(2,3*ncol(cor_matrix_D_M),by=3)]
names(cl)<-colnames(cor_matrix_D_M)
cl<-cl[order(cl)]
cl<-cl[cl %in% intersect(rw,cl)]
rw<-rw[rw %in% intersect(rw,cl)]
cor_matrix_D_M<-cor_matrix_D_M[names(rw),names(cl)]



#M_C<-names(clusters_location)[which(clusters_location !="Decidua")]
cor_matrix_M_C<-cor_matrix[M,C]

rw<-unlist(str_split(rownames(cor_matrix_M_C),"_"))[seq(2,3*nrow(cor_matrix_M_C),by=3)]
names(rw)<-rownames(cor_matrix_M_C)
rw<-rw[order(rw)]
cl<-unlist(str_split(colnames(cor_matrix_M_C),"_"))[seq(2,3*ncol(cor_matrix_M_C),by=3)]
names(cl)<-colnames(cor_matrix_M_C)
cl<-cl[order(cl)]
cl<-cl[cl %in% intersect(rw,cl)]
rw<-rw[rw %in% intersect(rw,cl)]
cor_matrix_M_C<-cor_matrix_M_C[names(rw),names(cl)]



dc<-cbind(rep("Decidua-Cervix",length(as.vector(cor_matrix_D_C))),as.vector(cor_matrix_D_C))
dm<-cbind(rep("Decidua-Myometrium",length(as.vector(cor_matrix_D_M))),as.vector(cor_matrix_D_M))
mc<-cbind(rep("Myometrium-Cervix",length(as.vector(cor_matrix_M_C))),as.vector(cor_matrix_M_C))
data_locations<-rbind(dc,dm,mc)
colnames(data_locations)<-c("Location","Correlation")
data_locations<-as.data.frame(data_locations)
data_locations$Correlation<-as.numeric(data_locations$Correlation)


pdf(paste0(outFolder,"boxplot_celltype_locations.pdf"),width=3,height = 3)
ggplot(data_locations,aes(x=Location, y=Correlation,fill=Location)) + 
geom_boxplot()+
xlab("")+
scale_fill_manual(values=c("Decidua-Cervix"="#74AF81","Decidua-Myometrium"="#B46FB5","Myometrium-Cervix"="#99BE7F") )+
theme_bw()+
theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.2))
dev.off()


########################################################
# within location across cell types

cor_matrix_M<-cor_matrix[M,M]
cor_matrix_D<-cor_matrix[D,D]
cor_matrix_C<-cor_matrix[C,C]


library(pheatmap)

fname=paste0(outFolder,"heatmap_celltype_cor_Myometrium.pdf");
pdf(fname,width=25,height=25)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#333399", "white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix_M), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix_M)/paletteLength, max(cor_matrix_M), length.out=floor(paletteLength/2)))
theme_bw()+
pheatmap(cor_matrix_M,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=30)
dev.off()


fname=paste0(outFolder,"heatmap_celltype_cor_Cervix.pdf");
pdf(fname,width=25,height=25)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#333399", "white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix_C), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix_C)/paletteLength, max(cor_matrix_C), length.out=floor(paletteLength/2)))
theme_bw()+
pheatmap(cor_matrix_C,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=30)
dev.off()


fname=paste0(outFolder,"heatmap_celltype_cor_Decidua.pdf");
pdf(fname,width=25,height=25)
paletteLength<-50
my_palette <- colorRampPalette(colors = c("#333399", "white", "#A50021"))(n = paletteLength)
myBreaks <- c(seq(min(cor_matrix_D), 0, length.out=ceiling(paletteLength/2) ), 
              seq(max(cor_matrix_D)/paletteLength, max(cor_matrix_D), length.out=floor(paletteLength/2)))
theme_bw()+
pheatmap(cor_matrix_D,cluster_rows=TRUE,color=my_palette,scale="none",breaks=myBreaks,fontsize=30)
dev.off()


########################################################

# Comparison within each location

########################################################


data<-array()
cell_types<-c()
for (i in 1: length(clust2Names))
{
  locs<-which(clusters_number==as.numeric(names(clust2Names)[i]))
  
  if (length(locs)==3)
  {
    y<-as.vector(cor_matrix[names(locs),names(locs)])
    cell_types<-c(cell_types,clust2Names[i])
    y<-y[c(2,3,6)]
    data<-cbind(data,y)
  }
    
  
}
data<-data[,-1]
colnames(data)<-cell_types

dat1<- melt(data)
colnames(dat1)<-c("var1","Celltype","Correlation")
dat1$Celltype<-as.character(dat1$Celltype)
dat1$Clusternumber<-as.numeric(names(clust2Name)[which(clust2Name %in% dat1$Celltype)])
dat1<-dat1[order(dat1$Clusternumber,decreasing = FALSE),]
dat1$Celltype <- factor(dat1$Celltype,levels=unique(dat1$Celltype))
dat1$Clusternumber <- factor(dat1$Clusternumber,levels=unique(dat1$Clusternumber))


pdf(paste0(outFolder,"boxplot_celltype.pdf"),width=5,height = 4)
ggplot(dat1,aes(x=Celltype, y=Correlation,fill=Celltype)) + 
  geom_boxplot()+
  scale_fill_manual(values=cluster.Colors ) +
  xlab("")+
  theme_bw()+
  #axis.text=element_text(size=20)
  theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.2))
dev.off()



cor_matrix_M<-cor_matrix[M,M]
m<-lower.tri(cor_matrix_M, diag = TRUE)
cor_matrix_M[m]<-NA
d<-as.vector(cor_matrix_M)
d<-d[!is.na(d)]

m_d<-cbind(rep("Myometrium",length(d)),d)

cor_matrix_D<-cor_matrix[D,D]

m<-lower.tri(cor_matrix_D, diag = TRUE)
cor_matrix_D[m]<-NA
d<-as.vector(cor_matrix_D)
d<-d[!is.na(d)]
d_d<-cbind(rep("Decidua",length(d)),d)

cor_matrix_C<-cor_matrix[C,C]
m<-lower.tri(cor_matrix_C, diag = TRUE)
cor_matrix_C[m]<-NA
d<-as.vector(cor_matrix_C)
d<-d[!is.na(d)]

c_d<-cbind(rep("Cervix",length(d)),d)

data2<-rbind(c_d,d_d,m_d)
colnames(data2)<-c("Location","Correlation")
data2<-as.data.frame(data2)
data2$Correlation<-as.numeric(data2$Correlation)

data2<-data2 %>% arrange(desc(Location))
data2$Location[which(data2$Location=="Myometrium")]<-"Uterus"
data2$Location <- factor(data2$Location, levels = unique(data2$Location))



pdf(paste0(outFolder,"boxplot_location.pdf"),width=5,height = 4)
ggplot(data2,aes(x=Location, y=Correlation,fill=Location)) + 
  geom_boxplot()+
  #scale_fill_manual("legend",values=c("Cervix"="#33CC33","Decidua"="#7030A0","Myometrium"="#FF6699") )+
  scale_fill_manual("legend",values=c("Uterus"="#8AD2CD" ,"Decidua"="#DDABD6","Cervix"="#F4B183") )+
  #scale_fill_manual(values=cluster.Colors ) +
  xlab("")+
  theme_bw()+
  #axis.text=element_text(size=20)
  theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.2))
#scale_color_manual(values=c("red", "green", "blue"))
#facet_wrap(~Location, scale="free")
dev.off()
