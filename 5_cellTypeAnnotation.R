
########################################################
# Cell type annotation I
# Mapping based on cell markers from Cellranger
########################################################

library(Matrix)
library(tidyverse)

library(future)



m2 = read_tsv("./5_harmonyClustersDGE-soupx-doubletfinder-newfilter_chrM_res0.5/ClusterDEG.tsv")


Hmax=log2(max(m2$cluster)+1)

m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
  group_by(gene) %>%
  mutate(H=log2(length(cluster))) %>%
  filter(H<=1) %>%
  ungroup()

table(m3$cluster)

top20 <- m3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% ungroup()
table(top20$cluster)


m2 <- markers %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)
top20 <- m2 %>% top_n(n = 20, wt = avg_log2FC)

sig_cellranger = read_csv("./cellranger_files/Singe-cell_analysis_mouse_cellranger.csv")
colnames(sig_cellranger)<-c("Cluster_cellranger","Celltype_cellranger","Gene","Name","Function")

sig_cellrangerj<-sig_cellranger %>% select(gene=Gene,Cluster_cellranger,Celltype_cellranger)
intersected<-top20 %>% inner_join(sig_cellrangerj) 
intersected<-intersected  %>% arrange(cluster) %>% group_by(cluster)
write.csv(intersected,file="intersected_Cellranger_Seurat.csv")

