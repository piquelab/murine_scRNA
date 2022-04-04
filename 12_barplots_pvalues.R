library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)


#######################################################################
#cell count / DEGs comparison between E.coli and control
#######################################################################



outFolder="12_comparison_DEGs_cellcounts/"
system(paste0("mkdir -p ", outFolder))


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

###########################################
# annotation
cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster
clust2Name<-cell.type.annotation$Cell.type
clust2Name<-paste0(c(0:30),"_",clust2Name)
names(clust2Name)<-c(0:30)


# seurat obj
sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")


#######################################
# cell count compariosn E.coli vs pbs 
#######################################
aa <- FetchData(sc,c("seurat_clusters","Location","Condition","status","Library")) 
aa$seurat_clusters <- clust2Name[aa$seurat_clusters]
aa<-aa %>% mutate(Pregnancy_ID=str_match(Library,"[0-9]+"))

# per location
ccdiff_all<-lapply( unique(aa$Location),function(locationx){
#locationx <-"Myometrium"
  aa2<-aa %>% filter(Location==locationx)
  cc <- aa2 %>% group_by(seurat_clusters,Condition,Pregnancy_ID) %>%
    summarize(n=n()) %>%
    group_by(seurat_clusters) %>% mutate(p0=sum(n)/nrow(aa2)) %>%
    group_by(Pregnancy_ID) %>%
    mutate(nt=sum(n),p=n/nt,z=(p-p0)/sqrt(p0*(1-p0)*(1/nt+1/nrow(aa2))))
  
  cluster_filter<-cc %>% group_by(seurat_clusters)%>% summarise(count_control=count(Condition=="Control"), count_ecoli=count(Condition=="E. coli"))%>% filter (count_control>=2 & count_ecoli>=2) %>% select (seurat_clusters) %>% unlist %>% unique()
  
  
  ccdiff <-  cc %>% filter(seurat_clusters%in% cluster_filter) %>% group_by(seurat_clusters) %>% summarize(wilcox.pval = wilcox.test(p ~ Condition)$p.value, 
                                                                                                           test.t = t.test(z ~ Condition)$statistic) %>% ungroup() %>% mutate(wilcox.padj = p.adjust(wilcox.pval))  
  
  
  ccdiff_filtered<-ccdiff %>% filter(wilcox.padj<0.1)
  #ccdiff$Location<-locationx
  #  return(ccdiff_filtered)
  write.csv(ccdiff,file=paste0(outFolder,"/wilcox_result_cellcount_",locationx,".csv"))
}
)




######################################################################
# up-regulated vs down-regulated 
######################################################################
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res<-res %>% filter(padj<0.1)
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
sc@meta.data$cluster_name <- clust2Names[sc@meta.data$seurat_clusters]

binom.test.locs<-lapply( unique(aa$Location),function(locationx){
  res_loc<-res %>% filter (Location==locationx)
  resDE<-res_loc %>% filter(padj<0.1 & !is.na(padj)) #single cell fdr 0.1
  total<-table(resDE$Cell_type)
  res_up<-resDE%>%filter(log2FoldChange>0)
  upregulated<-rep(0,length(total))
  names(upregulated)<-names(total)
  upregulated[names(table(res_up$Cell_type))]<-table(res_up$Cell_type)
  binom.test.res<-c()
  
  sc_loc<-subset(sc,Location==locationx)
  cell_counts<-table(sc_loc$cluster_name)
  
  for( x in unique(resDE$Cell_type))
  {
    
    btest<-binom.test(upregulated[x],total[x],0.5)
    if(is.na(upregulated[x]))
      upregulated[x]<-0
    rb<-as.numeric(c(cell_counts[x],upregulated[x], (total[x]-upregulated[x]),total[x],btest$p.value))
    binom.test.res<-rbind(binom.test.res,rb)
    print(x)
    print(btest$p.value)
  }
  
  rownames(binom.test.res)<-unique(resDE$Cell_type)
  colnames(binom.test.res)<-c("Cell-counts","Up-regulated","Down-regulated","Total","P-value")
  binom.test.res<-as.data.frame(binom.test.res)
  binom.test.res$padj<-p.adjust(binom.test.res$`P-value`,"fdr")
  binom.test.res<-binom.test.res[order(binom.test.res[,"padj"],decreasing = FALSE),]
  binom.test.res$Location<-locationx
  write.csv(binom.test.res,file=paste0(outFolder,"binom.test.DEGs.",locationx,".csv"))
  return(binom.test.res)
  })

binom.test.all<-do.call(rbind,binom.test.locs)
