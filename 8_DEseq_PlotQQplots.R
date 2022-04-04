########################################################
# QQ-plot

########################################################

library(tidyverse)
library(DESeq2)
library(qqman)
library(org.Hs.eg.db)
library(clusterProfiler)

outFolder <- paste0("./8_outputs_DESeq-newfilter-Plots_res0.5/")
system(paste0("mkdir -p ",outFolder))



sc <- readRDS("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(padj))


## cluster colors 
cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster
cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
names(cluster.Colors)<-clust2Names 
res$Cell_type<-clust2Names[res$Cell_type]
res$cluster.Colors<-cluster.Colors[res$Cell_type]





res2 <- res %>% filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    #group_by(Cell_type,Origin) %>%
    group_by(Cell_type) %>%
    mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))

####################################################################
# simple qq-plot  
####################################################################


fname=paste0(outFolder,"all.qqplot.png");
png(fname,width=800,height=800)
qq(res$pvalue)
dev.off()


####################################################################
# simple qq-plot / cell types with colors 
####################################################################

# qqplot to show the p-values splited by Origin and Location  
fname=paste0(outFolder,"split.qqplot.png");
p1 <- res2 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    #facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    facet_grid(~Location) +
    theme_bw()

ggsave(fname,p1,width=10,height=3)



#### some investigation
outFolder <- paste0("./11_results_insights/")
system(paste0("mkdir -p ",outFolder))

Locs<-unique(res$Location)
res<-res %>% filter(padj<0.1)
for (loc in Locs)
{
    res3<-res %>% filter(Location==loc)
    major_celltypes<-table(res3$Cell_type)
    major_celltypes<-major_celltypes[order(major_celltypes,decreasing=TRUE)]
    write.csv(major_celltypes,file=paste0(outFolder,"DEGs_major_celltypes_",loc,".csv"))
}


res2$logobserved<--log10(res2$pvalue)
res2$logexpected<--log10(res2$pexp)
Locs<-unique(res2$Location)
for (loc in Locs)
{
    print(loc)
    res3<-res2 %>% filter(Location==loc)

    res3<-res3 %>% filter( logobserved>logexpected )
    
    major_celltypes<-res3$logobserved
    names(major_celltypes)<-res3$Cell_type
    
    
    major_celltypes<-major_celltypes[order(major_celltypes,decreasing=TRUE)]
    print(unique(names(major_celltypes)[1:30]))
}
