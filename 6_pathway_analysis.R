###################################################
### pathway enrichment analysis

###################################################
library(tidyverse)
library(qqman)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(stringr)


outFolder <- paste0("6_pathway_enrichment/")
system(paste0("mkdir -p ",outFolder))

######################################################################################################


pathway_enrich<-function(res_gene=res3,cname_select,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
    cat ("=============== ",cname_select,"=============== ", "\n")
    pathway_enrich_cname_dir<-paste0(outFolder,cname_select,"/")
    system(paste0("mkdir -p ",pathway_enrich_cname_dir))
    result<-list()
    aux <- res_gene %>% filter(cname==cname_select)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    message(".................................")
    message("Number of DE genes: ",length(genes))
    #print(length(genes))
    
    message(".................................")
    message("enrichGO")
    ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Mm.eg.db,ont="BP",minGSSize=5)
    print(head(ego))
    result$enrichGO<-ego
    #save(ego,file=paste0(pathway_enrich_cname_dir,"ego.RData"))
    #write.csv(ego,file=paste0(pathway_enrich_cname_dir,"ego.csv"))
    
    print(".................................")
    print("enrichKEGG")
    ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="mmu",minGSSize=5)
    print(head(ekegg)) 
    result$enrichKEGG<-ekegg
    #save(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.RData"))
    #write.csv(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.csv"))
    
    
    message(".................................")
    message("enrichPathway")
    erpath <- enrichPathway(gene=genes,universe=geneUniv,organism="mouse",minGSSize=5)
    print(head(erpath))
    result$enrichPathway<-erpath
    #save(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.RData"))
    # write.csv(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.csv"))
    
    message(".................................")
    message("gseGO")
    # BP: biological_process, CC: cellular_component, MF: molecular_function
    gseGO.res <- gseGO(geneList,  OrgDb=org.Mm.eg.db,ont="BP",minGSSize=5)
    print(head(gseGO.res))
    result$gseGO<-gseGO.res
    #save(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.RData"))
    #write.csv(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.csv"))
    
    
    message(".................................")
    message("gsePathway")
    gseRPath.res <- gsePathway(geneList,organism="mouse",minGSSize=5)
    print(head(gseRPath.res))
    result$gsePathway<-gseRPath.res
    return (result)
}


# load DE genes
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")

# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cluster"),sep="_",remove=FALSE)
cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location)  

# Removing na pvalues
# Grouping pvalues based on the Location,Cell_type,and Origin
# Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 
# res2 <- res %>% filter(!is.na(pvalue)) %>%
#     arrange(pvalue) %>%
#     group_by(Cell_type) %>%
#     mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))


#ENTREZID id 
eg = bitr(res$kbid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

colnames(res)[1]<-"gene_name"

res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))




# # counting the number of DE genes per cname   
DE_per_cname<-sapply(unique(res3$cname), function(x,padj_cutoff=0.1,log2FoldChange_cutoff=0.5){
    aux <- res3 %>% filter(cname==x)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    length(genes)
    
})

cname_selected<-names(DE_per_cname)[which(DE_per_cname>5)]
result_pathway_en_list<-lapply(cname_selected, function(x) return(pathway_enrich(res3,x)))
names(result_pathway_en_list)<-cname_selected

save(result_pathway_en_list,file=paste0(outFolder,"pathwayEnrich_result.RData"))
which(DE_per_cname>0)


##########################################################################################  
#####                                  dot plot
##########################################################################################


load(paste0("6_pathway_enrichment/pathwayEnrich_result.RData"))
cname_selected<-names(result_pathway_en_list)

#gseGO
res_gseGO_list<-lapply(cname_selected, function(x)
{
    
    rs<-result_pathway_en_list[[x]]$gseGO@result %>% filter(qvalues<=0.05)
    dim1<-dim(rs)[1]
    
    if(min(dim1,5)>0)
    {
        res_en<-rs
        # to calculate GeneRatio=count/setSize
        #count
        gene_count<- res_en %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        
        ## merge with the original dataframe
        dot_df<- left_join(res_en, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
        
        dot_df<-dot_df[1:min(dim1,5),c("ID","Description" ,"enrichmentScore","p.adjust","GeneRatio")]
        dot_df$cname<-rep(x,min(dim1,5))
        dot_df
    }
})


res_df_gseGO <- do.call(rbind,res_gseGO_list)
res_df_gseGO<-res_df_gseGO %>% filter(p.adjust<0.1)

res_df_gseGO$Description
res_df_gseGO$cname
res_df_gseGO$p.adjust

mt<-matrix(nrow=length(unique(res_df_gseGO$Description)),ncol=length(unique(res_df_gseGO$cname)),0)
rownames(mt)<-unique(res_df_gseGO$Description)
colnames(mt)<-unique(res_df_gseGO$cname)

for ( i in unique(res_df_gseGO$Description))
{
    inx<-which(res_df_gseGO$Description==i)
    mt[i,res_df_gseGO$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_gseGO$orderpathways<-orderpathways[res_df_gseGO$Description]
res_df_gseGO$Location<-sapply(res_df_gseGO$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
})


res_df_gseGO <- res_df_gseGO %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()
res_df_gseGO$cname <- factor(res_df_gseGO$cname, levels = unique(res_df_gseGO$cname))

pdf(paste0(outFolder,"gseGO_cname_DotPlot.pdf"),width=20,height=15)
ggplot(res_df_gseGO, 
       aes(x = cname, y = reorder(Description,orderpathways))) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 14) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1)) + #text = element_text(size=30)
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
    ylab(NULL)+ 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+
    xlab(NULL) +
    theme(axis.text.y = element_text(hjust = 1))+
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
dev.off()


#enrichGO
res_enrichGO_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichGO)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichGO@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,5)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,5),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,5))
        }
        res_en }
}
)  



res_df_enrichGO <- do.call(rbind,res_enrichGO_list)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichGO<-res_df_enrichGO %>% filter(p.adjust<0.1) 

mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$cname)),0)
rownames(mt)<-unique(res_df_enrichGO$Description)
colnames(mt)<-unique(res_df_enrichGO$cname)


for ( i in unique(res_df_enrichGO$Description))
{
    inx<-which(res_df_enrichGO$Description==i)
    mt[i,res_df_enrichGO$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]


res_df_enrichGO$Location<-sapply(res_df_enrichGO$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
})


res_df_enrichGO <- res_df_enrichGO %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()

res_df_enrichGO$cname <- factor(res_df_enrichGO$cname, levels = unique(res_df_enrichGO$cname))

pdf(paste0(outFolder,"enrichGO_cname_DotPlot.pdf"),width=29,height=19)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = reorder(Description,orderpathways))) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    #theme_black()+
    theme_bw()+
    theme(axis.text.y = element_text(hjust = 1))+
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
dev.off()


####################

#enrichKEGG
res_enrichKEGG_list<-lapply(cname_selected, function(x)
{
    if(length(result_pathway_en_list[[x]]$enrichKEGG)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichKEGG@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,5)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,5),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,5))
        }
        res_en }
}
)  

res_df_enrichKEGG <- do.call(rbind,res_enrichKEGG_list)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})


res_df_enrichKEGG<-res_df_enrichKEGG %>% filter(p.adjust<0.1) 

mt<-matrix(nrow=length(unique(res_df_enrichKEGG$Description)),ncol=length(unique(res_df_enrichKEGG$cname)),0)
rownames(mt)<-unique(res_df_enrichKEGG$Description)
colnames(mt)<-unique(res_df_enrichKEGG$cname)


for ( i in unique(res_df_enrichKEGG$Description))
{
    inx<-which(res_df_enrichKEGG$Description==i)
    mt[i,res_df_enrichKEGG$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichKEGG$orderpathways<-orderpathways[res_df_enrichKEGG$Description]

res_df_enrichKEGG$Location<-sapply(res_df_enrichKEGG$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
})




res_df_enrichKEGG$Location[which(res_df_enrichKEGG$Location=="Myometrium")]<-"a"
res_df_enrichKEGG$Location[which(res_df_enrichKEGG$Location=="Decidua")]<-"b"
res_df_enrichKEGG$Location[which(res_df_enrichKEGG$Location=="Cervix")]<-"c"
res_df_enrichKEGG<-res_df_enrichKEGG %>% arrange(desc(Location))%>%  dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>% ungroup()
res_df_enrichKEGG$Location[which(res_df_enrichKEGG$Location=="a")]<-"Uterus"
res_df_enrichKEGG$Location[which(res_df_enrichKEGG$Location=="b")]<-"Decidua"
res_df_enrichKEGG$Location[which(res_df_enrichKEGG$Location=="c")]<-"Cervix"

res_df_enrichKEGG$cname <- factor(res_df_enrichKEGG$cname, levels = unique(res_df_enrichKEGG$cname))
res_df_enrichKEGG$Location <- factor(res_df_enrichKEGG$Location, levels = unique(res_df_enrichKEGG$Location))


pdf(paste0(outFolder,"enrichKEGG_cname_DotPlot.pdf"),width=20,height=12)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = reorder(Description,orderpathways))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    #theme_black()+
    theme_bw()+
    #facet_grid(.~Location)+
    #theme(axis.text.x = element_text(angle = 45))+
    theme(axis.text.y = element_text(hjust = 1))+
    #theme(text = element_text(size=40)) +
    #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
dev.off()


####################

#enrichPathway
res_enrichPathway_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichPathway)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichPathway@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,5)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,5),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,5))
        }
        res_en }
}
)  

res_df_enrichPathway <- do.call(rbind,res_enrichPathway_list)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})
res_df_enrichPathway<-res_df_enrichPathway %>% filter(p.adjust<0.1) 



mt<-matrix(nrow=length(unique(res_df_enrichPathway$Description)),ncol=length(unique(res_df_enrichPathway$cname)),0)
rownames(mt)<-unique(res_df_enrichPathway$Description)
colnames(mt)<-unique(res_df_enrichPathway$cname)


for ( i in unique(res_df_enrichPathway$Description))
{
    inx<-which(res_df_enrichPathway$Description==i)
    mt[i,res_df_enrichPathway$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichPathway$orderpathways<-orderpathways[res_df_enrichPathway$Description]



res_df_enrichPathway$Location<-sapply(res_df_enrichPathway$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
})

#res_df_enrichPathway<-res_df_enrichPathway[1:15,]

res_df_enrichPathway <- res_df_enrichPathway %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()

res_df_enrichPathway$cname <- factor(res_df_enrichPathway$cname, levels = unique(res_df_enrichPathway$cname))


pdf(paste0(outFolder,"enrichPathway_cname_DotPlot.pdf"),width=29,height=15)
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest 
       aes(x = cname, y = reorder(Description,orderpathways))) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(axis.text.y = element_text(hjust = 1))+
    theme(axis.text.x = element_text(angle = 90,  hjust=1)) 
dev.off()

############################################################

# Pathway analysis 
### DE genes combined

############################################################


# load DE genes
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res <- res %>% separate(cname,c("Location","Cluster"),sep="_",remove=FALSE)

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type 
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location)  

#ENTREZID id 
eg = bitr(res$kbid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

colnames(res)[1]<-"gene_name"

res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

genes <- filter(res3,padj<0.1) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
geneUniv <- res3 %>% dplyr::select(ENTREZID) %>% unlist %>% unique

geneList <- -log10(res3$pvalue)
names(geneList) <- res3$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

message(".................................")
message("enrichGO")
ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb="org.Mm.eg.db",ont="BP",minGSSize = 5)
print(head(ego))

res_df_enrichGO<-ego@result %>% filter(qvalue<=0.1)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichGO<-res_df_enrichGO %>% filter(p.adjust<0.1) 
res_df_enrichGO<-res_df_enrichGO[1:15,]


pdf(paste0(outFolder,"enrichGO_combined_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = reorder(Description,GeneRatio))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    #theme_black()+
    theme_bw()+
    theme(axis.text.y = element_text(hjust = 1))+
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
dev.off()




print(".................................")
print("enrichKEGG")
ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="mouse",minGSSize = 5)
print(head(ekegg)) 

res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.05)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichKEGG<-res_df_enrichKEGG %>% filter(p.adjust<0.1) 

res_df_enrichKEGG<-res_df_enrichKEGG[1:15,]

pdf(paste0(outFolder,"enrichKEGG.combined_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = reorder(Description,GeneRatio))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    #theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    #theme_black()+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 45))+
    theme(axis.text.y = element_text(hjust = 1))+
    #theme(text = element_text(size=30)) +
    #theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
dev.off()

message(".................................")
message("enrichPathway")
erpath <- enrichPathway(gene=genes,universe=geneUniv,organism="mouse",minGSSize = 5)
print(head(erpath))

res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.05)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichPathway<-res_df_enrichPathway[1:15,]
pdf(paste0(outFolder,"enrichPathway.combined_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = reorder(Description,GeneRatio))) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    #theme_black()+
    theme_bw()+
    #theme(axis.text.x = element_text(angle = 45))+
    theme(axis.text.y = element_text(hjust = 1))+
    #theme(text = element_text(size=30))+
    #theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
    theme(axis.text.x = element_text(angle = 90, hjust=1)) 
dev.off()




###########################################

# Boxplot 
###########################################

load(paste0("6_pathway_enrichment/pathwayEnrich_result.RData"))
cname_selected<-names(result_pathway_en_list)

cell.type.annotation<-read.delim("cell-type.annotation.txt")
clust2Names<-cell.type.annotation$Cell.type
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster


cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F","#FAFA33","#FFA6C9","#F4C2C2","#1034A6","#08E8DE","#00BFFF","#6F00FF")
names(cluster.Colors)<-clust2Names #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")


load(paste0("6_pathway_enrichment/pathwayEnrich_result.RData"))
cname_selected<-names(result_pathway_en_list)

cname_selected_number<-sapply(cname_selected,function(x){
    y=unlist(strsplit(x,"_"))
    return(y[1])}
)

cname_selected_location<-sapply(cname_selected,function(x){
    y=unlist(strsplit(x,"_"))
    return(y[length(y)])}
)

names(result_pathway_en_list)<-paste0(clust2Names[cname_selected_number],"_",cname_selected_location)
cname_selected<-names(result_pathway_en_list)


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res <- res %>% separate(cname,c("Location","Cluster"),sep="_",remove=FALSE)

cell.type.annotation<-read.delim("cell-type.annotation.txt")

clust2Names<-cell.type.annotation$Cell.type #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location)  

#ENTREZID id 
eg = bitr(res$kbid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

colnames(res)[1]<-"gene_name"

res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

#enrichGO
res_enrichGO_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichGO)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichGO@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,10)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust","geneID")]
            res_en$cname<-rep(x,min(dim1,10))
        }
        res_en }
}
)  

res_df_enrichGO <- do.call(rbind,res_enrichGO_list)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichGO<-res_df_enrichGO %>% filter(p.adjust<0.1) 
#res_df_enrichGO<-res_df_enrichGO[1:15,]


mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$cname)),0)
rownames(mt)<-unique(res_df_enrichGO$Description)
colnames(mt)<-unique(res_df_enrichGO$cname)


for ( i in unique(res_df_enrichGO$Description))
{
    inx<-which(res_df_enrichGO$Description==i)
    mt[i,res_df_enrichGO$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]


res_df_enrichGO$Location<-sapply(res_df_enrichGO$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
})


res_df_enrichGO <- res_df_enrichGO %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()

res_df_enrichGO$cname <- factor(res_df_enrichGO$cname, levels = unique(res_df_enrichGO$cname))


###  pathway: defense response to bacterium
# Cell type: 5_Epithelial−1 (Basal)_Cervix



select_path<-rownames(mt)[rowSums(mt)>=2]



outFolder <- paste0("6_pathway_enrichment/Boxplots/")
system(paste0("mkdir -p ",outFolder))

subFolder<-"GO"
system(paste0("mkdir -p ",outFolder,subFolder,"/"))

for (i in 1:length(select_path))
{
    path<-select_path[i]
    
    #path<-"defense response to bacterium"
    #path<-"regulation of smooth muscle cell migration"
    print(path)
    category<-names(which(mt[path,]==1))
    
    
        index<-which(res_df_enrichGO$Description==path)
        geneList <- unique(unlist(str_split(res_df_enrichGO$geneID[index],"/")))
        geneList<-e2g[geneList]
        
        
        res_boxplot<-res3 %>% filter(gene_name %in% geneList)
        
        exist_celltype<-names(which(table(res_boxplot$cname)==length(geneList)))
        
        print(exist_celltype)
        
        if (length(exist_celltype)>1)
        {
            res_boxplot<-res_boxplot %>%filter (cname %in% exist_celltype)
            
            
            res_boxplot<-res_boxplot[order(as.numeric(res_boxplot$Cluster),decreasing = FALSE),]
            
            res_boxplot$Cluster <- factor(res_boxplot$Cluster,levels=unique(res_boxplot$Cluster))
            res_boxplot$cname <- factor(res_boxplot$cname,levels=unique(res_boxplot$cname))
            
            
            
            
            #pdf(paste0(outFolder,subFolder,"/",gsub(" ", "_", path),"_",category,"_enGO.pdf"), width=30, height=10)
            pdf(paste0(outFolder,subFolder,"/",gsub(" ", "_", path),"_enGO.pdf"), width=10, height=8,onefile=TRUE)
            #paste0(outFolder,gsub(" ", "_", path),"_",category,"_enGO.pdf")
            ggplot(res_boxplot, aes(x=cname, y=log2FoldChange, fill=Cell_type))+
                geom_boxplot(outlier.shape=NA)+
                ylab("log2FoldChange")+
                xlab("")+
                scale_fill_manual(values=cluster.Colors) +
                theme_bw()+
                theme(legend.position="none",axis.text=element_text(size=20),axis.text.x = element_text(vjust = 0.2,angle = 90, hjust=1))
            dev.off()
        }
       
   
}


### Reactome 

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res <- res %>% separate(cname,c("Location","Cluster"),sep="_",remove=FALSE)

cell.type.annotation<-read.delim("cell-type.annotation.txt")

clust2Names<-cell.type.annotation$Cell.type #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,"_",clust2Names)
names(clust2Names)<-cell.type.annotation$Cluster #c(0:23)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Location)  

#ENTREZID id 
eg = bitr(res$kbid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

colnames(res)[1]<-"gene_name"

res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

##################################################
#enrichPathway

res_enrichPathway_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichPathway)>0)
    {
        print(x)
        rs<-result_pathway_en_list[[x]]$enrichPathway@result %>% filter(qvalue<0.1)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,20)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,20),c("ID","Description" ,"GeneRatio","p.adjust","geneID")]
            res_en$cname<-rep(x,min(dim1,20))
        }
        res_en }
}
)  

res_df_enrichPathway <- do.call(rbind,res_enrichPathway_list)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})
res_df_enrichPathway<-res_df_enrichPathway %>% filter(p.adjust<0.1) 



mt<-matrix(nrow=length(unique(res_df_enrichPathway$Description)),ncol=length(unique(res_df_enrichPathway$cname)),0)
rownames(mt)<-unique(res_df_enrichPathway$Description)
colnames(mt)<-unique(res_df_enrichPathway$cname)


for ( i in unique(res_df_enrichPathway$Description))
{
    inx<-which(res_df_enrichPathway$Description==i)
    mt[i,res_df_enrichPathway$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichPathway$orderpathways<-orderpathways[res_df_enrichPathway$Description]



res_df_enrichPathway$Location<-sapply(res_df_enrichPathway$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
})


res_df_enrichPathway <- res_df_enrichPathway %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()

res_df_enrichPathway$cname <- factor(res_df_enrichPathway$cname, levels = unique(res_df_enrichPathway$cname))






select_path<-rownames(mt)[rowSums(mt)>=2]

outFolder <- paste0("6_pathway_enrichment/Boxplots/")
system(paste0("mkdir -p ",outFolder))

subFolder<-"Reactome"
system(paste0("mkdir -p ",outFolder,subFolder,"/"))

for (i in 1:length(select_path))
{
    path<-select_path[i]
    
    #path<-"defense response to bacterium"
    #path<-"regulation of smooth muscle cell migration"
    print(path)
    category<-names(which(mt[path,]==1))
    
    
    index<-which(res_df_enrichPathway$Description==path)
    geneList <- unique(unlist(str_split(res_df_enrichPathway$geneID[index],"/")))
    geneList<-e2g[geneList]
    
    
    res_boxplot<-res3 %>% filter(gene_name %in% geneList)
    
    exist_celltype<-names(which(table(res_boxplot$cname)==length(geneList)))
    
    print(exist_celltype)
    
    if (length(exist_celltype)>1)
    {
        res_boxplot<-res_boxplot %>%filter (cname %in% exist_celltype)
        
        
        res_boxplot<-res_boxplot[order(as.numeric(res_boxplot$Cluster),decreasing = FALSE),]
        
        res_boxplot$Cluster <- factor(res_boxplot$Cluster,levels=unique(res_boxplot$Cluster))
        res_boxplot$cname <- factor(res_boxplot$cname,levels=unique(res_boxplot$cname))
        
        
        
        
        #pdf(paste0(outFolder,subFolder,"/",gsub(" ", "_", path),"_",category,"_enGO.pdf"), width=30, height=10)
        pdf(paste0(outFolder,subFolder,"/",gsub(" ", "_", path),"_enReactome.pdf"), width=20, height=10,onefile=TRUE)
        #paste0(outFolder,gsub(" ", "_", path),"_",category,"_enGO.pdf")
        ggplot(res_boxplot, aes(x=cname, y=log2FoldChange, fill=Cell_type))+
            geom_boxplot(outlier.shape=NA)+
            ylab("log2FoldChange")+
            xlab("")+
            scale_fill_manual(values=cluster.Colors) +
            theme_bw()+
            theme(legend.position="none",axis.text=element_text(size=20),axis.text.x = element_text(vjust = 0.2,angle = 90, hjust=1))
        dev.off()
    }
    
    
}



#####################
# KEGG
####################


res_enrichKEGG_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichKEGG)>0)
    {
        print(x)
        rs<-result_pathway_en_list[[x]]$enrichKEGG@result %>% filter(qvalue<0.1)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,20)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,20),c("ID","Description" ,"GeneRatio","p.adjust","geneID")]
            res_en$cname<-rep(x,min(dim1,20))
        }
        res_en }
}
)  

res_df_enrichKEGG <- do.call(rbind,res_enrichKEGG_list)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})
res_df_enrichKEGG<-res_df_enrichKEGG %>% filter(p.adjust<0.1) 



mt<-matrix(nrow=length(unique(res_df_enrichKEGG$Description)),ncol=length(unique(res_df_enrichKEGG$cname)),0)
rownames(mt)<-unique(res_df_enrichKEGG$Description)
colnames(mt)<-unique(res_df_enrichKEGG$cname)


for ( i in unique(res_df_enrichKEGG$Description))
{
    inx<-which(res_df_enrichKEGG$Description==i)
    mt[i,res_df_enrichKEGG$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichKEGG$orderpathways<-orderpathways[res_df_enrichKEGG$Description]



res_df_enrichKEGG$Location<-sapply(res_df_enrichKEGG$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
})

res_df_enrichKEGG <- res_df_enrichKEGG %>% dplyr::group_by(Location) %>% arrange(desc(orderpathways),.by_group = TRUE) %>%ungroup()

res_df_enrichKEGG$cname <- factor(res_df_enrichKEGG$cname, levels = unique(res_df_enrichKEGG$cname))






select_path<-rownames(mt)[rowSums(mt)>=2]

select_path<-"NF-kappa B signaling pathway" 
select_path<-"IL-17 signaling pathway"   
select_path<-"TNF signaling pathway"      
#select_path<-orderpathways[which(orderpathways==1)]
#select_cname<-colSums(mt)[which(colSums(mt)==1)]

outFolder <- paste0("6_pathway_enrichment/Boxplots/")
system(paste0("mkdir -p ",outFolder))

subFolder<-"KEGG"
system(paste0("mkdir -p ",outFolder,subFolder,"/"))

for (i in 1:length(select_path))
{
    path<-select_path[i]
    
    #path<-"defense response to bacterium"
    #path<-"regulation of smooth muscle cell migration"
    print(path)
    category<-names(which(mt[path,]==1))
    
    
    index<-which(res_df_enrichKEGG$Description==path)
    
    allgenes<-lapply(1:length(res_df_enrichKEGG$geneID),function(x){
        
        if (res_df_enrichKEGG$Description[x]==path) return(unique(unlist(str_split(res_df_enrichKEGG$geneID[x],"/"))))
        })
    
    geneList<-unique(unlist(allgenes))
    
    #geneList <- unique(unlist(str_split(res_df_enrichKEGG$geneID[index],"/")))
    
    
    geneList<-e2g[geneList]
    
    
    res_boxplot<-res3 %>% filter(gene_name %in% geneList)
    
    
    #res_boxplot<-res_boxplot %>% filter(Cell_type=="6_Neutrophil")
    #exist_celltype<-names(which(table(res_boxplot$cname)==length(geneList)))
    exist_celltype<-names(which(table(res_boxplot$cname)>=55))
    
    print(exist_celltype)
    
    if (length(exist_celltype)>1)
    {
        #res_boxplot<-res_boxplot %>%filter (cname %in% exist_celltype)
        
        
        res_boxplot<-res_boxplot[order(as.numeric(res_boxplot$Cluster),decreasing = FALSE),]
        
        res_boxplot$Cluster <- factor(res_boxplot$Cluster,levels=unique(res_boxplot$Cluster))
        res_boxplot$cname <- factor(res_boxplot$cname,levels=unique(res_boxplot$cname))
        
        
        
        
        #pdf(paste0(outFolder,subFolder,"/",gsub(" ", "_", path),"_",category,"_enGO.pdf"), width=30, height=10)
        #pdf(paste0(outFolder,subFolder,"/","6_Neutrophil_",gsub(" ", "_", path),"_KEGG.pdf"), width=6, height=8,onefile=TRUE)
        pdf(paste0(outFolder,subFolder,"/",gsub(" ", "_", path),"_KEGG.pdf"), width=20, height=10,onefile=TRUE)
        
        #paste0(outFolder,gsub(" ", "_", path),"_",category,"_enGO.pdf")
        ggplot(res_boxplot, aes(x=cname, y=log2FoldChange, fill=Cell_type))+
            
            geom_boxplot(outlier.shape=NA)+
            ylab("log2FoldChange")+
            xlab("")+
            scale_fill_manual(values=cluster.Colors) +
            theme_bw()+
            
            theme(legend.position="none",axis.text=element_text(size=20),axis.text.x = element_text(vjust = 0.2,angle = 90, hjust=1))
            #theme_update(plot.title = element_text(hjust = 0.5))+
            #ggtitle(path)
        dev.off()
    }
    
    
}




path<-"Innate Immune System"
#path<-"regulation of smooth muscle cell migration"
print(path)
category<-names(which(mt[path,]==1))

index<-which(res_df_enrichPathway$Description==path)
geneList <- unique(unlist(str_split(res_df_enrichPathway$geneID[index],"/")))
geneList<-e2g[geneList]


res_boxplot<-res3 %>% filter(gene_name %in% geneList)

exist_celltype<-names(which(table(res_boxplot$cname)==length(geneList)))
res_boxplot<-res_boxplot %>%filter (cname %in% exist_celltype)


res_boxplot<-res_boxplot[order(as.numeric(res_boxplot$Cluster),decreasing = FALSE),]

res_boxplot$Cluster <- factor(res_boxplot$Cluster,levels=unique(res_boxplot$Cluster))
res_boxplot$cname <- factor(res_boxplot$cname,levels=unique(res_boxplot$cname))



pdf(paste0(outFolder,gsub(" ", "_", path),"_",category,"_enPathway.pdf"), width=30, height=10)
paste0(outFolder,gsub(" ", "_", path),"_",category,"_enGO.pdf")
ggplot(res_boxplot, aes(x=cname, y=log2FoldChange, fill=Cell_type))+
    geom_boxplot(outlier.shape=NA)+
    ylab("log2FoldChange")+
    xlab("")+
    theme_bw()+
    scale_fill_manual(values=cluster.Colors) +
    theme(legend.position="none",axis.text=element_text(size=20),axis.text.x = element_text(vjust = 0.2,angle = 90, hjust=1)) 
dev.off()



##########################################################################################
# Myometrium specific
#result_pathway_en_list
##########################################################################################

load(paste0("./6_pathway_enrichment/Myometrium/pathwayEnrich_result.RData"))

#enrichGO
cname_selected<-names(result_pathway_en_list)
res_enrichGO_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichGO)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichGO@result %>% filter(qvalue<=0.1)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,20)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,20),c("ID","Description" ,"GeneRatio","p.adjust","geneID")]
            res_en$cname<-rep(x,min(dim1,20))
        }
        res_en }
}
)  

res_df_enrichGO <- do.call(rbind,res_enrichGO_list)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichGO<-res_df_enrichGO %>% filter(p.adjust<0.1) 
#res_df_enrichGO<-res_df_enrichGO[1:15,]


mt<-matrix(nrow=length(unique(res_df_enrichGO$Description)),ncol=length(unique(res_df_enrichGO$cname)),0)
rownames(mt)<-unique(res_df_enrichGO$Description)
colnames(mt)<-unique(res_df_enrichGO$cname)


for ( i in unique(res_df_enrichGO$Description))
{
    inx<-which(res_df_enrichGO$Description==i)
    mt[i,res_df_enrichGO$cname[inx]]<-1
}

orderpathways<-rowSums(mt)
orderpathways<-orderpathways[order(orderpathways,decreasing = TRUE)]

res_df_enrichGO$orderpathways<-orderpathways[res_df_enrichGO$Description]


res_df_enrichGO$Location<-sapply(res_df_enrichGO$cname, function(x){
    x<-unlist(strsplit(x,"_"))
    return(x[length(x)])
    
})


path<-"regulation of smooth muscle cell migration"
print(path)
category<-names(which(mt[path,]==1))

index<-which(res_df_enrichGO$Description==path)
geneList <- unique(unlist(str_split(res_df_enrichGO$geneID[index],"/")))
geneList<-e2g[geneList]


res_boxplot<-res3 %>% filter(gene_name %in% geneList)
res_boxplot<-res_boxplot[order(as.numeric(res_boxplot$Cluster),decreasing = FALSE),]

res_boxplot$Cluster <- factor(res_boxplot$Cluster,levels=unique(res_boxplot$Cluster))
res_boxplot$cname <- factor(res_boxplot$cname,levels=unique(res_boxplot$cname))
paste0(outFolder,gsub(" ", "_", path),"_",category,"_enGO.pdf")
exist_celltype<-names(which(table(res_boxplot$cname)==length(geneList)))
res_boxplot<-res_boxplot %>%filter (cname %in% exist_celltype)
pdf(paste0(outFolder,gsub(" ", "_", path),"_",category,"_enGO.pdf"), width=30, height=10)
ggplot(res_boxplot, aes(x=cname, y=log2FoldChange, fill=Cell_type))+
    geom_boxplot(outlier.shape=NA)+
    ylab("log2FoldChange")+
    xlab("")+
    theme_bw()+
    scale_fill_manual(values=cluster.Colors) +
    theme(legend.position="none",axis.text=element_text(size=20),axis.text.x = element_text(vjust = 0.2,angle = 90, hjust=1)) 
dev.off()
