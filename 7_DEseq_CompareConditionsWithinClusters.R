#############################################################
###  Calculating differentially expressed genes
### 
#############################################################


library(tidyverse)
##library(knitr)
library(DESeq2)
##library(annotables)
library(qqman)
library(tidyr)
library(Seurat)
library(reshape)
library(variancePartition)
library("BiocParallel")
 register(MulticoreParam(16))


# future::plan(strategy = 'multicore', workers = 16)
# options(future.globals.maxSize = 30 * 1024 ^ 3)

outFolder <- paste0("7_outputs_DESeq_ConditionsByCluster_res0.5/")
system(paste0("mkdir -p ",outFolder))

## Load cells
sc <- read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2021-06-28.rds")
md <- sc@meta.data
md<-md %>% mutate(Pregnancy_ID=str_match(Library,"[0-9]+"))

cSize <- md %>% dplyr::count(Location,seurat_clusters,Pregnancy_ID,Condition)

sSize <- cSize %>% filter(n>5) %>% dplyr::count(Location,seurat_clusters,Condition) %>% filter(n>1) #%>% filter(n>2)

kSize <- sSize %>% dplyr::count(Location,seurat_clusters) %>% filter(n>=2) %>%
    mutate(cname=paste(Location,seurat_clusters,sep="_"))

kSize

resList <- lapply(1:nrow(kSize),function(ii){
    ##ii <- 1
    cname <- kSize$cname[ii]
    cat("#",cname,"\n")
    ##
    md_i <- md %>%
        filter(
            Location==kSize$Location[ii],
            seurat_clusters==kSize$seurat_clusters[ii]) 
    
    ## samples (control/case)
    bti <- md_i %>% transmute(bti=paste(Condition,Pregnancy_ID,sep="_")) %>% unlist %>% factor
    
    ## data selected
    all_i <- sc@assays$RNA@data[,rownames(md_i)]
    ##
    
    ##
    X <- model.matrix( ~ 0 + bti)
    qr.X <- qr(X)
    qr.X$rank
    dim(X)
    YtX <- all_i %*% X
    YtX <- as.matrix(YtX)
    dim(YtX)
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    bti2 <- gsub("bti","",colnames(YtX))
    colnames(YtX) <- bti2
    ##
    cmat<-YtX
    ##
    ##
    # anno_i <- tibble(kbid=rownames(cmat),rs=rowSums(cmat),nz=rowSums(cmat>0)) %>%
    #     inner_join(anno %>% select(kbid,Chr,TSS,Strand,gene_name)) %>%
    #     filter(Chr %in% paste0("chr", 1:22), rs > 10, nz > 3) ## keep only autosomal
    # 
    
    anno_i <- tibble(kbid=rownames(cmat),rs=rowSums(cmat),nz=rowSums(cmat>0)) %>%
        filter( rs > 10, nz > 3) ## keep only autosomal
    ##
    ##
    #table(anno_i$Chr)
    ##
    genesel <- (unique(anno_i$kbid))
    cmat <- cmat[genesel,]
    dim(cmat)
    ##
    ##Create sample table
    cn<-colnames(cmat)
    x<-strsplit(cn,"_")
    ##
    cvt <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
    colnames(cvt)<-c("Group","Indiv")
    ##
    cvt$Group = relevel(factor(cvt$Group),"Control")
    ##
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~ Group)    
    dds <- DESeq(dds,parallel=TRUE)
    ##
    ## save DEseq object.    
    fname <- paste0(outFolder,cname,"_dds_",Sys.Date(),".rds")
    fname
    write_rds(dds,fname)
    ##
    ## Parse the results 
    ## ========================================================
    ##
    res <- results(dds)
    myres <- as.data.frame(res) %>%
        rownames_to_column("kbid") %>%
        left_join(anno_i)
    myres$cname <- cname
    ##
    write_tsv(myres,paste0(outFolder,cname,".",Sys.Date(),".txt"))
    ##
    nres<-nrow(myres %>% filter(padj<.1,abs(log2FoldChange)>0.5))
    cat("# Total sig for",cname,": ",nres,"\n")
    myres
})


res <- do.call(rbind,resList)
sum(res$padj<0.1,na.rm=TRUE)

Sys.Date()

system(paste0("mkdir -p ",outFolder))

res %>% filter(padj<0.1,abs(log2FoldChange)>0) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.",Sys.Date(),".tsv"))

res %>% write_tsv(paste0(outFolder,"ALL.combined.",Sys.Date(),".tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>0.0) %>%
    write_tsv(paste0(outFolder,"SIG.combined.2021-02-17.tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>1.0) %>% select(kbid,cname,log2FoldChange,padj) 

res %>% filter(padj<0.1,abs(log2FoldChange)>1) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.abslogfc_1",Sys.Date(),".tsv"))