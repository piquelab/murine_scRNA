library(tidyverse)
library(dplyr)
library(qqman)
library(org.Hs.eg.db)
library(clusterProfiler)
library(reshape2)
library(ggplot2)

########################################################
# comparison with bulk data

########################################################
outFolder="12_comparison_with_bulk/"
system(paste0("mkdir -p ", outFolder))


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


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_res0.5/ALL.combined.2021-06-29.tsv")
res <- res %>% separate(cname,c("Location","Cell_type"),sep="_",remove=FALSE)
res$Cell_type<-clust2Names[res$Cell_type]
celltype_DE<-table(res$Cell_type,res$Location)

eg = bitr(res$kbid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
colnames(res)[1]<-"gene_name"
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))





# bulk data 
datasets<-c("CERVIX_LPS_PregLPS-PregCtl_ALLList.txt","UTERUS_LPS_PregLPS-PregCtl_ALLList.txt")
names(datasets)<-c("Cervix","Uterus")

#threshold<-"Bulk"
threshold<-"single"
#dirfilter<-"Bulk_padj0.1/"
if(threshold=="single") dirfilter<-"Single_padj0.1/"
for (i in 1:length(datasets))
{
  
  bulk<-names(datasets)[i]
  system(paste0("mkdir -p ", outFolder,dirfilter,"Bulk_",bulk,"/"))
  dataset<-datasets[i]
  ref_data<-read_tsv(dataset)
  ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
  colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
  ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
  ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
  ref_data<-ref_data %>%filter(!is.na(ENTREZID) & !is.na(R.Log2FC) & !is.na(Rpvalue))
  
  cor_locations<-lapply(unique(res$Location), function(loc)
    {
   
    subFolder=loc
    system(paste0("mkdir -p ", outFolder,dirfilter,"/","Bulk_",bulk,"/",subFolder,"/"))
    finalFolder<-paste0(outFolder,dirfilter,"Bulk_",bulk,"/",subFolder,"/")
    
    cor_celltypes_list<-lapply (unique(res$Cell_type), function(cl ){
      
      res_filter<-res %>% filter(Cell_type ==cl & Location==loc)
      resJoin<-res_filter %>% inner_join(ref_data)
      #resJoin<-resJoin %>% filter(padj< 0.1 )
      #cat(table(resJoin$ref_color),sep=":")
      resJoin$t<-resJoin$log2FoldChange/resJoin$lfcSE
      
      
      
      resJoin_singlecell<-resJoin %>% filter(padj<0.1)
      resJoin_bulk<-resJoin %>% filter(Rpadj<0.1)
      
      #res_filter2<-res %>% filter(padj<0.1)
      #res_filter2<-res_filter2 %>% filter(Cell_type %in% cl)
      # resJoin2<-res_filter2 %>% inner_join(ref_data)
      # #resJoin<-resJoin %>% filter(padj< 0.1 )
      # #cat(table(resJoin$ref_color),sep=":")
      # resJoin2$t<-resJoin2$log2FoldChange/resJoin2$lfcSE
      
      gtitle<-paste0(cl,"-",tolower(loc), " vs Bulk ",tolower(bulk))
      
      spearman_all<-c()
      spearman_cor_pval<-c()
      
      if(nrow(resJoin_bulk)>=10 | nrow(resJoin_singlecell)>=10)
      {
        
        
        if(nrow(resJoin_bulk)>=10 & threshold!="single")
        {
          spearman_all<-cor.test(resJoin_bulk$t,resJoin_bulk$Rt,method="spearman",na.rm=TRUE)
        }
          
          
        
        if (nrow(resJoin_singlecell)>=10 & threshold=="single" )
        {
          spearman_all<-cor.test(resJoin_singlecell$t,resJoin_singlecell$Rt,method="spearman",na.rm=TRUE)
        }
          
        if (length(spearman_all)>0)
        {
          spearman_cor_pval<-c(spearman_all$estimate,spearman_all$p.value)
          names(spearman_cor_pval)<-c("spearman_cor","pvalue")
          #gtitle<-paste0(cl,"-",tolower(loc), " vs bulk ",tolower(bulk)," (cor= ",round(spearman_cor_pval[1],2),", p= ",formatC(spearman_cor_pval[2], format = "e", digits = 2),")")
          
        }
          
        }
      
      res_rest<-res %>% filter(padj<0.1)
     resJoin$ref_color <- "None"          
      #light blue
      resJoin$ref_color[resJoin$padj<0.1  &resJoin$Rpadj>0.1 ] <- "Only single cell"  
      #purpule  //res_rest union of all DEGs across all cell types
      resJoin$ref_color[(!resJoin$ENTREZID %in% unique(res_rest$ENTREZID)) &  resJoin$Rpadj<0.1  ] <- "Only bulk" 
      #blue
      resJoin$ref_color[resJoin$padj<0.1 & resJoin$Rpadj<0.1]="Single cell and bulk" 
      
      # // resJoin$ref_color<-as.fac
      resJoin$ref_color <- factor(resJoin$ref_color,levels=unique(resJoin$ref_color))
      
      
      table(resJoin$ref_color)
      
      cat(table(resJoin$ref_color),sep=":")
      #resJoin$t<-resJoin$log2FoldChange/resJoin$lfcSE
      if (nrow(resJoin)>0)
      {
        p2 <- resJoin %>% arrange(ref_color) %>%
          ggplot(aes(Rt,t,color=ref_color)) +
          geom_point()+ #aes(colour = ref_color)) +
          #geom_smooth(method=lm, se=FALSE,linetype = "dashed", color="black")+  # filtered way
          xlab(paste0(bulk," (bulk)"," standardized Log2FC"))+
          ylab("Standardized log2FC")+
          scale_color_manual(name="Differentially expressed",values=c("None"="#CCCCCC","Only single cell"="#BDD7EE","Only bulk"="#DD99DD","Single cell and bulk"="#0000EE"))+
          ggtitle(gtitle)+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5))
        
        cl<-gsub("\\(", "",cl)
        cl<-gsub("\\)", "",cl)
        cl<-gsub("\\ ", "-",cl)
        
        
        fname=paste0(finalFolder,paste0(cl,".sc.",loc,"_bulk.",bulk,".png"))
        ggsave(fname,p2,width=9,height=7)
        
      }
        
      if (length(spearman_cor_pval)==0)
        spearman_cor_pval<-c(NA,NA)
      return(spearman_cor_pval)
      
    }) 
    
    cor_celltypes<-do.call(rbind,cor_celltypes_list)
    
    #cor_celltypes<-t(cor_celltypes)
    colnames(cor_celltypes)<-c("spearman_cor","spearman_pvalue")
    #cor_celltypes<-cor_celltypes %>% arrange(desc(spearman_cor,-log(spearman_pvalue)))
    
    cor_celltypes<-as.data.frame(cor_celltypes)
    
    cor_celltypes$celltype<-unique(res$Cell_type)
    cor_celltypes$location<-loc
    
    cor_celltypes<-cor_celltypes %>% filter(!is.na(spearman_cor))
    cor_celltypes<-cor_celltypes %>% arrange(desc(spearman_cor))
    cor_celltypes$spearman_padj<-p.adjust(cor_celltypes$spearman_pvalue,"fdr")
    cor_celltypes2<-cor_celltypes %>% select(celltype,spearman_cor,spearman_pvalue,spearman_padj)
    
    fname=paste0(finalFolder,paste0("cor.sc.",loc,"_bulk.",bulk,".csv"))
    write.csv(cor_celltypes2,file=fname)
    
    
    cor_celltypes$sig<-sapply(cor_celltypes$spearman_padj, function(x){
      if (x<=0.0001) x="****"
      else if (x<=0.001)x="***" 
      else if (x<=0.01) x="**"
      else if (x<0.1) x="*"
      else if (x>=0.1) return ("ns")
    })
     
   
    
    cor_celltypes<- cor_celltypes%>% filter(spearman_cor!=0)
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".v.pdf"))
    pdf(fname,width=10,height=4.5)
    p2<-ggplot(data=cor_celltypes, aes(x=celltype, y=spearman_cor,fill=celltype)) +
      geom_bar(stat="identity",position="stack")+
      #geom_text(aes(label=sig,vjust = -sign(spearman_cor)), vjust=1.6, color="black", size=3.5)+
      theme_bw()+
      scale_fill_manual(values=cluster.Colors) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))+
      theme(legend.position="none")+
      xlab("")+
      ylab("Spearman correlation")
    p2
    dev.off() 
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".v.png"))
    ggsave(fname,p2,width=10,height=4)
    
    
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".v.pdf"))
    pdf(fname,width=10,height=4)
    p2<-ggplot(data=cor_celltypes, aes(x=celltype, y=spearman_cor,fill=celltype)) +
      geom_bar(stat="identity",position="stack")+
      #geom_text(aes(label=sig,vjust = -sign(spearman_cor)), vjust=1.6, color="black", size=3.5)+
      theme_bw()+
      scale_fill_manual(values=cluster.Colors) +
      theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="none")+
      #theme(legend.position="none")+
      xlab("")+
      ylab("Spearman correlation")
    p2
    dev.off()
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".withstars.v.png"))
    ggsave(fname,p2,width=10,height=4)
    
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".h.pdf"))
    pdf(fname,width=12,height=6)
    p2<-ggplot(cor_celltypes, aes(x=reorder(celltype,-spearman_cor), y=spearman_cor,fill=celltype))+#,color=mycolor)+#,fill=DE) +
      geom_bar(stat='identity') +
      geom_text(aes(label=sig,hjust = -sign(spearman_cor)), vjust=0.7, color="black", size=3.5)+
      theme_bw()+
      scale_fill_manual(values=cluster.Colors) +
      #facet_grid(.~Location,scales="free") + 
      coord_flip() +
      theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="none")+
      #theme(legend.title=element_blank())+
      ylab("Spearman correlation")+
      xlab("")
    p2
    dev.off() 
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".h.png"))
    ggsave(fname,p2,width=8,height=10)
    
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".withstars.h.pdf"))
    pdf(fname,width=8,height=10)
    p2<-ggplot(cor_celltypes, aes(x=reorder(celltype,-spearman_cor), y=spearman_cor,fill=celltype))+#,color=mycolor)+#,fill=DE) +
      geom_bar(stat='identity') +
      geom_text(aes(label=sig,hjust = -sign(spearman_cor)), vjust=0.7, color="black", size=3.5)+
      theme_bw()+
      theme(legend.position="none")+
      scale_fill_manual(values=cluster.Colors) +
      #facet_grid(.~Location,scales="free") + 
      coord_flip() +
      #theme(legend.title=element_blank())+
      ylab("Spearman correlation")+
      xlab("")
    p2
    dev.off() 
    
    fname=paste0(finalFolder,paste0("barplot.sc.",loc,"_bulk.",bulk,".withstars.h.png"))
    ggsave(fname,p2,width=8,height=10)
    
    
    return(cor_celltypes)
    
  })
  
  
  total_cor_locations<-do.call(rbind,cor_locations)
  finalFolder<-paste0(outFolder,dirfilter,"/","Bulk_",bulk,"/")
  fname=paste0(finalFolder,paste0("cor_all_with_bulk.",bulk,".csv"))
  #total_cor_locations<-total_cor_locations[,-1]
  total_cor_locations$bulk<-bulk
  #total_cor_locations<-total_cor_locations %>% arrange(desc(spearman_cor,-log(spearman_pvalue)))
  write.csv(total_cor_locations,file=fname)
}






