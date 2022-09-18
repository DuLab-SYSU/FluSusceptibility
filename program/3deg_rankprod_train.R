
#======================================[ Differential expression analysis use the rankprod package ]==========================================
library(RankProd)
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  

exprs_h3n2<-read.csv('exprs_train.csv',row.names = 1)
h3n2_group <- read.csv("h3n2_group_train.csv",row.names = 1)


exprs.gnames<-rownames(exprs_h3n2)
h3n2_group<-h3n2_group[colnames(exprs_h3n2),]

exprs.cl <- as.numeric(as.factor(h3n2_group$group))
exprs.cl[exprs.cl==1]<-0
exprs.cl[exprs.cl==2]<-1

RP.adv.out <- RP.advance(exprs_h3n2,exprs.cl,h3n2_group$batch, logged=TRUE,na.rm=FALSE,  gene.names=exprs.gnames,
                         plot = FALSE,rand = 123,MinNumOfValidPairs =1)     ###class1:asym,class2:sym

degs<-topGene(RP.adv.out,cutoff = 0.05,method="pfp",logged=TRUE,logbase=2,gene.names=exprs.gnames) 
down<-degs$Table1
up<-degs$Table2
commgene <- intersect(rownames(up),rownames(down))

up<-up[!rownames(up)%in%commgene,]
down<-down[!rownames(down)%in%commgene,]
write.csv(up,'updeg_h3n2_train.csv')  
write.csv(down,'downdeg_h3n2_train.csv')  
geneFC <- RP.adv.out$AveFC

stable_gene <- data.frame(gene=setdiff(setdiff(exprs.gnames,rownames(up)),rownames(down)),type="unchanged")
upgene <- data.frame(gene=setdiff(rownames(up),rownames(down)),type="up")
downgene <- data.frame(gene=setdiff(rownames(down),rownames(up)),type="down")
#stable_gene1 <- data.frame(gene=intersect(rownames(up),rownames(down)),type="unchanged")
allgene <- rbind(upgene,downgene,stable_gene)#,stable_gene1)

colnames(geneFC)<-c("AveFC")
write.csv(geneFC,"gene_fc_train.csv",row.names = T)


#================================[ make heatmap of DEG-rankprod ]==============================

  library(ggplot2)
  library(pheatmap)
  library(ComplexHeatmap)
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  

  
  up_deg<-read.csv('updeg_h3n2_train.csv',row.names = 1)  ##477
  down_deg<-read.csv('downdeg_h3n2_train.csv',row.names = 1)  ##273
  degs<-rbind(up_deg,down_deg)
  exprs<-read.csv('exprs_train.csv',row.names = 1)  ##8286*72
  gene_pfp <- RP.adv.out$pfp
  gene_pfp<- as.data.frame(gene_pfp)
  gene_pfp$gene <- rownames(gene_pfp)
  gene_pval <- RP.adv.out$pval
  gene_pval<- as.data.frame(gene_pval)
  gene_pval$gene <- rownames(gene_pval)
  
  geneFC <- read.csv("gene_fc_train.csv",header = T)
  colnames(geneFC)<- c("gene","aveFC")
  
  alldata <- cbind(gene_pfp,gene_pval,geneFC)
  myfun <- function(x){
    y1<-x[c(3,1,4,8)]
    y2<-x[c(3,2,5,8)]
    if(as.numeric(x[8])<0){
      y=y1
    }else{
      y=y2
    }
    return(t(y))}
  tmp <- as.data.frame(t(apply(alldata,1,myfun)))
  colnames(tmp)<- c("gene","pfp","pval","aveFC")
  tmp$pfp <- as.numeric(tmp$pfp)
  tmp$pval <- as.numeric(tmp$pval)
  tmp$aveFC <- as.numeric(tmp$aveFC)
  tmp<- merge(tmp,allgene,by=c("gene"))
  
  ggplot()+geom_point(aes(x= tmp$aveFC,
                          y=-log10(tmp$pfp),color=tmp$type))+#,color=CB7vsCB13_select$State
    mytheme+geom_hline(aes(yintercept=-log10(0.05)),colour="#526373",linetype="dashed")+
    #geom_hline(aes(yintercept=-0.58),colour="#526373",linetype="dashed")+
    geom_vline(aes(xintercept=0),colour="#526373",linetype="dashed")+
    scale_color_manual(values=c("#085A9C","#526373","#EF0808"))+
    labs(x="Log2FoldChange",y="-Log10pfp",colour="",fill="") +ylim(0,20)+xlim(-1.2,1.2)
  
  