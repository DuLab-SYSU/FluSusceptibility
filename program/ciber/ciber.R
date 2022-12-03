library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(ggsci)

mytheme<-theme_bw()+theme(legend.position="right",
                          #panel.border=element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          plot.title=element_text(size=16,
                                                  colour="black",
                          ), #family="CA"),
                          axis.title =element_text(size=16,
                                                   colour="black",
                          ),#family="CA"),
                          legend.text=element_text(size=16,colour="black",
                          ),#family="CA"),
                          legend.key=element_blank(),
                          axis.text=element_text(size=16,colour="black",
                          ),#family="CA"),
                          strip.text=element_text(size=16,colour="#085A9C"),
                          legend.title = element_text(face = "bold",size = 16),#family="CA"),
                          strip.background=element_blank()
                         )
#https://www.jianshu.com/p/03a7440c0960
#https://cloud.tencent.com/developer/article/1622907
setwd("/home/dulab/Documents/wrok/flu_paper/data/CODE_FI/ciber")  

source("CIBERSORT.R")
LM22.file <- "LM22.txt"

ciber_fun <- function(expression,group,filename){
  source("CIBERSORT.R")
  LM22.file <- "LM22.txt"
  TME.results = CIBERSORT(LM22.file,expression,perm = 1000, QN = TRUE)
  TME.results <-TME.results[,-c(8,12,17,21)]
  write.table(TME.results, paste("TME.results.output",paste(filename,".txt",sep=""),sep = "_"), 
              sep = "\t", row.names = T, col.names = T, quote = F)
  dd1 <- TME.results %>% 
    as.data.frame() %>% 
    rownames_to_column("gseID") %>% 
    pivot_longer(cols=2:19,
                 names_to= "celltype",
                 values_to = "Proportion")
  
  dd1<- as.data.frame(dd1)
  dd2 <- merge(dd1,h3n2_group,by="gseID")
  write.csv(dd2,paste(filename,"_ciber_result.csv",sep=""),sep = "\t", row.names = T, col.names = T, quote = F)
  cols =c("#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#990000","#9900cc","#66FF66","#663300","#0000FF","#CC0033","#FF0000","#000099","#660066","#333333") 
  #"#999999",
  pp<-ggplot(dd2,aes(gseID,Proportion,fill = celltype)) +  scale_fill_manual(values = cols)+
    geom_bar(stat = "identity",alpha=0.8) +
    labs(fill = "",x = "",y = "Estiamted Proportion") + 
    scale_y_continuous(expand = c(0.01,0)) +mytheme+theme(axis.text.x=element_text(size=12,angle = 90))+facet_wrap(~group,nrow=2,scales = "free_x")
  pdf(paste(filename,"proportion.pdf",sep = "_"),height=9,width=16)
  print(pp)
  dev.off()
  
 # ggplot(dd1,aes(celltype,Proportion,fill = celltype))+geom_boxplot()+mytheme+theme(axis.text.x=element_text(size=12,angle = 90))+
   # guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")
  
  pp<-ggplot(dd2,aes(celltype,Proportion,fill=group))+geom_boxplot()+mytheme+
    theme(axis.text.x=element_text(size=12,angle = 90))+
    scale_fill_manual(values=c("#085A9C","#EF0808"))+
    guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")+
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test")#kruskal.test
  
  pdf(paste(filename,"sig.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
  
  pp<-ggplot(dd2,aes(celltype,Proportion,fill=group))+geom_boxplot()+mytheme+
    theme(axis.text.x=element_text(size=12,angle = 90))+
    scale_fill_manual(values=c("#085A9C","#EF0808"))+
    guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")+
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test",method.args = list(alternative = "greater"))#kruskal.test
  pdf(paste(filename,"greater.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
  pp<-ggplot(dd2,aes(celltype,Proportion,fill=group))+geom_boxplot()+mytheme+
    theme(axis.text.x=element_text(size=12,angle = 90))+
    scale_fill_manual(values=c("#085A9C","#EF0808"))+
    guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")+
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test",method.args = list(alternative = "less"))#kruskal.test
  pdf(paste(filename,"less.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
  sink(paste(filename,"log.txt"), split =T)
  naive_b <- subset(dd2,celltype=="B cells naive")
  tt<-kruskal.test(Proportion~group,naive_b)#Kruskal-Wallis chi-squared = 12.227, df = 1, p-value = 0.0004711
  #Kruskal-Wallis chi-squared = 13.356, df = 1, p-value = 0.0002576
  print(filename)
  print("B cells naive")
  print(tt)
  tt<-wilcox.test(Proportion~group,naive_b)#p-value = 0.0002619
  print(tt)
  
  
  pp<- ggplot(naive_b,aes(group,Proportion,fill=group))+geom_violin(width = 0.8)+mytheme+geom_boxplot(width=0.1,fill="white")+
    geom_jitter(width = 0.1)+
    theme(axis.text.x=element_text(size=12))+
    scale_fill_manual(values=c("#CC3333","#526373","#FF9418","#219431","#9C52AD"))+
    guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")+
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test",method.args = list(alternative = "greater"))
  pdf(paste(filename,"naive_b.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
  
  mon <- subset(dd2,celltype=="Monocytes")#Kruskal-Wallis chi-squared = 4.0445, df = 1, p-value = 0.04432
  tt<-kruskal.test(Proportion~group,mon)#Kruskal-Wallis chi-squared = 4.8281, df = 1, p-value = 0.028
  
  print(filename)
  print("Monocytes")
  print(tt)
  tt<-wilcox.test(Proportion~group,mon) #p-value = 0.02779
  print(tt)
  pp<-ggplot(mon,aes(group,Proportion,fill=group))+geom_violin(width = 0.8)+mytheme+geom_boxplot(width=0.1,fill="white")+
    geom_jitter(width = 0.1)+
    theme(axis.text.x=element_text(size=12))+
    scale_fill_manual(values=c("#CC3333","#526373","#FF9418","#219431","#9C52AD"))+
    guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")+
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test",method.args = list(alternative = "less"))
  pdf(paste(filename,"monocytes.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
  
  tcell <- subset(dd2,celltype=="T cells CD4 memory activated")
  tt<-kruskal.test(Proportion~group,tcell)#Kruskal-Wallis chi-squared = 8.0555, df = 1, p-value = 0.004536
  print(filename)
  print("T cells CD4 memory activated")
   print(tt)
  #Kruskal-Wallis chi-squared = 7.4728, df = 1, p-value = 0.006264
  tt<- wilcox.test(Proportion~group,tcell)#p-value = 0.006451
  print(tt)
  
  pp<-ggplot(tcell,aes(group,Proportion,fill=group))+geom_violin(width = 0.8)+mytheme+geom_boxplot(width=0.1,fill="white")+
    geom_jitter(width = 0.1)+
    theme(axis.text.x=element_text(size=12))+
    scale_fill_manual(values=c("#CC3333","#526373","#FF9418","#219431","#9C52AD"))+
    guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")+
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test",method.args = list(alternative = "greater"))
  pdf(paste(filename,"tcell.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
  
  neutrophils <- subset(dd2,celltype=="Neutrophils")
  tt<-kruskal.test(Proportion~group,neutrophils)#Kruskal-Wallis chi-squared = 3.9041, df = 1, p-value = 0.04817
  print(filename)
  print("Neutrophils")
  print(tt)
  tt<-wilcox.test(Proportion~group,neutrophils)#p-value = 0.0483
  print(tt)
  sink()
  
  pp<-ggplot(neutrophils,aes(group,Proportion,fill=group))+geom_violin(width = 0.8)+mytheme+geom_boxplot(width=0.1,fill="white")+
    geom_jitter(width = 0.1)+
    theme(axis.text.x=element_text(size=12))+
    scale_fill_manual(values=c("#CC3333","#526373","#FF9418","#219431","#9C52AD"))+
    guides(fill=guide_legend(ncol=1))+ labs(x="",y="Proportion",fill="")+
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test",method.args = list(alternative = "less"))
  pdf(paste(filename,"neutrophils.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
}
##############
#############################
exprs00_h3n2<-read.csv('/home/dulab/Documents/wrok/flu_paper/data/result_fi/h3n2_combat_data.csv',row.names = 1)  #87*8286
exprs_h3n2 <- exprs00_h3n2
write.table(exprs_h3n2,'h3n2_combat_data.txt',row.names = T,quote = F,sep = "\t")

exp.file <-"h3n2_combat_data.txt"
h3n2_group <- read.csv("/home/dulab/Documents/wrok/flu_paper/data/result_fi/group_h3n2.csv",row.names = 1)
ciber_fun(exp.file,h3n2_group,"all_combat1")
setwd("/home/dulab/Documents/wrok/flu_paper/data/CODE_FI/ciber")  

#################################################train

exprs00_h3n2<-read.csv('/home/dulab/Documents/wrok/flu_paper/data/result_fi/h3n2_combat_data_train.csv',row.names = 1)  #87*8286

exprs_h3n2 <- exprs00_h3n2
write.table(exprs_h3n2,'exprs00_h3n2_train_combat.txt',row.names = T,quote = F,sep = "\t")

exp.file <-"exprs00_h3n2_train_combat.txt"#h3n2_combat_data.txt
h3n2_group <- read.csv("/home/dulab/Documents/wrok/flu_paper/data/result_fi/group_h3n2_train.csv",row.names = 1)
ciber_fun(exp.file,h3n2_group,"train_combat1")

############################BH
