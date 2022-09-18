setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  
library(ggplot2)
library(reshape)
library(reshape2)
library(VennDiagram)
library(plotROC)

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
                          strip.text=element_text(size=16,colour="#085A9C",
                          ),#family="CA"),
                          strip.background=element_blank())
mytheme2 <- theme(
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 13),
  legend.title = element_blank(),
  title = element_text(size = 20))+
  theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
panel.background = element_blank(),panel.border=element_rect(fill='transparent',
color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12),
                                         axis.text.y=element_text(size=12),
                                         legend.background = element_blank())
########
GOPlotfun <- function(goresult,picturename1,picturename2){
  #goresult <- goresult[1:10,]
  #goresult<-term_numfun(goresult)
  goresult<-as.data.frame(t(apply(goresult,1,term_lengthfun)))
  goresult$logPadjust<-as.numeric(as.character(goresult$logPadjust))
  goresult$Count<-as.numeric(as.character(goresult$Count))
  goresult<-goresult[order(goresult$logPadjust,decreasing = T),]
  goresult$Description <- factor(goresult$Description,levels = rev(goresult$Description))
  goresult$geneRatio <- apply(goresult,1,generationfun)
  goresult$geneRatio<-as.numeric(as.character(goresult$geneRatio))
  
  pp<-ggplot(goresult,aes(logPadjust,Description)) + 
    geom_point(aes(color=geneRatio,size=Count)) +
    labs(x="- Log10Padjust",y="",fill="") +
    scale_colour_gradient(low="green",high="red")+
    theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
          panel.background = element_blank(),
          panel.border=element_rect(fill='transparent', color='black'), 
          axis.ticks.x = element_blank(),
          axis.title.x  = element_text(size = 30),
          axis.text.x=element_text(size=24), axis.text.y=element_text(size=30),
          legend.background = element_blank())
  
  
  pdf(paste(picturename1,picturename2,"bubble.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
  
  pp<- ggplot(goresult,aes(x=Description,y=logPadjust,fill=ONTOLOGY)) + 
    geom_bar(position=position_dodge(0.5), stat="identity",width = 0.45) +
    labs(x = " ",y="- Log10Padjust",fill="")+
    scale_fill_manual(values=c("#085A9C","#D4C009","#CC3333"))+
    theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"), 
          panel.background = element_blank(),
          panel.border=element_rect(fill='transparent', color='black'), 
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 30),
          axis.text.x=element_text(size=24), axis.text.y=element_text(size=30),
          legend.position = "right")+
    coord_flip()
  
  pdf(paste(picturename1,picturename2,"bar.pdf",sep = "_"),height=6,width=12)
  print(pp)
  dev.off()
}
##################ROC curve
plotROC <- function(.data, predict_col, target, group, picturename,positive=1, all=TRUE){
  if(!(require(tidyverse) & require(plotROC))){
    stop("--> tidyverse and plotROC packages are required..")
  } 
  
  predict_col <- enquo(predict_col)
  target <- enquo(target)
  group  <- enquo(group)
  
  predictN <- quo_name(predict_col)
  groupN   <- quo_name(group)
  
  df <- .data %>% dplyr::select(!! predict_col, !! target, !! group) %>%
    mutate(targetN = ifelse(!! target == positive, 1, 0)) %>% as.data.frame()
  if (all){
    df2 <- df 
    df2[, groupN] <- "ALL"
    
    df <- rbind(df, df2)
  }
  p  <- df %>%  ggplot(aes_string(m = predictN, 
                                  d = "targetN",
                                  color = groupN)) + geom_roc(show.legend = TRUE, labels=FALSE)
  p <- p + mytheme2+theme(legend.position = "none")+scale_color_manual(values = c("#085A9C","#CC3333","#D4C009"))
  
  ng <- levels(factor(df[, groupN]))
  if(length(ng) == 3){
    auc <- calc_auc(p)$AUC
    names(auc) <- ng
    auc <- base::sort(auc, decreasing = TRUE)
    p <- p + annotate("text", x = .75, y = .1, 
                      label = paste(names(auc)[1], " AUC =", round(auc[1], 3), "\n"),
                      size = 6)
  }
  
  pp <- p + xlab("1 - Specificity") + ylab("Sensitivity") + 
    scale_x_continuous(expand = c(0, 0),limits = c(0, 1.1)) + scale_y_continuous(expand = c(0, 0),limits = c(0, 1.1))
  pdf(paste(picturename,"roc.pdf",sep = "_"),height=5,width=5)
  print(pp)
  dev.off()
}

####all_data
group<-c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0)
name <- ifelse(group==0,"Asymptomatic","Symptomatic")
score <-c(0.78, 0.99 ,0.88, 0.76, 0.95, 1.,   0.99, 0.8,  0.28, 0.09, 0.04, 0.25, 0.1,  0.21,
          0.19 ,0.11 ,0.6  ,0.98, 0.86, 0.91, 0.98, 0.98, 0.98, 0.92, 0.89, 0.26, 0.08, 0.07,
          0.11 ,0.24, 0.04 ,0.16, 0.16, 0.84, 0.93, 0.74, 0.94, 0.96, 0.01, 0.05, 0.17)
roc_df <- data.frame(group,score,name)


plotROC(roc_df, predict_col = score, target = group, group = name,"all_data", positive = 1)
########## data 4
#######
group<-c(1, 1, 0, 1, 0, 1, 1, 0)
name <- ifelse(group==0,"Asymptomatic","Symptomatic")
score <-c(0.35, 0.64, 0.13, 0.26, 0.09, 0.14, 0.95, 0.47)
roc_df <- data.frame(group,score,name)


plotROC(roc_df, predict_col = score, target = group, group = name,"data4", positive = 1)

########## data 4 _GSE30550
#######
group<-c(1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0)
name <- ifelse(group==0,"Asymptomatic","Symptomatic")
score <-c(0.78723404, 0.44680851, 0.04255319, 0.29787234, 0.76595745, 0.89361702,
          1.         ,0.5106383,  0.4893617,  0.93617021, 0.04255319, 0.95744681,
          0.85106383 ,0.10638298, 0.76595745, 0.12765957, 0.25531915)
roc_df <- data.frame(group,score,name)


plotROC(roc_df, predict_col = score, target = group, group = name,"data4_GSE30550", positive = 1)

######################## importence
###all_data
hub_importance <- read.csv("/home/dulab/Documents/wrok/C/deg_module__inter_importance_all.csv")
hub_importance<-hub_importance[,-1]
hub_importance$type <- "importance of module and deg genes"
hub_importance<-hub_importance[order(hub_importance$importance,decreasing = F),]
hub_importance$gene <- factor(hub_importance$gene,levels = hub_importance$gene)
rownames(hub_importance)<-hub_importance$gene
pp<-ggplot(hub_importance,aes(x=gene,y=importance,fill=-importance))+geom_bar(stat="identity")+ labs(x = " ",y="Importance",fill="")+
  mytheme+coord_flip()
pdf(paste("importance","all.pdf",sep = "_"),height=5,width=5)
print(pp)
dev.off()
backgene<-read.csv('h3n2_combat_data1.csv',row.names = 1)   #72*8286
allgofun(hub_importance,backgene,"deg_module_inter_importance_all")
goresult<-read.csv("deg_module_inter_importance_all_go.csv",row.names = 1)   ##24
goresult$type <- "all"
#goresult<-goresult[goresult$ONTOLOGY=='BP',]   #19
keggresult<-read.csv("deg_module_inter_importance_all_kegg.csv",row.names = 1)  ##0
GOPlotfun(goresult,'dm_me_inter','all')
keggPlotfun(keggresult,'KEGG_dm_me_inter','all')

#########data4
hub_importance <- read.csv("/home/dulab/Documents/wrok/C/deg_module__inter_importance_dataset4.csv")
hub_importance<-hub_importance[,-1]
hub_importance$type <- "importance of module and deg genes"
hub_importance<-hub_importance[order(hub_importance$importance,decreasing = F),]
hub_importance$gene <- factor(hub_importance$gene,levels = hub_importance$gene)
rownames(hub_importance)<-hub_importance$gene
pp<-ggplot(hub_importance,aes(x=gene,y=importance,fill=-importance))+geom_bar(stat="identity")+ labs(x = " ",y="Importance",fill="")+
  mytheme+coord_flip()
pdf(paste("importance","data4.pdf",sep = "_"),height=5,width=5)
print(pp)
dev.off()



