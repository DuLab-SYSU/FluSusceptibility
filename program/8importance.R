setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  
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
group<-c(1, 1 ,1 ,1 ,1 ,1 ,1 ,1 ,1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
name <- ifelse(group==0,"Asymptomatic","Symptomatic")
score <-c(0.91489362, 0.93617021, 0.91489362, 0.80851064, 0.91489362,
          0.93617021, 0.70212766, 0.95744681, 0.95744681, 0.14893617,
          0.04255319, 0.29787234, 0.19148936, 0.31914894, 0.19148936,
          0.14893617, 0.19148936, 0.87234043, 0.89361702, 0.87234043,
          0.85106383, 0.93617021, 0.89361702, 0.85106383, 0.9787234 ,
          0.82978723, 0.82978723, 0.95744681, 0.82978723, 0.91489362,
          0.27659574, 0.08510638, 0.27659574, 0.21276596, 0.08510638,
          0.29787234, 0.21276596, 0.12765957, 0.17021277, 0.14893617,
          0.42553191, 0.74468085, 0.93617021, 0.89361702, 0.80851064,
          0.70212766, 0.91489362, 0.95744681, 0.91489362)
roc_df <- data.frame(group,score,name)


plotROC(roc_df, predict_col = score, target = group, group = name,"all_data", positive = 1)
########## 
#######
group<-c(1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1)
name <- ifelse(group==0,"Asymptomatic","Symptomatic")
score <-c(0.79032258, 0.82258065 ,0.61290323, 0.62903226 ,0.59677419, 0.32258065,
          0.53225806, 0.66129032, 0.75806452 ,0.48387097, 0.74193548)
roc_df <- data.frame(group,score,name)


plotROC(roc_df, predict_col = score, target = group, group = name,"test", positive = 1)

########## data 4 _GSE30550


######################## importence
###all_data
hub_importance <- read.csv("deg_module__inter_importance_all.csv")
hub_importance<-hub_importance[,-1]
hub_importance$type <- "importance of module and deg genes"
hub_importance<-hub_importance[order(hub_importance$importance,decreasing = F),]
hub_importance$gene <- factor(hub_importance$gene,levels = hub_importance$gene)
rownames(hub_importance)<-hub_importance$gene
pp<-ggplot(hub_importance,aes(x=gene,y=importance,fill=-importance))+geom_bar(stat="identity")+ labs(x = " ",y="Importance",fill="")+
  mytheme+coord_flip()
pdf(paste("importance","_all.pdf",sep = ""),height=5,width=5)
print(pp)
dev.off()



#########data4
hub_importance <- read.csv("deg_module__inter_importance_train.csv")
hub_importance<-hub_importance[,-1]
hub_importance$type <- "importance of module and deg genes"
hub_importance<-hub_importance[order(hub_importance$importance,decreasing = F),]
hub_importance$gene <- factor(hub_importance$gene,levels = hub_importance$gene)
rownames(hub_importance)<-hub_importance$gene
pp<-ggplot(hub_importance,aes(x=gene,y=importance,fill=-importance))+geom_bar(stat="identity")+ labs(x = " ",y="Importance",fill="")+
  mytheme+coord_flip()
pdf(paste("importance","_train.pdf",sep = ""),height=5,width=5)
print(pp)
dev.off()


