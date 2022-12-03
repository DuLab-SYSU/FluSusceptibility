
library(stringr)
library(ggplot2)
library(grid)
library(futile.logger)
library(reshape)
library(reshape2)
library(VennDiagram)
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
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  

###all data
exprs<-read.csv('h3n2_combat_data.csv',row.names = 1)
h3n2_group <- read.csv("group_h3n2.csv",row.names = 1)
#exprs<-exprs[!(rownames(exprs)%in%c("EIF1AY","KDM5D","DDX3Y","PRS4Y1","PRKY")),]

#exprs_train<-read.csv('exprs00_h3n2_train.csv',row.names = 1)
#h3n2_group_train <- read.csv("group_h3n2_train.csv",row.names = 1)

###############
h3n2_group_train<-subset(h3n2_group,!(batch=="GSE61754"))#32, #38
exprs_train <- exprs[,h3n2_group_train$gseID]

h3n2_group_train$gseID <- factor(h3n2_group_train$gseID,levels = colnames(exprs_train))
group_train <-h3n2_group_train[colnames(exprs_train),]
group_train<-as.numeric(factor(group_train$group))
group_train[group_train==1]<-0
group_train[group_train==2]<-1

#1 1 0 1 0 1 1 0 1 0 0 0 1 1 1 1 0 1 0 1 1 0 1 0 0 1 1 1 0 0 0 1 0 1 0 0 0 1 0 0 1 1 1 1 1 0 1 1 0 1 1 1 0
#1 0 0 0 1 1 1 1 0 1 0 1 1 0 1 0 0
write.csv(exprs_train,"exprs_train.csv")
write.csv(group_train,"group_train.csv")
write.csv(h3n2_group_train,"h3n2_group_train.csv")

h3n2_group_test<-subset(h3n2_group,batch=="GSE61754")#11, #14
exprs_test <- exprs[,h3n2_group_test$gseID]

#exprs_test<-read.csv('gse61754_sym_exprs',row.names = 1)#11, #14
#h3n2_group_test <- subset(h3n2_group,batch=="GSE61754")

h3n2_group_test$gseID <- factor(h3n2_group_test$gseID,levels = colnames(exprs_test)[-12])
group_test <-h3n2_group_test[colnames(exprs_test)[-12],]
group_test<-as.numeric(factor(group_test$group))
group_test[group_test==1]<-0
group_test[group_test==2]<-1
#1 0 0 0 1 1 1 1 0 1 0 1 1 0 0 0 0
write.csv(exprs_test[,-12],"exprs_test.csv")
write.csv(group_test,"group_test.csv")
write.csv(h3n2_group_test,"h3n2_group_test.csv")


##########all 87 data set 
h3n2_group$gseID <- factor(h3n2_group$gseID,levels = colnames(exprs))
group_all <-h3n2_group[colnames(exprs),]
group_all<-as.numeric(factor(group_all$group))
group_all[group_all==1]<-0
group_all[group_all==2]<-1


#sym_group <- h3n2_group[h3n2_group$group=="Symptomatic",]
#asym_group <- h3n2_group[h3n2_group$group=="Asymptomatic",]

###############extract deg
updegs<-read.csv('updeg_h3n2_all.csv',row.names = 1)  ##477
downdegs<-read.csv('downdeg_h3n2_all.csv',row.names = 1)  ##273

#updegs<-read.csv('updeg_h3n2_all_nosex.csv',row.names = 1)  ##477
#downdegs<-read.csv('downdeg_h3n2_all_nosex.csv',row.names = 1)  ##273
degs<-rbind(updegs,downdegs)
deg_expression <- exprs[rownames(degs),]
#write.csv(deg_expression,"deg_expression_all_nosex.csv")
write.csv(deg_expression,"deg_expression_all.csv")

##################extract the expression of gene in the co-expression modules
lightgreen<-read.csv('lightgreengenes_all.csv',row.names = 1)
lightcyan<-read.csv('lightcyangenes_all.csv',row.names = 1)
#cyangenes<-read.csv('cyangenes_all.csv',row.names = 1)
#bluegenes <- read.csv("bluegenes_all.csv",row.names = 1)
purplegenes <- read.csv("purplegenes_all.csv",row.names = 1)
#salmongenes <- read.csv("salmongenes_all.csv",row.names = 1)
#greenyellowgenes <- read.csv("greenyellowgenes_all.csv",row.names = 1)
#magentagenes <- read.csv("magentagenes_all.csv",row.names = 1)
#browngenes <- read.csv("browngenes_all.csv",row.names = 1)
#yellowgenes <- read.csv("yellowgenes_all.csv",row.names = 1)

#greengenes <- read.csv("greengenes.csv",row.names = 1)
#midnightbluegenes <- read.csv("midnightbluegenes.csv",row.names = 1)

#redgenes <- read.csv("redgenes.csv",row.names = 1)
#tangenes <- read.csv("tangenes_all.csv",row.names = 1)
blackgenes <- read.csv("blackgenes_all.csv",row.names = 1)
#pinkgenes <- read.csv("pinkgenes.csv",row.names = 1)
#turquoisegenes <- read.csv("turquoisegenes.csv",row.names = 1)
#greygenes <- read.csv("greygenes.csv",row.names = 1)

expression_fun <- function(x,module_name){
  colnames(x)<- c("gene")
  x$modeuleType <-module_name
  x<- cbind(x,exprs[rownames(exprs) %in% x$gene,])
  rownames(x)<-x$gene
  write.csv(x,paste(module_name,"_expression_all.csv"),row.names = T)
  return(x)}

lightgreengenes<-expression_fun(lightgreen,"MELightgreen")
lightcyangenes<-expression_fun(lightcyan,"MELightcyan")
#bluegenes<-expression_fun(bluegenes,"MEBlue")

purplegenes <- expression_fun(purplegenes,"MEPurple")
#salmongenes <- expression_fun(salmongenes,"MESalmon")
#greenyellowgenes <- expression_fun(greenyellowgenes,"MEGreenyellow")
#magentagenes <- expression_fun(magentagenes,"MEMagenta")
#browngenes <- expression_fun(browngenes,"MEbrown")
#yellowgenes <- expression_fun(yellowgenes,"MEYellow")
#greengenes <- expression_fun(greengenes,"MEGreen")
#midnightbluegenes <- expression_fun(midnightbluegenes,"MEMidnightblue")
#redgenes <- expression_fun(redgenes,"MERed")
#tangenes <- expression_fun(tangenes,"MEtan")
#cyangenes <- expression_fun(cyangenes,"MEcyan")

#blackgenes <- expression_fun(blackgenes,"MEBlack")
#pinkgenes <- expression_fun(pinkgenes,"MEPink")
#turquoisegenes <- expression_fun(turquoisegenes,"METurquoise")
#greygenes <- expression_fun(greygenes,"MEGrey")


#trait_module <- rbind(magentagenes,salmongenes,lightgreengenes)#yellowgenes,
#trait_module <- rbind(browngenes,bluegenes,yellowgenes,tangenes,salmongenes,greenyellowgenes)
#trait_module <- rbind(browngenes,tangenes,cyangenes)#yellowgenes,
#trait_module <- rbind(lightgreengenes,lightcyangenes,purplegenes,blackgenes)#yellowgenes,

#write.csv(trait_module,"trait_module_all.csv",row.names = T)

trait_module <- rbind(lightgreengenes,lightcyangenes,purplegenes)#yellowgenes,

write.csv(trait_module,"trait_module_all_noblack.csv",row.names = T)


#dm_inter_sym <- dm_inter[,sym_group$gseID]
#dm_inter_asym <- dm_inter[,asym_group$gseID]
#write.csv(dm_inter_sym,"dm_inter_sym_all.csv",row.names = T)
#write.csv(dm_inter_asym,"dm_inter_asym_all.csv",row.names = T)

module_gene <-read.csv("trait_module_all_noblack.csv",row.names = 1)
dm_intersection <- intersect(rownames(degs),rownames(module_gene))
dm_merge <- unique(union(rownames(degs),rownames(module_gene)))

dm_inter <- exprs[dm_intersection,]
dm_me <- exprs[dm_merge,]

write.csv(dm_inter,"dm_inter_all_noblack.csv",row.names = T)
write.csv(dm_me,"dm_me_all_noblack.csv",row.names = T)

dm_model<-module_gene[dm_intersection,][,c(1,2)]
write.csv(dm_model,"dm_inter_gene_module_all_noblack.csv",row.names = T)


############################################################4 dataset 
exprs<-read.csv('exprs_train.csv',row.names = 1)
h3n2_group <- read.csv("h3n2_group_train.csv",row.names = 1)

exprs1<-read.csv('exprs_test.csv',row.names = 1)
h3n2_group1 <- read.csv("h3n2_group_test.csv",row.names = 1)

###############extract deg
updegs<-read.csv('updeg_h3n2_train.csv',row.names = 1)  ##477
downdegs<-read.csv('downdeg_h3n2_train.csv',row.names = 1)  ##273
degs<-rbind(updegs,downdegs)
deg_expression <- exprs[rownames(degs),]
write.csv(deg_expression,"deg_expression_train.csv")
##################extract the expression of gene in the co-expression modules

cyangenes <- read.csv("cyangenes_train.csv",row.names = 1)
bluegenes <- read.csv("bluegenes_train.csv",row.names = 1)
#purplegenes <- read.csv("purplegenes_train.csv",row.names = 1)
#salmongenes <- read.csv("salmongenes_train.csv",row.names = 1)
greenyellowgenes <- read.csv("greenyellowgenes_train.csv",row.names = 1)
#magentagenes <- read.csv("magentagenes_train.csv",row.names = 1)
#browngenes <- read.csv("browngenes.csv",row.names = 1)
#yellowgenes <- read.csv("yellowgenes_train.csv",row.names = 1)

#greengenes <- read.csv("greengenes_train.csv",row.names = 1)
midnightbluegenes <- read.csv("midnightbluegenes_train.csv",row.names = 1)

#redgenes <- read.csv("redgenes_train.csv",row.names = 1)
#tangenes <- read.csv("tangenes_train.csv",row.names = 1)
#blackgenes <- read.csv("blackgenes_train.csv",row.names = 1)
#pinkgenes <- read.csv("pinkgenes.csv",row.names = 1)
#turquoisegenes <- read.csv("turquoisegenes_train.csv",row.names = 1)
#greygenes <- read.csv("greygenes.csv",row.names = 1)
grey60genes <- read.csv("grey60genes_train.csv",row.names = 1)

expression_fun <- function(x,module_name){
  colnames(x)<- c("gene")
  rownames(x) <-x$gene
  x$modeuleType <-module_name
  x<- x[intersect(x$gene,rownames(exprs)),]
  x<- cbind(x,exprs[rownames(exprs) %in% x$gene,])
  rownames(x)<-x$gene
  write.csv(x,paste(module_name,"_expression_train.csv"),row.names = T)
  return(x)}

cyangenes<-expression_fun(cyangenes,"MECyan")
bluegenes<-expression_fun(bluegenes,"MEBlue")

#purplegenes <- expression_fun(purplegenes,"MEPurple")
#salmongenes <- expression_fun(salmongenes,"MESalmon")
greenyellowgenes <- expression_fun(greenyellowgenes,"MEGreenyellow")
#magentagenes <- expression_fun(magentagenes,"MEMagenta")
#browngenes <- expression_fun(browngenes,"MEBrown")
#yellowgenes <- expression_fun(yellowgenes,"MEYellow")
#greengenes <- expression_fun(greengenes,"MEGreen")
midnightbluegenes <- expression_fun(midnightbluegenes,"MEMidnightblue")
#redgenes <- expression_fun(redgenes,"MERed")
#tangenes <- expression_fun(tangenes,"METan")

#blackgenes <- expression_fun(blackgenes,"MEBlack")
#pinkgenes <- expression_fun(pinkgenes,"MEPink")
#turquoisegenes <- expression_fun(turquoisegenes,"METurquoise")
#greygenes <- expression_fun(greygenes,"MEGrey")
grey60genes <- expression_fun(grey60genes,"MEGrey60")


trait_module <- rbind(cyangenes,bluegenes,greenyellowgenes,midnightbluegenes,grey60genes)

write.csv(trait_module,"trait_module_train.csv",row.names = T)

#trait_module1 <-rbind(magentagenes,greengenes,blackgenes,midnightbluegenes,tangenes)
#write.csv(trait_module1,"most_revelent_train.csv",row.names = T)

################################################
module_gene <-read.csv("trait_module_train.csv",row.names = 1)
dm_intersection <- intersect(rownames(degs),rownames(module_gene))
dm_merge <- unique(union(rownames(degs),rownames(module_gene)))

#dm_inter <- exprs[dm_intersection,]
#dm_me <- exprs[dm_merge,]

#write.csv(na.omit(dm_inter),"dm_inter_train.csv",row.names = T)
#write.csv(dm_me,"dm_me_train.csv",row.names = T)
##################train1111

###################################################
#dm_intersection1 <- intersect(rownames(exprs1),dm_intersection)#TXNDC5 TXNDC5    MEGrey60
#dm_intersection2<- setdiff(dm_intersection1,"TXNDC5")
module_gene <-read.csv("trait_module_train.csv",row.names = 1)
dm_intersection <- intersect(rownames(degs),rownames(module_gene))
dm_merge <- unique(union(rownames(degs),rownames(module_gene)))
dm_inter <- exprs1[dm_intersection,]
dm_me <- exprs1[dm_merge,]

write.csv(dm_inter,"dm_inter_test1113.csv",row.names = T)
write.csv(dm_me,"dm_me_test1113.csv",row.names = T)


dm_inter <- exprs[dm_intersection,]
dm_me <- exprs[dm_merge,]


write.csv(na.omit(dm_inter),"dm_inter_train1113.csv",row.names = T)
write.csv(dm_me,"dm_me_train1113.csv",row.names = T)


dm_model<-module_gene[dm_intersection,][,c(1,2)]
write.csv(dm_model,"dm_inter_gene_module_train1113.csv",row.names = T)


