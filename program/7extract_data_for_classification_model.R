
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
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  

###all data
exprs<-read.csv('h3n2_combat_data1.csv',row.names = 1)
h3n2_group <- read.csv("group_h3n21.csv",row.names = 1)
#exprs<-exprs[!(rownames(exprs)%in%c("EIF1AY","KDM5D","DDX3Y","PRS4Y1","PRKY")),]


###############
h3n2_group_4_dataset<-subset(h3n2_group,!(batch=="GSE61754"))#32, #38
exprs_4_dataset <- exprs[,h3n2_group_4_dataset$subject]

h3n2_group_4_dataset$subject <- factor(h3n2_group_4_dataset$subject,levels = colnames(exprs_4_dataset))
group_4_dataset <-h3n2_group_4_dataset[colnames(exprs_4_dataset),]
group_4_dataset<-as.numeric(factor(group_4_dataset$group))
group_4_dataset[group_4_dataset==1]<-0
group_4_dataset[group_4_dataset==2]<-1

#1 1 0 1 0 1 1 0 1 0 0 0 1 1 1 1 0 1 0 1 1 0 1 0 0 1 1 1 0 0 0 1 0 1 0 0 0 1 0 0 1 1 1 1 1 0 1 1 0 1 1 1 0
#1 0 0 0 1 1 1 1 0 1 0 1 1 0 1 0 0
write.csv(exprs_4_dataset,"exprs_4_dataset.csv")
write.csv(group_4_dataset,"group_4_dataset.csv")
write.csv(h3n2_group_4_dataset,"h3n2_group_4_dataset.csv")

h3n2_group_1_dataset<-subset(h3n2_group,batch=="GSE61754")#11, #14
exprs_1_dataset <- exprs[,h3n2_group_1_dataset$subject]

h3n2_group_1_dataset$subject <- factor(h3n2_group_1_dataset$subject,levels = colnames(exprs_1_dataset))
group_1_dataset <-h3n2_group_1_dataset[colnames(exprs_1_dataset),]
group_1_dataset<-as.numeric(factor(group_1_dataset$group))
group_1_dataset[group_1_dataset==1]<-0
group_1_dataset[group_1_dataset==2]<-1
#1 0 0 0 1 1 1 1 0 1 0 1 1 0 0 0 0
write.csv(exprs_1_dataset,"exprs_1_dataset.csv")
write.csv(group_1_dataset,"group_1_dataset.csv")
write.csv(h3n2_group_1_dataset,"h3n2_group_1_dataset.csv")


##########all 87 data set 
h3n2_group$subject <- factor(h3n2_group$subject,levels = colnames(exprs))
group_all <-h3n2_group[colnames(exprs),]
group_all<-as.numeric(factor(group_all$group))
group_all[group_all==1]<-0
group_all[group_all==2]<-1

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
#lightgreen<-read.csv('lightgreengenes_all.csv',row.names = 1)
#lightcyan<-read.csv('lightcyangenes_all.csv',row.names = 1)
cyangenes<-read.csv('cyangenes_all.csv',row.names = 1)
#bluegenes <- read.csv("bluegenes_all.csv",row.names = 1)
#purplegenes <- read.csv("purplegenes_all.csv",row.names = 1)
#salmongenes <- read.csv("salmongenes_all.csv",row.names = 1)
#greenyellowgenes <- read.csv("greenyellowgenes_all.csv",row.names = 1)
#magentagenes <- read.csv("magentagenes_all.csv",row.names = 1)
browngenes <- read.csv("browngenes_all.csv",row.names = 1)
#yellowgenes <- read.csv("yellowgenes_all.csv",row.names = 1)

#greengenes <- read.csv("greengenes.csv",row.names = 1)
#midnightbluegenes <- read.csv("midnightbluegenes.csv",row.names = 1)

#redgenes <- read.csv("redgenes.csv",row.names = 1)
tangenes <- read.csv("tangenes_all.csv",row.names = 1)
#blackgenes <- read.csv("blackgenes.csv",row.names = 1)
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

#lightgreengenes<-expression_fun(lightgreen,"MELightgreen")
#bluegenes<-expression_fun(bluegenes,"MEBlue")

#purplegenes <- expression_fun(purplegenes,"MEPurple")
#salmongenes <- expression_fun(salmongenes,"MESalmon")
#greenyellowgenes <- expression_fun(greenyellowgenes,"MEGreenyellow")
#magentagenes <- expression_fun(magentagenes,"MEMagenta")
browngenes <- expression_fun(browngenes,"MEbrown")
#yellowgenes <- expression_fun(yellowgenes,"MEYellow")
#greengenes <- expression_fun(greengenes,"MEGreen")
#midnightbluegenes <- expression_fun(midnightbluegenes,"MEMidnightblue")
#redgenes <- expression_fun(redgenes,"MERed")
tangenes <- expression_fun(tangenes,"MEtan")
cyangenes <- expression_fun(cyangenes,"MEcyan")

#blackgenes <- expression_fun(blackgenes,"MEBlack")
#pinkgenes <- expression_fun(pinkgenes,"MEPink")
#turquoisegenes <- expression_fun(turquoisegenes,"METurquoise")
#greygenes <- expression_fun(greygenes,"MEGrey")


#trait_module <- rbind(magentagenes,salmongenes,lightgreengenes)#yellowgenes,
#trait_module <- rbind(browngenes,bluegenes,yellowgenes,tangenes,salmongenes,greenyellowgenes)
trait_module <- rbind(browngenes,tangenes,cyangenes)#yellowgenes,

write.csv(trait_module,"trait_module_all.csv",row.names = T)

################################################
module_gene <-read.csv("trait_module_all.csv",row.names = 1)
dm_intersection <- intersect(rownames(degs),rownames(module_gene))
dm_merge <- unique(union(rownames(degs),rownames(module_gene)))

dm_inter <- exprs[dm_intersection,]
dm_me <- exprs[dm_merge,]

write.csv(dm_inter,"dm_inter_all.csv",row.names = T)
write.csv(dm_me,"dm_me_all.csv",row.names = T)

dm_model<-module_gene[dm_intersection,][,c(1,2)]
write.csv(dm_model,"dm_inter_gene_module_all.csv",row.names = T)



############################################################4 dataset 
exprs<-read.csv('exprs_4_dataset.csv',row.names = 1)
h3n2_group <- read.csv("h3n2_group_4_dataset.csv",row.names = 1)

exprs1<-read.csv('exprs_1_dataset.csv',row.names = 1)
h3n2_group1 <- read.csv("h3n2_group_1_dataset.csv",row.names = 1)

###############extract deg
updegs<-read.csv('updeg_h3n2_4_dataset.csv',row.names = 1)  ##477
downdegs<-read.csv('downdeg_h3n2_4_dataset.csv',row.names = 1)  ##273
degs<-rbind(updegs,downdegs)
deg_expression <- exprs[rownames(degs),]
write.csv(deg_expression,"deg_expression_4_dataset.csv")
##################extract the expression of gene in the co-expression modules

#cyangenes <- read.csv("cyangenes_4dataset.csv",row.names = 1)
#bluegenes <- read.csv("bluegenes_4dataset.csv",row.names = 1)
#purplegenes <- read.csv("purplegenes_4dataset.csv",row.names = 1)
#salmongenes <- read.csv("salmongenes_4dataset.csv",row.names = 1)
#greenyellowgenes <- read.csv("greenyellowgenes_4dataset.csv",row.names = 1)
magentagenes <- read.csv("magentagenes_4dataset.csv",row.names = 1)
#browngenes <- read.csv("browngenes.csv",row.names = 1)
#yellowgenes <- read.csv("yellowgenes_4dataset.csv",row.names = 1)

greengenes <- read.csv("greengenes_4dataset.csv",row.names = 1)
midnightbluegenes <- read.csv("midnightbluegenes_4dataset.csv",row.names = 1)

#redgenes <- read.csv("redgenes_4dataset.csv",row.names = 1)
tangenes <- read.csv("tangenes_4dataset.csv",row.names = 1)
blackgenes <- read.csv("blackgenes_4dataset.csv",row.names = 1)
#pinkgenes <- read.csv("pinkgenes.csv",row.names = 1)
#turquoisegenes <- read.csv("turquoisegenes_4dataset.csv",row.names = 1)
#greygenes <- read.csv("greygenes.csv",row.names = 1)

expression_fun <- function(x,module_name){
  colnames(x)<- c("gene")
  x$modeuleType <-module_name
  x<- cbind(x,exprs[rownames(exprs) %in% x$gene,])
  rownames(x)<-x$gene
  write.csv(x,paste(module_name,"_expression_4dataset.csv"),row.names = T)
  return(x)}

#cyangenes<-expression_fun(cyangenes,"MECyan")
#bluegenes<-expression_fun(bluegenes,"MEBlue")

#purplegenes <- expression_fun(purplegenes,"MEPurple")
#salmongenes <- expression_fun(salmongenes,"MESalmon")
#greenyellowgenes <- expression_fun(greenyellowgenes,"MEGreenyellow")
magentagenes <- expression_fun(magentagenes,"MEMagenta")
#browngenes <- expression_fun(browngenes,"MEBrown")
#yellowgenes <- expression_fun(yellowgenes,"MEYellow")
greengenes <- expression_fun(greengenes,"MEGreen")
midnightbluegenes <- expression_fun(midnightbluegenes,"MEMidnightblue")
#redgenes <- expression_fun(redgenes,"MERed")
tangenes <- expression_fun(tangenes,"METan")

blackgenes <- expression_fun(blackgenes,"MEBlack")
#pinkgenes <- expression_fun(pinkgenes,"MEPink")
#turquoisegenes <- expression_fun(turquoisegenes,"METurquoise")
#greygenes <- expression_fun(greygenes,"MEGrey")


trait_module <- rbind(magentagenes,tangenes)

write.csv(trait_module,"trait_module_4dataset.csv",row.names = T)


################################################
module_gene <-read.csv("trait_module_4dataset.csv",row.names = 1)
dm_intersection <- intersect(rownames(degs),rownames(module_gene))
dm_merge <- unique(union(rownames(degs),rownames(module_gene)))

dm_inter <- exprs[dm_intersection,]
dm_me <- exprs[dm_merge,]

write.csv(dm_inter,"dm_inter_4dataset.csv",row.names = T)
write.csv(dm_me,"dm_me_4dataset.csv",row.names = T)

dm_inter <- exprs1[dm_intersection,]
dm_me <- exprs1[dm_merge,]

write.csv(dm_inter,"dm_inter_1dataset.csv",row.names = T)
write.csv(dm_me,"dm_me_1dataset.csv",row.names = T)

dm_model<-module_gene[dm_intersection,][,c(1,2)]
write.csv(dm_model,"dm_inter_gene_module_4dataset.csv",row.names = T)

###################module compare between all and train
library(pheatmap)

module_gene_all <- read.csv("dm_inter_gene_module_all.csv",row.names = 1)
module_gene_train <- read.csv("dm_inter_gene_module_4dataset.csv",row.names = 1)
df_all<-data.frame(number=numeric(),len_all=numeric(),len_train=numeric(),len_union=numeric()
                   ,type=character())
for (i in unique(module_gene_all$modeuleType)){
  for (j in unique(module_gene_train$modeuleType)){
    tmp1 <- subset(module_gene_all,modeuleType==i)
    tmp2<- subset(module_gene_train,modeuleType==j)
    
    df=data.frame(number=length(intersect(tmp1$gene,tmp2$gene)),len_all=nrow(tmp1),len_train=nrow(tmp2),len_union=length(union(tmp1$gene,tmp2$gene))
                 ,type=paste(i,j,sep = '_'))
    df_all<-rbind(df_all,df)
  }
}
df_all<-na.omit(df_all)
df_all$dis <- df_all$number/df_all$len_union

module_gene_train <- read.csv("dm_inter_gene_module_4dataset_most.csv",row.names = 1)
df_all1<-data.frame(number=numeric(),len_all=numeric(),len_train=numeric(),len_union=numeric()
                   ,type=character())
for (i in unique(module_gene_all$modeuleType)){
  for (j in unique(module_gene_train$modeuleType)){
    tmp1 <- subset(module_gene_all,modeuleType==i)
    tmp2<- subset(module_gene_train,modeuleType==j)
    df=data.frame(number=length(intersect(tmp1$gene,tmp2$gene)),len_all=nrow(tmp1),len_train=nrow(tmp2),len_union=length(union(tmp1$gene,tmp2$gene))
                  ,type=paste(i,j,sep = '_'))
    df_all1<-rbind(df_all1,df)
  }
}
df_all<-na.omit(df_all)
df_all1$dis <- df_all1$number/df_all1$len_union

sim_df <- data.frame(matrix(df_all1$dis,nrow = 5))
colnames(sim_df)<-c("MEbrown","MEtan","MEcyan")
rownames(sim_df)<-c("MEmagenta","MEmidnightblue","MEblack","MEgreen","MEtan")

anno_row <- data.frame("Module_train"=rownames(sim_df))
rownames(anno_row)<-rownames(sim_df)

anno_col <- data.frame("Module_all"=colnames(sim_df))
rownames(anno_col)<-colnames(sim_df)


ann_colors <- list(Module_all=c(MEbrown="#7C3F07",MEmagenta="#FF80FFFF",
                                MEblack="black",MEcyan="cyan",
                                MEmidnightblue="midnightblue",MEgreen="green" ,
                                MEtan="tan"))

pheatmap(sim_df,cluster_rows=F,cluster_cols=F,
         display_numbers=T,number_format="%.2f",
         border="white",
         #scale="row",
         fontsize_number=8,
         fontsize_col = 8,
         fontsize_row = 8,
         #angle_col = 90,
         legend_breaks=c(-0.35,-0.2,0,0.2,0.35),
         legend_labels=c("-0.35","-0.2","0","0.2","0.35"),
         #color = mycol,
         fontsize=8,
         annotation_row=anno_row,
         annotation_col=anno_col,
         annotation_colors=ann_colors
)

