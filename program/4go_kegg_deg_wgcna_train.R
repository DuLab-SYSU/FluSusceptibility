#============================================[ GO and KEGG for degs of rankprod and key module genes of WGCNA ]=========================================
rm(list = ls())


#===========================[ Load packges ]==============================
{
  library(clusterProfiler)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(topGO)
}
#============================[ build function to import GO and KEGG function ]======================
{
  ##
  allgofun <- function(intregene,backgene,picturename){
    dir.create(picturename)
    path=paste(path1,picturename,sep = "/")
    setwd(path)
    intregenes <- bitr(rownames(intregene),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
    backgenes <- bitr(rownames(backgene),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
    ##GO
    goann <- GOfun(intregenes,backgenes)
    write.csv(goann,paste(picturename,"go.csv",sep = "_"))
    ##KEGG
    keggann <- keggfun(intregenes,backgenes)
    write.csv(keggann,paste(picturename,"kegg.csv",sep = "_"))
    #keggpp <- keggPlotfun(keggann,picturename,"Kegg")
  }
  ##
  allgofun1 <- function(intregene,backgene,picturename){
    dir.create(picturename)
    path=paste(path1,picturename,sep = "/")
    setwd(path)
    intregenes <- bitr(intregene$x,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
    backgenes <- bitr(rownames(backgene),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
    ##GO
    goann <- GOfun(intregenes,backgenes)
    write.csv(goann,paste(picturename,"go.csv",sep = "_"))
    ##KEGG
    keggann <- keggfun(intregenes,backgenes)
    write.csv(keggann,paste(picturename,"kegg.csv",sep = "_"))
  }
}





#====================================[ build function to do GO and KEGG ]====================================
{
  term_lengthfun <- function(term,k=50){
    x = as.character(term[3])
    if(nchar(x)>k){x=substr(x,start=1,stop=k)}
    term[3]<-x
    return(term)}
  ##GO
  GOfun <- function(intregenes,backgenes){
    resultgo <- enrichGO(intregenes$ENTREZID, 'org.Hs.eg.db', keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, 
                         pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 5, 
                         maxGSSize = 500, readable = T,universe=as.character(backgenes$ENTREZID))##universe=as.character(backgenes$ENTREZID)
    resultgo <- as.data.frame(resultgo)
    resultgo$logPadjust <- -log10(as.numeric(resultgo$p.adjust))
    
    # resultgo$Description <- factor(resultgo$Description,levels = rev(resultgo$Description))  
    return(resultgo)
  }
  ##KEGG
  keggfun <- function(intregenes,backgenes){
    resultgo <- enrichKEGG(intregenes$ENTREZID,  organism = 'hsa', keyType = 'kegg', universe=as.character(backgenes$ENTREZID),
                           pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                           minGSSize = 5,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
    resultgo <- as.data.frame(resultgo)
    resultgo$logPadjust <- -log10(as.numeric(resultgo$p.adjust))
    
    #resultgo$Description <- factor(resultgo$Description,levels = rev(resultgo$Description))
    
    
    #resultgo$geneRatio <- apply(resultgo,1,generationfun)
    return(resultgo)
  }
}

#==================================[ draw bar picture and dot picture ]==================================
###GO PLOT
##adjust geneRatio value to decimal
generationfun <- function(x){
  a<- strsplit(x[4],split = "/")
  a<- round(as.numeric(a[[1]][1])/as.numeric(a[[1]][2]),digits =2)#paste(,"%",sep = "")
  return(a)
}

##GO PLOT
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
  
  
  pdf(paste(picturename1,picturename2,"bubble.pdf",sep = "_"),height=16,width=20)
  print(pp)
  dev.off()
  png(paste(picturename1,picturename2,"bubble.png",sep = "_"),width = 1300, height = 600)
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
  
  pdf(paste(picturename1,picturename2,"bar.pdf",sep = "_"),height=20,width=16)
  print(pp)
  dev.off()
  png(paste(picturename1,picturename2,"bar.png",sep = "_"),width = 1500, height = 500)
  print(pp)
  dev.off()
}



##KEGG PLOT
generationfun1 <- function(x){
  a<- strsplit(x[3],split = "/")
  a<- round(as.numeric(a[[1]][1])/as.numeric(a[[1]][2]),digits =2)#paste(,"%",sep = "")
  return(a)
}
keggPlotfun <- function(keggresult,picturename1,picturename2){
  #goresult <- ifelse(nrow(goresult)<=20,goresult,goresult[1:20,])
  #goresult<-term_numfun(goresult)
  keggresult<-as.data.frame(t(apply(keggresult,1,term_lengthfun)))
  keggresult$logPadjust<-as.numeric(keggresult$logPadjust)
  keggresult$Count<-as.numeric(keggresult$Count)
  keggresult<-keggresult[order(keggresult$logPadjust,decreasing = T),]
  keggresult$Description <- factor(keggresult$Description,levels = rev(keggresult$Description))
  keggresult$geneRatio <- apply(keggresult,1,generationfun1)
  keggresult$geneRatio<-as.numeric(as.character(keggresult$geneRatio))
  
  pp<-ggplot(keggresult,aes(logPadjust,Description)) + geom_point(aes(color=geneRatio,size=Count)) +
    labs(x="- Log10Padjust",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                                         panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                                      color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=18),
                                         axis.text.y=element_text(size=18),
                                         legend.background = element_blank())+scale_colour_gradient(low="green",high="red")+labs(x = "- Log10Padjust",y="",fill="")
  
  pdf(paste(picturename1,picturename2,"bubble.pdf",sep = "_"),height=6,width=8)
  print(pp)
  dev.off()
  png(paste(picturename1,picturename2,"bubble.png",sep = "_"),width = 800, height = 400)
  print(pp)
  dev.off()
  
  
  pp<- ggplot(keggresult,aes(x=Description,y=logPadjust,fill="class")) + geom_bar(position=position_dodge(), stat="identity",width = 0.45) +
    scale_fill_manual(values=c("#085A9C"))+theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                                                 panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                                              color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=18),
                                                 axis.text.y=element_text(size=18),
                                                 legend.position = "none")+labs(x = " ",y="- Log10Padjust",fill="")+coord_flip()
  pdf(paste(picturename1,picturename2,"bar.pdf",sep = "_"),height=6,width=8)
  print(pp)
  dev.off()
  png(paste(picturename1,picturename2,"bar.png",sep = "_"),width = 800, height = 400)
  print(pp)
  dev.off()
}


#############################################
term_lengthfun <- function(term,k=50){
  x = as.character(term[3])
  if(nchar(x)>k){x=substr(x,start=1,stop=k)}
  term[3]<-x
  return(term)}






#=============================[ setting default and load data ]==================================

##degs
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  
#setwd("~/result1/1deg/")
path1<-getwd()
gene_up<-read.csv('updeg_h3n2_train.csv',row.names = 1)  ##478
gene_down<-read.csv('downdeg_h3n2_train.csv',row.names = 1)  ##267
backgene<-read.csv('exprs_train.csv',row.names = 1)   #72*8286
pic1<-'up_dataset4'
pic2<-'down_dataset4'



#============================[ import function to do GO and KEGG and draw ]=============================
###deg
{
  ##rankprod up
  allgofun(gene_up,backgene,pic1)
  goresult<-read.csv("up_dataset4_go.csv",row.names = 1)   ##51
  goresult<-goresult[goresult$ONTOLOGY=='BP',]   #35
  keggresult<-read.csv("up_dataset4_kegg.csv",row.names = 1)  ##0
  
  GOPlotfun(goresult,'up_dataset4','BP')
  keggPlotfun(keggresult,'KEGG','up_dataset4')
  
  ##rankprod down
  setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  
  path1<-getwd()
  allgofun(gene_down,backgene,pic2)
  goresultdown<-read.csv("down_dataset4_go.csv",row.names = 1)   ###32
  goresultdown<-goresultdown[goresultdown$ONTOLOGY=='BP',]
  keggresult1<-read.csv("down_dataset4_kegg.csv",row.names = 1)
  GOPlotfun(goresultdown,'down_dataset4','BP') 
  keggPlotfun(keggresult1,'KEGG','down_dataset4')
}


#########up and down
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  

goresult$Type <- "UP"
goresultdown$Type <- "DOWN"

up_down <- rbind(goresult[1:10,],goresultdown)
up_down$Type <- factor(up_down$Type,levels = c("UP","DOWN"))
up_down<-up_down[order(up_down$Type),]

up_down$Description <- factor(up_down$Description,levels = rev(up_down$Description))

# up_down$Description <- factor(up_down$Description,levels =rev(up_down$Description))
allGOPlotfun <- function(goresult,picturename){
  
  pdf(paste(picturename,"bubble.pdf",sep = "_"),height=5,width=12)
  
  #resultgo$GeneRatio <- as.numeric(resultgo$GeneRatio)
  pp<-ggplot(goresult,aes(Type,Description)) + geom_point(aes(color=Type,size=Count)) +
    labs(x="",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                           panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                        color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                           axis.text.y=element_text(size=12),
                           legend.background = element_blank())+scale_colour_manual(values=c(
                             "#CC3333",  "#085A9C","#D4C009"))
  
  print(pp)
  dev.off()
  png(paste(picturename,"bubble.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
  
  
  
  pp<- ggplot(goresult,aes(x=Description,y=logPadjust,fill=Type)) + geom_bar(position=position_dodge(), stat="identity",width = 0.6) +
    scale_fill_manual(values=c("#CC3333","#085A9C","#D4C009"))+theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                                                                     panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                                                                  color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                                                                     legend.position = "right")+labs(x = " ",y="- Log10Padjust",fill="")+coord_flip()
  pdf(paste(picturename,"bar.pdf",sep = "_"),height=5,width=12)
  print(pp)
  dev.off()
  png(paste(picturename,"bar.png",sep = "_"),width = 800, height = 400)
  print(pp)
  dev.off()
} 


allGOPlotfun(up_down,"up_down_dataset4")
## module genes
##key module genes

setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  

allGOPlotfun <- function(goresult,picturename){
  #resultgo$GeneRatio <- as.numeric(resultgo$GeneRatio)
  pp<-ggplot(goresult,aes(Module,Description)) + geom_point(aes(color=Module,size=logPadjust)) +
    labs(x="",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                           panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                        color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                           axis.text.y=element_text(size=12),
                           legend.background = element_blank())+scale_colour_manual(values=c(
                             "#FF80FFFF","tan"))
  pdf(paste(picturename,"bubble.pdf",sep = "_"),height=8,width=8)
  print(pp)
  dev.off()
  png(paste(picturename,"bubble.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
  
  
  
  pp<- ggplot(goresult,aes(x=Description,y=logPadjust,fill=Module)) + geom_bar(position=position_dodge(), stat="identity",width = 0.6) +
    scale_fill_manual(values=c("#FF80FFFF","tan"))+
    theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
    panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
     color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    legend.position = "right")+labs(x = " ",y="- Log10Padjust",fill="")+coord_flip()
  pdf(paste(picturename,"bar.pdf",sep = "_"),height=8,width=13)
  print(pp)
  dev.off()
  png(paste(picturename,"bar.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
} 

pic3<-'salmon'; pic4<-'yellow';pic5='magenta';pic6='tan';pic7='brown';pic8='green';pic9='greenyellow';pic10='turquoise';pic11='purple';
pic12='red';pic13='cyan';pic14='lightcyan';pic15='pink'
backgene<-read.csv('exprs_train.csv',row.names = 1)   #72*8286

tangenes<-read.csv('tangenes_train.csv',row.names = 1) 
midnightblue<-read.csv('midnightbluegenes_train.csv',row.names = 1) 
greengenes<-read.csv('greengenes_train.csv',row.names = 1) 
#yellowgenes<-read.csv('yellowgenes_train.csv',row.names = 1)
#salmongenes<-read.csv('salmongenes_train.csv',row.names = 1)
magentagenes<-read.csv('magentagenes_train.csv',row.names = 1)
#purplegenes<-read.csv('purplegenes_train.csv',row.names = 1)
#redgenes<-read.csv('redgenes_train.csv',row.names = 1)
#bluegenes<-read.csv("bluegenes_train.csv",row.names = 1)
#greenyellow<-read.csv('greenyellowgenes_train.csv',row.names = 1)
#cyangenes<-read.csv('cyangenes.csv',row.names = 1)
#lightcyan<-read.csv('lightcyangenes.csv',row.names = 1)
#pinkgenes<-read.csv('pinkgenes.csv',row.names = 1)
#backgene2<-as.data.frame(t(read.csv("~/result1/2WGCNA/h3n2_wgcnaexprs_4143.csv",row.names = 1)))
blackgenes<-read.csv('blackgenes_train.csv',row.names = 1) 

path1<-getwd()
allgofun1(tangenes,backgene,"tan_train")
setwd(path1)
allgofun1(magentagenes,backgene,"magenta_train")
setwd(path1)
#allgofun1(purplegenes,backgene,"purple_train")
#setwd(path1)
#allgofun1(greenyellow,backgene,"greenyellow_train")
#setwd(path1)
allgofun1(midnightblue,backgene,"midightblue_train")
setwd(path1)

allgofun1(blackgenes,backgene,"black_train")
setwd(path1)
allgofun1(greengenes,backgene,"green_train")
setwd(path1)
#allgofun1(cyangenes,backgene,"cyan_train")
#setwd(path1)

tangoresult<-read.csv('tan_train/tan_train_go.csv',row.names = 1)   ##13
tangoresult<- tangoresult[tangoresult$ONTOLOGY=='BP' | tangoresult$ONTOLOGY=='MF',]
midnightbluegoresult<-read.csv('midightblue_train/midightblue_train_go.csv',row.names = 1)   ##13
midnightbluegoresult<- midnightbluegoresult[midnightbluegoresult$ONTOLOGY=='BP' | midnightbluegoresult$ONTOLOGY=='MF',]

magentagoresult<-read.csv('magenta_train/magenta_train_go.csv',row.names = 1)   ##1
magentagoresult<- magentagoresult[magentagoresult$ONTOLOGY=='BP' | magentagoresult$ONTOLOGY=='MF',]

#purplegoresult<-read.csv('purple_train/purple_train_go.csv',row.names = 1)   ##52
#purplegoresult<- purplegoresult[purplegoresult$ONTOLOGY=='BP' | purplegoresult$ONTOLOGY=='MF',]
blackgoresult<-read.csv('black_train/black_train_go.csv',row.names = 1)   ##1
blackgoresult<- blackgoresult[blackgoresult$ONTOLOGY=='BP' | blackgoresult$ONTOLOGY=='MF',]
greengoresult<-read.csv('green_train/green_train_go.csv',row.names = 1)   ##1
greengoresult<- greengoresult[greengoresult$ONTOLOGY=='BP' | greengoresult$ONTOLOGY=='MF',]


#greenyellowgoresult<-read.csv('greenyellow_train/greenyellow_train_go.csv',row.names = 1)   ##16
#greenyellowgoresult<-greenyellowgoresult[greenyellowgoresult$ONTOLOGY=='BP' |greenyellowgoresult$ONTOLOGY=='MF',]



tangoresult$Module <- "MEtan"
magentagoresult$Module <- "MEmagenta"
#purplegoresult$Module <- "MEpurple"
#greenyellowgoresult$Module <- "MEgreenyellow"

all_me <- rbind(magentagoresult,tangoresult[1:10,])#purplegoresult,,tangoresult[1:5,]magentagoresult[1:5,],
#all_me$Module <- factor(all_me$Module,levels = c("MEpurple","MEmagenta","MEgreenyellow","MEtan"))
all_me<-all_me[order(all_me$Module),]
all_me$Description<-factor(all_me$Description,levels = all_me$Description)
# all_me$Description <- factor(all_me$Description,levels =rev(all_me$Description))
allGOPlotfun(all_me,"dataset4_me")


pp<-ggplot(all_me,aes(Module,Description)) + geom_point(aes(color=Module,size=logPadjust)) +
  labs(x="",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                         panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                      color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                         axis.text.y=element_text(size=12),
                         legend.background = element_blank())+scale_colour_manual(values=c(
                          "#FF80FFFF","tan"))

pdf(paste("dataset4_me","bubble.pdf",sep = "_"),height=6,width=8)
print(pp)
dev.off()


tangoresult$Module <- "MEtan"
magentagoresult$Module <- "MEmagenta"
#redgoresult$Module <- "MEred"
#cyangoresult$Module <- "MEcyan"
blackgoresult$Module <- "MEblack"
#greenyellowgoresult$Module <- "MEgreenyellow"
greengoresult$Module <- "MEgreen"
midnightbluegoresult$Module <- "MEmidnightblue"
all_me <- rbind(magentagoresult,midnightbluegoresult,blackgoresult[1:5,],greengoresult[1:5,],tangoresult[1:5,])

#all_me <- rbind(cyangoresult[1:5,],blackgoresult[1:5,],redgoresult[1:5,],greenyellowgoresult[1:5,])
all_me$Module <- factor(all_me$Module,levels = c("MEmagenta","MEmidnightblue","MEblack","MEgreen","MEtan"))
all_me<-all_me[order(all_me$Module),]
all_me$Description<-factor(all_me$Description,levels = all_me$Description)
# all_me$Description <- factor(all_me$Description,levels =rev(all_me$Description))
allGOPlotfun1 <- function(goresult,picturename){
  #resultgo$GeneRatio <- as.numeric(resultgo$GeneRatio)
  pp<-ggplot(goresult,aes(Module,Description)) + geom_point(aes(color=Module,size=logPadjust)) +
    labs(x="",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                           panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                        color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                           axis.text.y=element_text(size=12),
                           legend.background = element_blank())+scale_colour_manual(values=c(
                             "#FF80FFFF","midnightblue","black","green","tan"))
  pdf(paste(picturename,"bubble.pdf",sep = "_"),height=6,width=8)
  print(pp)
  dev.off()
  png(paste(picturename,"bubble.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
  
  
  
  pp<- ggplot(goresult,aes(x=Description,y=logPadjust,fill=Module)) + geom_bar(position=position_dodge(), stat="identity",width = 0.6) +
    scale_fill_manual(values=c("#FF80FFFF","midnightblue","black","green","tan"))+
    theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
          panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                       color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.position = "right")+labs(x = " ",y="- Log10Padjust",fill="")+coord_flip()
  pdf(paste(picturename,"bar.pdf",sep = "_"),height=6,width=8)
  print(pp)
  dev.off()
  png(paste(picturename,"bar.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
} 
allGOPlotfun1(all_me,"dataset4_me1")

#tangoresult<-read.csv('tan_train/tan_train_kegg.csv',row.names = 1)   ##0

magentagoresult<-read.csv('magenta_train/magenta_train_kegg.csv',row.names = 1)   ##0
magentagoresult$Module <- "MEmagenta"
#greenyellowgoresult<-read.csv('greenyellow_train/greenyellow_train_kegg.csv',row.names = 1)   ##4
#greenyellowgoresult$Module <- "MEgreenyellow"
tangoresult<-read.csv('tan_train/tan_train_kegg.csv',row.names = 1)   ##0
tangoresult$Module <- "MEtan"
#purplegoresult<-read.csv('purple_train/purple_train_kegg.csv',row.names = 1)   ##0



all_me <- rbind(magentagoresult[1:10,],tangoresult)#,),purplegoresult,greenyellowgoresult
all_me<-all_me[order(all_me$Module),]
all_me$Description<-factor(all_me$Description,levels = all_me$Description)
allGOPlotfun(all_me,"dataset4_me_kegg")


####################################################333
magentagoresult<-read.csv('magenta_train/magenta_train_kegg.csv',row.names = 1)   ##0
magentagoresult$Module <- "MEmagenta"
#greenyellowgoresult<-read.csv('greenyellow_train/greenyellow_train_kegg.csv',row.names = 1)   ##4
#greenyellowgoresult$Module <- "MEgreenyellow"
greengoresult<-read.csv('green_train/green_train_kegg.csv',row.names = 1)   ##4
greengoresult$Module <- "MEgreen"
tangoresult<-read.csv('tan_train/tan_train_kegg.csv',row.names = 1)   ##0
tangoresult$Module<-"MEtan"
midnightbluegoresult<-read.csv('midightblue_train/midightblue_train_kegg.csv',row.names = 1)   ##0
midnightbluegoresult$Module<-"MEmidnightblue"

#redgoresult<-read.csv('red_train/red_train_kegg.csv',row.names = 1)   ##0
#redgoresult$Module<-"MEred"
#cyangoresult<-read.csv('cyan_train/cyan_train_kegg.csv',row.names = 1)   ##0
#cyangoresult$Module<-"MEcyan"
blackgoresult<-read.csv('black_train/black_train_kegg.csv',row.names = 1)   ##0
blackgoresult$Module<-"MEblack"

all_me <- rbind(magentagoresult[1:5,],midnightbluegoresult[1:5,],blackgoresult[1:5,],greengoresult,tangoresult)#,)tangoresult,purplegoresult
all_me$Module <- factor(all_me$Module,levels = c("MEmagenta","MEmidnightblue","MEblack","MEgreen","MEtan"))

all_me<-all_me[order(all_me$Module),]
all_me$Description<-factor(all_me$Description,levels = all_me$Description)
allGOPlotfun <- function(goresult,picturename){
  #resultgo$GeneRatio <- as.numeric(resultgo$GeneRatio)
  pp<-ggplot(goresult,aes(Module,Description)) + geom_point(aes(color=Module,size=logPadjust)) +
    labs(x="",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                           panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                        color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                           axis.text.y=element_text(size=12),
                           legend.background = element_blank())+scale_colour_manual(values=c(
                             "#FF80FFFF","midnightblue","black","green","tan"))
  pdf(paste(picturename,"bubble.pdf",sep = "_"),height=7,width=7)
  print(pp)
  dev.off()
  png(paste(picturename,"bubble.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
  
  
  
  pp<- ggplot(goresult,aes(x=Description,y=logPadjust,fill=Module)) + geom_bar(position=position_dodge(), stat="identity",width = 0.6) +
    scale_fill_manual(values=c("#FF80FFFF","midnightblue","black","green","tan"))+
    theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
          panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                       color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.position = "right")+labs(x = " ",y="- Log10Padjust",fill="")+coord_flip()
  pdf(paste(picturename,"bar.pdf",sep = "_"),height=8,width=8)
  print(pp)
  dev.off()
  png(paste(picturename,"bar.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
}
allGOPlotfun(all_me,"dataset4_me_kegg1")




