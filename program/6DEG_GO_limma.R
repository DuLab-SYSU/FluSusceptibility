library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(topGO)
library(reshape)
library(ReactomePA)

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
mytheme1 <- theme(panel.grid.major = element_line(colour=brewer.pal(9,"Pastel1")[9],
                                                 linetype = "longdash"),
                 panel.background = element_rect(fill='transparent', color="#000000"),
                 panel.border=element_rect(fill='transparent', color='black'),
                 axis.text = element_text(size = 16),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 13),
                 legend.title = element_text(size = 16),
                 legend.background = element_blank())
mytheme2 <- theme(
                  axis.text = element_text(size = 16),
                  axis.title = element_text(size = 20),
                  legend.text = element_text(size = 13),
                  legend.title = element_blank(),
                  title = element_text(size = 20))+theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                                                         panel.background = element_blank(),panel.border=element_rect(fill='transparent',
                                                                                                                      color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                                                         axis.text.y=element_text(size=12),
                                                         legend.background = element_blank())
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  
allgofun <- function(intregene,backgene,picturename){
  #dir.create(picturename)
  #path=paste(path1,picturename,sep = "/")
  #setwd(path)
  intregenes <- bitr(intregene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
  backgenes <- bitr(backgene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
  ##GO
  goann <- GOfun(intregenes,backgenes)
  if(nrow(goann)>0){
    goann$type <-picturename
  }else{
    print(paste(picturename,"NO ENRICHMENT"))
  }
  write.csv(goann,paste(picturename,"go.csv",sep = "_"))
  ##KEGG
  #keggann <- keggfun(intregenes,backgenes)
  #write.csv(keggann,paste(picturename,"kegg.csv",sep = "_"))
  #keggpp <- keggPlotfun(keggann,picturename,"Kegg")
  return(goann)}


allkeggfun <- function(intregene,backgene,picturename){
  #dir.create(picturename)
  #path=paste(path1,picturename,sep = "/")
  #setwd(path)
  intregenes <- bitr(intregene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
  backgenes <- bitr(backgene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = T)
  ##GO
  keggann <- keggfun(intregenes,backgenes)
  if(nrow(keggann)>0){
    keggann$type <-picturename
  }else{
    print(paste(picturename,"NO ENRICHMENT"))
  }
  write.csv(keggann,paste(picturename,"kegg.csv",sep = "_"))
  ##KEGG
  #keggann <- keggfun(intregenes,backgenes)
  #write.csv(keggann,paste(picturename,"kegg.csv",sep = "_"))
  #keggpp <- keggPlotfun(keggann,picturename,"Kegg")
  return(keggann)}##

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

#setwd("~/Documents/Computational_Biology/Analysis/H3N2_symptomatic_DEG/")
### Data (2021-08-09)
## pfp < 0.05
degs <- read.csv("deg_limma_alltime_summary.csv", header = TRUE,row.names = 1)
degs2 <- dplyr::filter(degs, degtype != "NotSig")%>%
                         dplyr::rename(time = type) 
#degs2<- subset(degs2,degtype=="DOWN"|degtype=="UP")
#%>%
 # dplyr::filter(time %in% c("00h","021h","045h", "069h", "093h", "108h"))

#("22h", "46h", "70h", "94h", "118h", "142h", "166h")

degs3 <- degs2 %>%
  group_by(degtype, time, group) %>%
  dplyr::summarise(number = n())
degs3 <- unite(degs3, "type", group, degtype, sep = "_", remove = FALSE)
degs3$degtype<-factor(degs3$degtype,levels=c("UP","DOWN"))

p0.05 <- ggplot(degs3, aes(time, number)) +
  geom_point(aes(color = group, shape = degtype), size = 4) +
  geom_line(aes(color = group, group = type), size = 1) +
  scale_x_discrete(limits = stringr::str_sort(unique(degs0.01_2$time), 
                                              numeric = TRUE)) +
  scale_color_manual(values=c("#085A9C","#EF0808")) +
  labs(x = "Time point", y = "Number", title = "BH<= 0.05") +
  #theme_classic() +
  mytheme2 +
  theme(legend.position = c(.15, .80))+facet_grid(~degtype)
pdf(paste("all_time_number_0.05",".pdf",sep = "_"),height=4,width=6)
print(p0.05)
dev.off()
## pfp < 0.01
degs0.01 <- dplyr::filter(degs2,padj_BH <= 0.01)

degs0.01_2 <- degs0.01 %>%
  group_by(degtype, time, group) %>%
  dplyr::summarise(number = n())
degs0.01_2 <- unite(degs0.01_2, "type", group, degtype, sep = "_", remove = FALSE)
degs0.01_2$degtype<-factor(degs0.01_2$degtype,levels=c("UP","DOWN"))
p0.01 <- ggplot(degs0.01_2, aes(time, number)) +
  geom_point(aes(color = group, shape = degtype), size = 4) +
  geom_line(aes(color = group, group = type), size = 1) +
  scale_x_discrete(limits = stringr::str_sort(unique(degs0.01_2$time), 
                                              numeric = TRUE)) +
  scale_color_manual(values=c("#085A9C","#EF0808")) +
  labs(x = "Time point", y = "Number", title = "BH<= 0.01") +
  #theme_classic() +
  mytheme2 +
  theme(legend.position = c(.15, .80))+facet_grid(~degtype)
pdf(paste("all_time_number",".pdf",sep = "_"),height=4,width=6)
print(p0.01)
dev.off()
### GO enrichment (2021-08-18)
degs4 <- degs0.01[,c(7,1:3)]
#aggregate(degs0.01[7], degs0.01[1:3],
 #                  FUN = function(X) paste(unique(X), collapse=","))
degs4 <- unite(degs4, "type", group, degtype, time, sep = "_", remove = FALSE)
exprs_h3n2<-read.csv('dee2_sym_exprs_alltime.csv',row.names = 1)

backgene <-rownames(exprs_h3n2)
all_term <- data.frame(ONTOLOGY=character(),ID=character(), Description = character(), 
                       GeneRatio = character(), BgRatio = character(), pvalue = numeric(), 
                       p.adjust = numeric(), qvalue=numeric(),geneID = character(),
                       Count = numeric(), logPadjust=numeric(),Type = character())

all_kegg <- data.frame(ID=character(), Description = character(), 
                       GeneRatio = character(), BgRatio = character(), pvalue = numeric(), 
                       p.adjust = numeric(), qvalue=numeric(),geneID = character(),
                       Count = numeric(), logPadjust=numeric(),Type = character())
for (i in unique(degs4$type)){
  #if (length(intersect(str_split(degs4[i,5], ",")[[1]], rownames(exprs_h3n2))) > 1){
  tt <- subset(degs4,type==i)  
  gene1<-tt$gene
    pic<-tt$type[1]
  tmp <- allgofun(gene1,backgene,pic)
  tmp2 <- allkeggfun(gene1,backgene,pic)
  if(nrow(tmp)>=1){
    all_term <-rbind(all_term,tmp)
    all_kegg<-rbind(all_kegg,tmp2)
  }
#}
}
write.csv(all_term,"all_time_degterm_limma.csv")
write.csv(all_kegg,"all_time_degkeggterm_limma.csv")

all_term<-read.csv("all_time_degterm_limma.csv",row.names = 1)
all_term$group <- apply(all_term,1,function(x) strsplit(x[12],split = "_")[[1]][1])
all_term$deg <- apply(all_term,1,function(x) strsplit(x[12],split = "_")[[1]][2])
all_term$deg <-factor(all_term$deg,levels = c("UP","DOWN"))
all_term$time <- apply(all_term,1,function(x) strsplit(x[12],split = "_")[[1]][3])
all_term$time <- factor(all_term$time,levels = c( "h36","h46","h53","h60","h70","h77"
                                                  ,"h84","h94", "h101", "h108", "h118", "h125", "h132", "h142", "h166"))
all_term1<-all_term %>% arrange(desc(logPadjust)) %>% group_by(ONTOLOGY,type) %>% filter(row_number() <= 5)

all_term1 <- subset(all_term1,ONTOLOGY=="BP")
all_term1$nlength <- apply(all_term1,1,function(x)nchar(x[3]))
all_term1<-all_term1[order(all_term1$nlength,decreasing = F),]
all_term1<-as.data.frame(all_term1)
shunxu <-unique(all_term1$Description)
shunxu <-factor(shunxu,levels = shunxu[order(apply(data.frame(shunxu),1,nchar))])
all_term1$Description<-factor(all_term1$Description,levels =c(shunxu))

allGOPlotfun <- function(goresult,picturename){
  #shunxu <-unique(goresult$Description)
  #shunxu <-factor(shunxu,levels = shunxu[order(apply(data.frame(shunxu),1,nchar))])
  goresult<-goresult[order(goresult$time),]
  goresult$Description<-factor(goresult$Description,levels =unique(goresult$Description))

  pdf(paste(picturename,"bubble.pdf",sep = "_"),height=12,width=11)
  
  #resultgo$GeneRatio <- as.numeric(resultgo$GeneRatio)
  pp<-ggplot(goresult,aes(time,Description)) + geom_point(aes(color=group,size=logPadjust)) +
    #scale_y_discrete(limits = rev(shunxu)) +
    labs(x="",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                           panel.background = element_blank(),panel.border=element_rect(fill='transparent',
                        color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                           axis.text.y=element_text(size=12),
                           legend.background = element_blank())+scale_colour_manual(values=c(
                             "#085A9C","#EF0808"))+facet_grid(~deg)
  
  print(pp)
  dev.off()
} 
allGOPlotfun(all_term1,"all_time_go")



all_term<-read.csv("all_time_degkeggterm_limma.csv",row.names = 1)
all_term$group <- apply(all_term,1,function(x) strsplit(x[11],split = "_")[[1]][1])
all_term$deg <- apply(all_term,1,function(x) strsplit(x[11],split = "_")[[1]][2])
all_term$deg <-factor(all_term$deg,levels = c("UP","DOWN"))
all_term$time <- apply(all_term,1,function(x) strsplit(x[11],split = "_")[[1]][3])
all_term$time <- factor(all_term$time,levels = c( "h36","h46","h53","h60","h70","h77"
                                                  ,"h84","h94", "h101", "h108", "h118", "h125", "h132",  "h166"))#"h142",
all_term1<-all_term %>% arrange(desc(logPadjust)) %>% group_by(type) %>% filter(row_number() <= 5)

#all_term1 <- subset(all_term1,ONTOLOGY=="BP")
all_term1$nlength <- apply(all_term1,1,function(x)nchar(x[3]))
all_term1<-all_term1[order(all_term1$nlength,decreasing = F),]
all_term1<-as.data.frame(all_term1)
shunxu <-unique(all_term1$Description)
shunxu <-factor(shunxu,levels = shunxu[order(apply(data.frame(shunxu),1,nchar))])
all_term1$Description<-factor(all_term1$Description,levels =c(shunxu))

keggPlotfun <- function(goresult,picturename){
  #goresult <- ifelse(nrow(goresult)<=20,goresult,goresult[1:20,])
  #goresult<-term_numfun(goresult)
  goresult<-goresult[order(goresult$time),]
  goresult$Description<-factor(goresult$Description,levels =unique(goresult$Description))
  
  pdf(paste(picturename,"kegg_bubble.pdf",sep = "_"),height=8,width=11)
  
  
  pp<-ggplot(goresult,aes(time,Description)) + geom_point(aes(color=group,size=logPadjust)) +
    #scale_y_discrete(limits = rev(shunxu)) +
    labs(x="",y="") +theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                           panel.background = element_blank(),panel.border=element_rect(fill='transparent',
                                                                                        color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=12,angle = 90),
                           axis.text.y=element_text(size=12),
                           legend.background = element_blank())+scale_colour_manual(values=c(
                             "#085A9C","#EF0808"))+facet_grid(~deg)
  
  print(pp)
  dev.off()
  
  pp<- ggplot(goresult,aes(x=Description,y=logPadjust,fill="group")) + geom_bar(position=position_dodge(), stat="identity",width = 0.45) +
    scale_fill_manual(values=c( "#085A9C","#EF0808"))+theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
                                                 panel.background = element_blank(),panel.border=element_rect(fill='transparent', 
                                                                                                              color='black'), axis.ticks.x = element_blank(),axis.text.x=element_text(size=18),
                                                 axis.text.y=element_text(size=18),
                                                 legend.position = "none")+labs(x = " ",y="- Log10Padjust",fill="")+coord_flip()
  pdf(paste(picturename,"kegg_bar.pdf",sep = "_"),height=8,width=8)
  print(pp)
  dev.off()
}

keggPlotfun(all_term1,"all_time")
