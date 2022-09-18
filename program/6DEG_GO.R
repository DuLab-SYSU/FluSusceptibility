library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(topGO)
library(reshape)
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
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  

#setwd("~/Documents/Computational_Biology/Analysis/H3N2_symptomatic_DEG/")
### Data (2021-08-09)
## pfp < 0.05
degs <- read.csv("deg_all.csv", header = TRUE)
degs2 <- dplyr::filter(degs, type != "Baseline") %>%
  mutate(type = str_replace(type, "(Hour) (\\d+)", "\\2h")) %>%
  dplyr::select(degtype:gene) %>%
  dplyr::rename(time = type) %>%
  dplyr::filter(time %in% c("00h","021h","045h", "069h", "093h", "108h"))

#("22h", "46h", "70h", "94h", "118h", "142h", "166h")

degs3 <- degs2 %>%
  group_by(degtype, time, group) %>%
  dplyr::summarise(number = n())
degs3 <- unite(degs3, "type", group, degtype, sep = "_", remove = FALSE)

p0.05 <- ggplot(degs3, aes(time, number)) +
  geom_point(aes(color = degtype, shape = group), size = 4) +
  geom_line(aes(color = degtype, group = type, linetype = group), size = 1) +
  scale_x_discrete(limits = stringr::str_sort(unique(degs3$time), 
                                              numeric = TRUE)) +
  scale_color_manual(values = brewer.pal(9,"Set1")[c(2,1)]) +
  labs(x = "Time point", y = "Number", title = "pfp < 0.05") +
  theme_classic() +
  mytheme2 +
  theme(legend.position = c(.85, .85))

## pfp < 0.01
degs0.01 <- dplyr::filter(degs, type != "Baseline" & pfp < 0.01) %>%
  mutate(type = str_replace(type, "(Hour) (\\d+)", "\\2h")) %>%
  dplyr::select(degtype:gene) %>%
  dplyr::rename(time = type) %>%
  dplyr::filter(time %in% c("00h","021h","045h", "069h", "093h", "108h"))#("012h","036h", "060h", "084h", "108h

degs0.01_2 <- degs0.01 %>%
  group_by(degtype, time, group) %>%
  dplyr::summarise(number = n())
degs0.01_2 <- unite(degs0.01_2, "type", group, degtype, sep = "_", remove = FALSE)
degs0.01_2$degtype<-factor(degs0.01_2$degtype,levels=c("up","down"))
p0.01 <- ggplot(degs0.01_2, aes(time, number)) +
  geom_point(aes(color = group, shape = degtype), size = 4) +
  geom_line(aes(color = group, group = type), size = 1) +
  scale_x_discrete(limits = stringr::str_sort(unique(degs0.01_2$time), 
                                              numeric = TRUE)) +
  scale_color_manual(values=c("#085A9C","#EF0808")) +
  labs(x = "Time point", y = "Number", title = "pfp < 0.01") +
  #theme_classic() +
  mytheme2 +
  theme(legend.position = c(.15, .80))+facet_grid(~degtype)
pdf(paste("all_time_number",".pdf",sep = "_"),height=4,width=6)
print(p0.01)
dev.off()
### GO enrichment (2021-08-18)
degs4 <- aggregate(degs0.01[4], degs0.01[1:3],
                   FUN = function(X) paste(unique(X), collapse=","))
degs4 <- unite(degs4, "type", group, degtype, time, sep = "_", remove = FALSE)
exprs_h3n2<-read.csv('gse30550_sym_shedding_exprs_alltime.csv',row.names = 1)

backgene <-rownames(exprs_h3n2)
all_term <- data.frame(ONTOLOGY=character(),ID=character(), Description = character(), 
                       GeneRatio = character(), BgRatio = character(), pvalue = numeric(), 
                       p.adjust = numeric(), qvalue=numeric(),geneID = character(),
                       Count = numeric(), logPadjust=numeric(),Type = character())

all_kegg <- data.frame(ID=character(), Description = character(), 
                       GeneRatio = character(), BgRatio = character(), pvalue = numeric(), 
                       p.adjust = numeric(), qvalue=numeric(),geneID = character(),
                       Count = numeric(), logPadjust=numeric(),Type = character())
for (i in 1:nrow(degs4)){
  #if (length(intersect(str_split(degs4[i,5], ",")[[1]], rownames(exprs_h3n2))) > 1){
    gene1<-str_split(degs4[i,5], ",")[[1]]
    pic<-degs4$type[i]
  tmp <- allgofun(gene1,backgene,pic)
  tmp2 <- allkeggfun(gene1,backgene,pic)
  if(nrow(tmp)>=1){
    all_term <-rbind(all_term,tmp)
    all_kegg<-rbind(all_kegg,tmp2)
  }
#}
}
write.csv(all_term,"all_time_degterm.csv")
write.csv(all_kegg,"all_time_degkeggterm.csv")

all_term<-read.csv("all_time_degterm.csv",row.names = 1)
all_term$group <- apply(all_term,1,function(x) strsplit(x[12],split = "_")[[1]][1])
all_term$deg <- apply(all_term,1,function(x) strsplit(x[12],split = "_")[[1]][2])
all_term$deg <-factor(all_term$deg,levels = c("up","down"))
all_term$time <- apply(all_term,1,function(x) strsplit(x[12],split = "_")[[1]][3])
all_term$time <- factor(all_term$time,levels = c("021h","045h", "069h", "093h", "108h"))
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

  pdf(paste(picturename,"bubble.pdf",sep = "_"),height=9,width=12)
  
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
  png(paste(picturename,"bubble.png",sep = "_"),width = 800, height = 800)
  print(pp)
  dev.off()
} 
allGOPlotfun(all_term1,"all_time_go")



all_term<-read.csv("all_time_degkeggterm.csv",row.names = 1)
all_term$group <- apply(all_term,1,function(x) strsplit(x[11],split = "_")[[1]][1])
all_term$deg <- apply(all_term,1,function(x) strsplit(x[11],split = "_")[[1]][2])
all_term$deg <-factor(all_term$deg,levels = c("up","down"))
all_term$time <- apply(all_term,1,function(x) strsplit(x[11],split = "_")[[1]][3])
all_term$time <- factor(all_term$time,levels = c("021h","045h", "069h", "093h", "108h"))
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
  
  pdf(paste(picturename,"kegg_bubble.pdf",sep = "_"),height=8,width=8)
  
  
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
  png(paste(picturename,"kegg_bubble.png",sep = "_"),width = 800, height = 800)
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
  png(paste(picturename,"kegg_bar.png",sep = "_"),width = 800, height = 400)
  print(pp)
  dev.off()
}

keggPlotfun(all_term1,"all_time")
