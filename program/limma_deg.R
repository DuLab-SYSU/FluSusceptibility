library(limma)
exprs_h3n2<-read.csv('D:/wrok/flu_paper/data/result2/result_new/h3n2_train.csv',row.names = 1)   #72*8286
#exprs<-exprs_h3n2[1:50,]
exprs.gnames<-rownames(exprs_h3n2)

##Group by outcome
exprs.cl<-c(1 ,1, 1, 1, 1, 1, 1, 1,0, 0, 0, 0,0,0,1, 1, 1, 1 ,1 ,0 ,0, 0 ,0 ,1 ,1, 1 ,1 ,1 ,1 ,1, 1, 0, 0, 0, 0, 0 ,0, 1, 1 ,1 ,1, 1 ,1,1, 0 ,0 ,0 ,0 ,0,0,0,0, 1, 1, 1  ,0 ,0, 0, 1 ,1 ,1, 1, 1 ,1  ,0 ,0, 0 ,0, 0 ,0,0,0)
##Group by batch
exprs.origin<-c(1, 1 ,1, 1, 1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,2 ,2 ,2, 2 ,2 ,2, 2 ,2, 2 ,3, 3, 3, 3 ,3 ,3, 3 ,3 ,3 ,3, 3, 3, 3, 3 ,4 ,4 ,4, 4 ,4 ,4 ,4 ,4 ,4, 4, 4 ,4 ,4, 4 ,4, 5 ,5 ,5 ,5 ,5 ,5 ,6, 6 ,6, 6 ,6 ,6 ,6 ,6, 6, 6, 6,6,6, 6)
exprs_h3n2<-as.matrix(exprs_h3n2) 
exp=exprs_h3n2
coldata<-data.frame(row.name=colnames(exp),condition=exprs.cl)
coldata$condition[coldata$condition==1]<-"sym"
coldata$condition[coldata$condition==0]<-"asym"
coldata$condition<- factor(coldata$condition)
design <- model.matrix(~0+coldata$condition)
colnames(design)=levels(coldata$condition)
rownames(design)=colnames(exp)
constrasts = paste0(unique(coldata$condition),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit <- lmFit(exprs_h3n2, design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
#logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
#logFC_cutoff=1
#k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
#k2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(DEG$logFC>0&DEG$adj.P.Val<=0.05,"UP",ifelse(DEG$logFC<0&DEG$adj.P.Val<=0.05,"DOWN","NOT"))
table(DEG$change)
head(DEG)
limma_voom_DEG <- DEG
write.csv(limma_voom_DEG,file = "D:/wrok/flu_paper/data/result2/result_new/limma_voom_DEG.csv")

limma_deg <- subset(limma_voom_DEG,change=="UP"|change=="DOWN")
deg_expression <- exprs_h3n2[rownames(limma_deg),]
write.csv(deg_expression,"deg_expression_limma.csv")
limma_up <- subset(limma_voom_DEG,change=="UP")
limma_down <- subset(limma_voom_DEG,change=="DOWN")


length(intersect(rownames(degs),rownames(limma_deg)))
#243
length(intersect(rownames(up_deg),rownames(limma_deg[limma_deg$change=="UP",])))
#154
length(intersect(rownames(up_deg),rownames(limma_deg[limma_deg$change=="DOWN",])))
#0
length(intersect(rownames(down_deg),rownames(limma_deg[limma_deg$change=="DOWN",])))
#89

up_deg_20<-limma_up[order(limma_up$logFC,decreasing = TRUE),][1:20,]
down_deg_20<-limma_down[order(limma_down$logFC,decreasing = TRUE),][1:20,]
up_deg_20<-exprs_h3n2[rownames(up_deg_20),]
down_deg_20<-exprs_h3n2[rownames(down_deg_20),]
deg_heatmap<-rbind(up_deg_20,down_deg_20)

deg_heatmap1<-cbind(deg_heatmap[,exprs.cl=='1'],deg_heatmap[,exprs.cl=='0'])
annotation_col<-as.factor(c(rep('sym',37),rep('asym',35)))
annotation_col<-data.frame(Traittype=annotation_col)
rownames(annotation_col)<-colnames(deg_heatmap1)
gene_anno<-data.frame(Gene_type=factor(c(rep('up',20),rep('down',20))))
rownames(gene_anno)<-rownames(deg_heatmap)
#ggplot(deg_heatmap,aes())
pheatmap(log(deg_heatmap1+1),annotation_col = annotation_col,show_colnames = FALSE,#annotation_row = gene_anno,
         cluster_cols = FALSE,cluster_rows = T,scale = 'row',#gaps_row = 20,
         annotation_names_row = FALSE,annotation_names_col = FALSE,color = colorRampPalette(c( "#085A9C", "white","#CC3333"))(8),#border=FALSE,
         annotation_colors = list(Traittype=c(asym="#085A9C",sym="#CC3333")))

###########go and kegg

#backgene<-read.csv('D:/wrok/flu_paper/data/result2/result_new/h3n2_train.csv',row.names = 1)   #72*8286
pic1<-'up_limma'
pic2<-'down_limma'
path1<-getwd()
  allgofun(limma_up,exprs_h3n2,pic1)
  goresult<-read.csv("up_limma_go.csv",row.names = 1)   ##24
  goresult<-goresult[goresult$ONTOLOGY=='BP',]   #19
  keggresult<-read.csv("up_limma_kegg.csv",row.names = 1)  ##0
  
  GOPlotfun(goresult,'up','BP')
  keggPlotfun(keggresult,'KEGG','up-2')
  
  ##rankprod down
  setwd('D:/wrok/flu_paper/data/result2/result_new/')
  path1<-getwd()
  allgofun(limma_down,exprs_h3n2,pic2)
  goresultdown<-read.csv("down_limma_go.csv",row.names = 1)   ###32
  goresultdown<-goresultdown[goresultdown$ONTOLOGY=='BP',][1:15,]
  keggresult1<-read.csv("down_limma_kegg.csv",row.names = 1)
  GOPlotfun(goresultdown,'down_limma','BP') 
  keggPlotfun(keggresult1,'KEGG_limma','down')



#########up and down
setwd('D:/wrok/flu_paper/data/result2/result_new/')

goresult$Type <- "UP"
goresultdown$Type <- "DOWN"

up_down <- rbind(goresult[1:5,],goresultdown[1:5,])

#up_down$Type <- factor(up_down$Type,levels = c("UP","DOWN"))
# up_down$Description <- factor(up_down$Description,levels =rev(up_down$Description))
allGOPlotfun(up_down,"up_down")






###########
purplegenes<-read.csv("D:/wrok/flu_paper/data/result2/result_new/purple_hub.csv",row.names = 1)
tangenes<-read.csv("tan_hub.csv",row.names = 1)
redgenes<-read.csv('red_hub.csv',row.names = 1)
browngenes<-read.csv('brown_hub.csv',row.names = 1)
cyangenes<-read.csv('cyan_hub.csv',row.names = 1)
hubgene0<-c(rownames(purplegenes),rownames(tangenes),rownames(redgenes),rownames(browngenes),rownames(cyangenes))  ##a total of 28
#write.table(hubgene0,'hubgene_wgcna.txt')


venn.plot <- venn.diagram(
  x = list( "Hub genes" = hubgene0,"Limma_DEGs" = rownames(limma_deg)),
  filename = NULL,
  height = 450, 
  width = 450,
  resolution =300, 
  #col = "transparent",      
  fill = c("#008B4533", "#63187933"),  #"#008B4533", "#63187933"'#008B4599'
  #alpha = 0.50,                                      
  #label.col = c("orange", "white", "darkorchid4", "white",
  #             "white", "darkgreen", "white"),
  cex = 2,    
  fontfamily = "serif",  
  fontface = "bold",    
  cat.col = c('black','black'),  #'#EE000099','#008B4599'
  cat.cex = 2,      
  cat.pos = c(100, 260),#0,        
  cat.dist = c(0.07, 0.07),#, 0.05    
  cat.fontfamily = "serif",     
  rotation.degree =180,        
  margin = 0.2               
)

grid.draw(venn.plot);

dev.off()


venn.plot <- venn.diagram(
  x = list( "RankPord" = rownames(degs),"Limma" = rownames(limma_deg)),
  filename = NULL,
  height = 450, 
  width = 450,
  resolution =300, 
  #col = "transparent",      
  fill = c("#008B4533", "#63187933"),  #"#008B4533", "#63187933"'#008B4599'
  #alpha = 0.50,                                      
  #label.col = c("orange", "white", "darkorchid4", "white",
  #             "white", "darkgreen", "white"),
  cex = 2,    
  fontfamily = "serif",  
  fontface = "bold",    
  cat.col = c('black','black'),  #'#EE000099','#008B4599'
  cat.cex = 2,      
  cat.pos = c(100, 260),#0,        
  cat.dist = c(0.07, 0.07),#, 0.05    
  cat.fontfamily = "serif",     
  rotation.degree =180,        
  margin = 0.2               
)

grid.draw(venn.plot);

dev.off()




venn.plot <- venn.diagram(
  x = list( "UP_RankPord" = rownames(up_deg),"DOWN_RankPord" = rownames(down_deg),"UP_Limma" = rownames(limma_up),"DOWN_Limma" = rownames(limma_down)),
  filename = NULL,
  height = 450, 
  width = 450,
  resolution =300, 
  #col = "transparent",      
  fill = c("#008B4533", "#63187933",'#008B4599','#EE000099'),  #"#008B4533", "#63187933"'#008B4599'
  #alpha = 0.50,                                      
  #label.col = c("orange", "white", "darkorchid4", "white",
  #             "white", "darkgreen", "white"),
  cex = 2,    
  fontfamily = "serif",  
  fontface = "bold",    
  cat.col = c('black','black','black','black'),  #'#EE000099','#008B4599'
  cat.cex = 2,      
  #cat.pos = c(100, 260,0,0),#0,        
  #cat.dist = c(0.07, 0.07,0.05,0),#, 0.05    
  cat.fontfamily = "serif",     
  #rotation.degree =180,        
  margin = 0.2               
)

grid.draw(venn.plot);

dev.off()
