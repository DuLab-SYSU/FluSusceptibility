#===============================[ WGCNA for top SD  genes ]=====================================
rm(list=ls())
gc()
options(stringsAsFactors = FALSE)
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  
#===================================[ load pkg and data ]======================================

  library(WGCNA)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(futile.logger)
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
  ###all data
  exprs_h3n2<-read.csv('h3n2_combat_data.csv',row.names = 1)
  h3n2_group <- read.csv("group_h3n2.csv",row.names = 1)
  # h3n2_group<-c(1 ,1, 1, 1, 1, 1, 1, 1,0, 0, 0, 0,0,0,1, 1, 1, 1 ,1 ,0 ,0, 0 ,0 ,1 ,1, 1 ,1 ,1 ,1 ,1, 1, 0, 0, 0, 0, 0 ,0, 1, 1 ,1 ,1, 1 ,1,1, 0 ,0 ,0 ,0 ,0,0,0,0, 1, 1, 1  ,0 ,0, 0, 1 ,1 ,1, 1, 1 ,1  ,0 ,0, 0 ,0, 0 ,0,0,0)

# exprs_h3n2<-exprs_h3n2[!(rownames(exprs_h3n2)%in%c("EIF1AY","KDM5D","DDX3Y","PRS4Y1","PRKY")),]
#write.csv(exprs_h3n2,"exprs_h3n2_data_nosexgene.csv")
#==========================================[ data prepare ]===========================================

  #use the top sd gene 
  exprs<-apply(exprs_h3n2,1,var)
  #exprs<-as.data.frame(t(exprs_h3n2[which(exprs>quantile(exprs, probs = seq(0, 1, 0.25))[4]),]))  ##2027*72
  exprs<-as.data.frame(t(exprs_h3n2[which(exprs>quantile(exprs, probs = seq(0, 1, 0.2))[5]),]))  ##1696*49
  
  #exprs<-as.data.frame(t(exprs_h3n2[which(exprs>quantile(exprs, probs = seq(0, 1, 0.15))[10]),]))  ##1696*49
  #exprs<-as.data.frame(t(exprs_h3n2[which(exprs>quantile(exprs, probs = 0.85)),]))  ##1696*49
  
  
  #check the outlier sample
  Traits <- subset(h3n2_group,select=c("gseID","group"))
  Traits<- Traits[rownames(exprs),]
  #  Traits<- data.frame(sample=rownames(exprs),type=h3n2_group)
  sampletree<-hclust(dist(exprs),method = 'average')
  traitColors = numbers2colors(as.numeric(as.factor(Traits$group)), signed = TRUE,centered=TRUE)
  plotDendroAndColors(sampletree, traitColors,groupLabels = Traits$type,cex.dendroLabels = 0.8,marAll = c(1,4,3,1),cex.rowText = 0.02,
                      main = "Sample dendrogram and trait heatmap")

  write.csv(exprs,'h3n2_wgcnaexprs_1696_all.csv')
  
#===================================[ build network ]===============================
    powers<-c(c(1:10), seq(from = 12, to=20, by=2))
    #powers<-c(1:20)
    
    sft = pickSoftThreshold(exprs, powerVector = powers, RsquaredCut = 0.9,verbose = 5,)  ###9
    sizeGrWindow(9, 5)
    par(mfrow = c(1,2))
    cex1 = 0.90
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red")
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.9,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  
  ####One-step building blocks:14 modules
    cor <- WGCNA::cor
    net = blockwiseModules(exprs, power =6#sft$powerEstimate
                           ,maxBlockSize = 5000, TOMType = "unsigned", minModuleSize = 15,reassignThreshold = 0, 
                           mergeCutHeight = 0.2,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,
                           saveTOMFileBase = "TOM",verbose = 5,deepSplit = 4)#sft$powerEstimate 
    table(net$colors) ##minModuleSize=20:12;15:16
    cor<-stats::cor
    mergedColors = labels2colors(net$colors)
    write.csv(table(mergedColors),'mergedColors.csv')
    number_module <- data.frame(module=names(table(mergedColors)),number=table(mergedColors))
    number_module$module_color <- sort(unique(mergedColors[net$blockGenes[[1]]]))
    
   pp<- ggplot(number_module,aes(x=module,y=number.Freq,fill=module))+geom_bar(stat = "identity",width = 0.8
                                                                           , alpha = 0.8)+
      scale_fill_manual(values = number_module$module_color)+mytheme+
      geom_text(aes(label = number.Freq),color="#E6E8FA", position = position_dodge(0.1))+
      theme(axis.text.x=element_text(size=12,angle = 90)) 
   pdf("module_number_all",height=4,width=9)
   print(pp)
   dev.off()
    gene_color <-mergedColors[net$blockGenes[[1]]]
    
    plotDendroAndColors(net$dendrograms[[1]],gene_color, "Module colors",
                        dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
    ###save variable
    moduleLabels = net$colors
    moduleColors = labels2colors(net$colors)
    MEs = net$MEs
    geneTree = net$dendrograms[[1]]
    save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna_median-networkConstruction-auto.RData")

#=======================================[ Screening modules related to phenotype ]===================================
    
  nGenes = ncol(exprs)
  nSamples = nrow(exprs)
  design=model.matrix(~0+ as.factor( Traits$group))
  colnames(design)=c('asym','sym')
  ##The first principal component of each module
  MEs0 = moduleEigengenes(exprs, moduleColors)$eigengenes
  MEs=orderMEs(MEs0)
  moduleTraitCor = cor(MEs, design,use = "p",method = "spearman")
  moduleTraitCor1<-moduleTraitCor[order(moduleTraitCor[,"sym"],decreasing = TRUE),]
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor1, nSamples)
  textMatrix = paste(signif(moduleTraitCor1, 2), "\n(",signif(moduleTraitPvalue, 2), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor1)
  sizeGrWindow(10,18)
  par(mar = c(4, 10, 3, 2))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor1,
                 xLabels = colnames(design),
                 yLabels = rownames(moduleTraitCor1),
                 ySymbols = rownames(moduleTraitCor1),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(100),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  MET = orderMEs(cbind(MEs, design))
  #################
  library(reshape)
  library(reshape2)
  cor<-as.data.frame(moduleTraitCor1)
  cor$Module <- rownames(cor)
  cor_melt <-melt(cor,id="Module")
  colnames(cor_melt)<-c("Module","group","R2")
  
  pp<-as.data.frame(moduleTraitPvalue)
  pp$Module <- rownames(pp)
  pp_melt <-melt(pp,id="Module")
  colnames(pp_melt)<-c("Module","group","Pvalue")
  cor_pp <-cbind(cor_melt,pp_melt)
  cor_pp<-cor_pp[,-c(4,5)]
  cor_pp<-cor_pp[order(cor_pp$Module),]
  cor_pp$mm_color <- rep(number_module$module,each=2) 
  ggplot(cor_pp,aes(R2,-log10(Pvalue),color=mm_color)) + 
    geom_point(aes(size=-log10(Pvalue),shape=group)) +
    labs(x="R2",y="- Log10Pvalue",fill="") +
    scale_color_manual(values =unique(cor_pp$mm_color) )+
    theme(panel.grid.major = element_line(colour="grey",linetype = "dotted"),
          panel.background = element_blank(),
          panel.border=element_rect(fill='transparent', color='black'), 
          axis.ticks.x = element_blank(),
          axis.title.x  = element_text(size = 16),
          axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
          legend.background = element_blank())+
    geom_hline(aes(yintercept = -log10(0.05)),colour="grey",linetype=2,size=1)
  library(pheatmap)
  pp<-as.data.frame(moduleTraitCor1)
  pp<-pp[order(pp$asym),]
  anno_row <- data.frame("Module"=rownames(moduleTraitCor1))
  rownames(anno_row)<-rownames(pp)
  ann_colors <- list(Module=c(MEbrown="#7C3F07",MEmagenta="#FF80FFFF",MEyellow="#FFFF00FF",
                              MEblue="#0080FFFF",MEblack="black",MEpurple="purple",MEcyan="cyan",
                              MEmidnightblue="midnightblue",MEgrey60="grey60",MElightgreen="lightgreen",
                              MEred="firebrick",MEgreen="green" ,MElightcyan="#e2fcfb",
                              MEgrey="grey",MEpink= "pink", MEturquoise="#1B9E77",
                              MEtan="tan",MEsalmon="#FF8B8B",MEgreenyellow="#A0D600"))
  
  ann_colors <- list(Module=c(MEbrown="brown",MEmagenta="magenta",MEyellow="yellow",
                              MEblue="blue",MEblack="black",MEpurple="purple",MEcyan="cyan",MEgreen="green",
                              MEgreenyellow='greenyellow',MElightyellow="lightyellow",MEroyalblue="royalblue",
                              MEmidnightblue="midnightblue",MEgrey60="grey60",MElightgreen="lightgreen",
                              MEred="red",MEgreen="green" ,MElightcyan="lightcyan",
                              MEgrey="grey",MEpink= "pink", MEturquoise="turquoise",
                              MEtan="tan",MEsalmon="salmon"))
  pheatmap(pp,cluster_rows=F,cluster_cols=F,
           display_numbers=T,number_format="%.2f",
           border="white",
           fontsize_number=8,
           fontsize_col = 8,
           fontsize_row = 8,
           #angle_col = 90,
           legend_breaks=c(-0.5,-0.2,0,0.2,0.5),
           legend_labels=c("-0.5","-0.2","0","0.2","0.5"),
           #color = mycol,
           fontsize=8,
          annotation_row=anno_row,
          annotation_colors=ann_colors
           )
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(9,8)
  par(cex = 1)
  plotEigengeneNetworks(MET, '', marDendro = c(0,4,1,6), marHeatmap = c(5,6,2,1), cex.lab = 1.2, xLabelsAngle  = 90)

  cell_type <- read.table("/home/dulab/Documents/wrok/flu_paper/data/CODE_FI/ciber/TME.results.output_all.txt",header = T
                          ,row.names = 1,sep ="\t" )
 # cell_type <- read.table("D:/wrok/flu_paper/data/result2/cibersort/TME.results.output_all_nosex.txt",header = T
  #                        ,row.names = 1,sep ="\t" )
 # cell_type <- cell_type[,-c(17,18)]
  cell_type<-cell_type[,c(1,10,18)]
  cell_type <- cbind(design,cell_type)
  library(dplyr)
  group_by(cell_type,asym)%>%summarize_each(funs(mean))
  ##The first principal component of each module
  MEs0 = moduleEigengenes(exprs, moduleColors)$eigengenes
  MEs=orderMEs(MEs0)
  moduleTraitCor = cor(MEs, cell_type,use = "p",method = "spearman")
  moduleTraitCor1<-moduleTraitCor[order(moduleTraitCor[,"sym"],decreasing = TRUE),]#B.cells.naive
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor1, nSamples)
  textMatrix = paste(signif(moduleTraitCor1, 2), "\n(",signif(moduleTraitPvalue, 2), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor1)
  sizeGrWindow(10,18)
  par(mar = c(10, 12, 3,2))#mar = c(4, 10, 3, 2)
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor1,
                 xLabels = colnames(cell_type),
                 yLabels = rownames(moduleTraitCor1),
                 ySymbols = rownames(moduleTraitCor1),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(100),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 #main = paste("Module-trait relationships")
  )
  
  MET = orderMEs(cbind(MEs, cell_type))
  ###################
  pp<-as.data.frame(moduleTraitCor1)
  pp<-pp[order(pp$asym),]
  anno_row <- data.frame("Module"=rownames(moduleTraitCor1))
  rownames(anno_row)<-rownames(pp)
 
  pheatmap(pp,cluster_rows=F,cluster_cols=F,
           display_numbers=T,number_format="%.2f",
           border="white",
           fontsize_number=8,
           fontsize_col = 8,
           fontsize_row = 8,
           #angle_col = 90,
           legend_breaks=c(-0.7,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.7),
           legend_labels=c("-0.7","-0.6","-0.4","-0.2","0","0.2","0.4","0.6","0.7"),
           #color = mycol,
           fontsize=8,
           annotation_row=anno_row,
           annotation_colors=ann_colors
  )
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(9,8)
  par(cex = 1)
  plotEigengeneNetworks(MET, '', marDendro = c(0,4,1,6), marHeatmap = c(5,6,2,1), cex.lab = 1.2, xLabelsAngle  = 90)


#===============================[ Extract all module genes and save matrix to cytoscape ]===========================

  modNames = substring(names(MEs), 3)
  probes = colnames(exprs) 
  for (module in modNames){
    inModule = (moduleColors==module)
    modProbes = probes[inModule]    
    ###output the module genes
    write.csv(modProbes,paste(module,'genes_all.csv',sep = ''))
    
    
    ##output tom matrix :genes  to vis and cytoscpae
    TOM = TOMsimilarityFromExpr(exprs, power = 9)
    modTOM = TOM[inModule, inModule]
    dimnames(modTOM) = list(modProbes, modProbes)
    
    vis = exportNetworkToVisANT(modTOM,
                                file = paste("VisANTInput-median-", module, ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0)
    
    
    
    cyt = exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("CytoscapeInput-edges-median-", paste(module, collapse="-"), ".txt", sep=""),
      nodeFile = paste("CytoscapeInput-nodes-median", paste(module, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.02,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule]
    )
  }
