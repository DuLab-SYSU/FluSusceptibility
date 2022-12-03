#==============[ Probe ID transform to symbol-->chip normalization --> PCA --> remove batch effects --> PCA ]==========
rm(list = ls())

#===========================[ load packages ]========================
{
  library(hgu133a2.db)
  library(hgu95av2.db)
  library(org.Hs.eg.db)
  library(illuminaHumanv4.db)
  library(GEOquery)
  library(tidyr)
  library(sva)
  library(bladderbatch)
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(reshape)
  library(reshape2)
  library(RColorBrewer)
  library(ggsci)
  
  setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  
  
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
  }

#================================[ load data ]=================================
{
  
  gse61754_sym<-read.csv('gse61754_sym.csv',row.names = 1)
  gse61754_asym<-read.csv('gse61754_asym.csv',row.names = 1)
  
  
  dee2_sym<-read.csv('dee2_sym_mean_matrix.csv',row.names = 1)
  dee2_asym<-read.csv('dee2_asym_mean_matrix.csv',row.names = 1)
  
 
  dee5_sym<-read.csv('dee5_sym_mean_matrix.csv',row.names = 1)
  dee5_asym<-read.csv('dee5_asym_mean_matrix.csv',row.names = 1)
  
 
}

#=================================[ pre-processing data ]==========================
{

  
  gse61754_sym_exprs0<-as.data.frame(cbind(gse61754_sym,gse61754_asym))
 # gse61754_sym_shedding_group0<-c(rep(1,5),rep(0,3))
  gse61754_group <- read.csv("gse61754_group.csv",row.names = 1)
  
 
  
  dee2_sym_exprs0<-as.data.frame(cbind(dee2_sym,dee2_asym))
  gse73072_group <- read.csv("gse73072_group.csv",row.names = 1)
  
  
  dee5_sym_exprs0<-as.data.frame(cbind(dee5_sym,dee5_asym))
  dee5_group <- read.csv("dee5_group.csv",row.names = 1)
 
}


#===============================[ transform ID ]=================================

    
    probe1symbol<-function(pak,exprs){
      ids<-toTable(pak)
      ids$probe_id<-paste0(ids$gene_id,"_at")
      exprs<-exprs[row.names(exprs) %in% ids$probe_id,]
      ids<-ids[ids$probe_id %in% row.names(exprs),]
      exprs$mean <-apply(exprs, 1, mean)
      ids$mean<-exprs[ids$probe_id,'mean']
      ids<-ids[order(ids$symbol,ids$mean,decreasing = T),]
      ids<-ids[!duplicated(ids$symbol),]
      exprs1<-exprs[ids$probe_id,]
      row.names(exprs1)<-ids$symbol
      exprs1<-na.omit(exprs1)
      exprs1$symbol<-row.names(exprs1)
      exprs1<-exprs1[,-ncol(exprs)]
      return(exprs1)
    }
    
    
    probe2symbol<-function(pak,exprs){
      ids<-toTable(pak)
      #ids$probe_id<-paste0('X',ids$probe_id)
      exprs<-exprs[row.names(exprs) %in% ids$probe_id,]
      ids<-ids[ids$probe_id %in% row.names(exprs),]
      exprs$mean <-apply(exprs, 1, mean)
      ids$mean<-exprs[ids$probe_id,'mean']
      ids<-ids[order(ids$symbol,ids$mean,decreasing = T),]
      ids<-ids[!duplicated(ids$symbol),]
      exprs1<-exprs[ids$probe_id,]
      row.names(exprs1)<-ids$symbol
      exprs1<-na.omit(exprs1)
      exprs1$symbol<-row.names(exprs1)
      exprs1<-exprs1[,-ncol(exprs)]
      return(exprs1)
    }
    
  
  ##apply the function to all of data
  {
    gse61754_sym_exprs<-probe2symbol(illuminaHumanv4SYMBOL,gse61754_sym_exprs0)

    dee2_sym_exprs<-probe1symbol(org.Hs.egSYMBOL,dee2_sym_exprs0)
    
    dee5_sym_exprs<-probe1symbol(org.Hs.egSYMBOL,dee5_sym_exprs0)
  }
  
  ##export the trans-ID expression matrix
  {
    write.csv(gse61754_sym_exprs,'gse61754_sym_exprs')

    write.csv(dee2_sym_exprs,'dee2_sym_exprs')
    
    write.csv(dee5_sym_exprs,'dee5_sym_exprs')
  }
  

  


#[================================================ combine all the expression matrix ===================================] 
{
  dee2_sym_exprs<- read.csv("dee2_sym_exprs",row.names = 1)
  dee5_sym_exprs<- read.csv("dee5_sym_exprs",row.names = 1)
  gse61754_sym_exprs<- read.csv("gse61754_sym_exprs",row.names = 1)

  
  
  #####
  exprs00_h3n2<-Reduce(function(x,y) merge(x,y,by='symbol'),list(dee2_sym_exprs,dee5_sym_exprs,
                                                                 gse61754_sym_exprs),accumulate = FALSE)
  row.names(exprs00_h3n2)<-exprs00_h3n2$symbol
  exprs00_h3n2<-exprs00_h3n2[,-1]
  batch_h3n2<-data.frame(gseID=colnames(exprs00_h3n2),batch=factor(c(rep("GSE7072_DEE2",17),rep("GSE7072_DEE5",21),rep("GSE61754",11))))
  group_h3n2 <- rbind(gse73072_group,dee5_group,gse61754_group)
  group_h3n2 <- merge(batch_h3n2,group_h3n2,by="gseID",all.x=T)
  rownames(group_h3n2)<-group_h3n2$gseID
  setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  
  
  
  write.csv(exprs00_h3n2,'exprs00_h3n2.csv')
  write.csv(group_h3n2,'group_h3n2.csv')

  
  dee2_sym_exprs<- read.csv("dee2_sym_exprs",row.names = 1)
  dee5_sym_exprs<- read.csv("dee5_sym_exprs",row.names = 1)
  gse61754_sym_exprs<- read.csv("gse61754_sym_exprs",row.names = 1)

  
  
  #####
  exprs00_h3n2<-Reduce(function(x,y) merge(x,y,by='symbol'),list(dee2_sym_exprs,dee5_sym_exprs
                                                                 ),accumulate = FALSE)
  row.names(exprs00_h3n2)<-exprs00_h3n2$symbol
  exprs00_h3n2<-exprs00_h3n2[,-1]
  batch_h3n2<-data.frame(gseID=colnames(exprs00_h3n2),batch=factor(c(rep("GSE7072_DEE2",17),rep("GSE7072_DEE5",21))))
  group_h3n2 <- rbind(gse73072_group,dee5_group)
  group_h3n2 <- merge(batch_h3n2,group_h3n2,by="gseID",all.x=T)
  rownames(group_h3n2)<-group_h3n2$gseID
  setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  
  
  
  write.csv(exprs00_h3n2,'exprs00_h3n2_train.csv')
  write.csv(group_h3n2,'group_h3n2_train.csv')
  
}
  


#=================[ check the overall expression and do the chip normalization and do PCA to check the quality of the data ]===================

    exprs00_h3n2<-read.csv('exprs00_h3n2.csv',row.names = 1)  #87*8286
    h3n2_group <- read.csv("group_h3n2.csv",row.names = 1)
   
    exprs00_h3n2_1 <- as.data.frame(exprs00_h3n2)
    exprs00_h3n2_1$gene <- rownames(exprs00_h3n2_1)
    exprs00_h3n2_2 <- melt(exprs00_h3n2_1,id="gene")
    colnames(exprs00_h3n2_2)<- c("gene","gseID","expression")
    exprs00_h3n2_2 <- merge(exprs00_h3n2_2,h3n2_group,by="gseID")
    
    pdf('H3N2_before1.pdf',width = 15, height = 6 )
    pp<-ggplot(exprs00_h3n2_2,aes(x=gseID,y=expression))+geom_boxplot(aes(fill=group)) +#,outlier.colour = NA
      scale_fill_manual(values = c("#085A9C","#EF0808"))+mytheme+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      guides(fill=guide_legend(title=NULL))
    print(pp)
     dev.off()
     png('H3N2_before1.png',width = 1500, height = 800)
     print(pp)
     dev.off()
     
    exprs00_h3n21<-normalizeBetweenArrays(exprs00_h3n2)   ##use the limma to do the chip normalization
   # boxplot(exprs00_h3n21,outline=FALSE,notch=T,col=group_h3n2,las=2,main='the distribution of H3N2_group after correction')
    exprs00_h3n2_1 <- as.data.frame(exprs00_h3n21)
    exprs00_h3n2_1$gene <- rownames(exprs00_h3n2_1)
    exprs00_h3n2_2 <- melt(exprs00_h3n2_1,id="gene")
    colnames(exprs00_h3n2_2)<- c("gene","gseID","expression")
    exprs00_h3n2_2 <- merge(exprs00_h3n2_2,h3n2_group,by="gseID")
    
    pdf('H3N2_after1.pdf',width = 15, height = 6 )
    pp<-ggplot(exprs00_h3n2_2,aes(x=gseID,y=expression))+geom_boxplot(aes(fill=group)) +#,outlier.colour = NA
      scale_fill_manual(values = c("#085A9C","#EF0808"))+mytheme+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      guides(fill=guide_legend(title=NULL))
    print(pp)
    dev.off()
    png('H3N2_after1.png',width = 1500, height = 800)
    print(pp)
    dev.off()
    
    
    
  
  ##do PCA
  
      do_pca<-function(exprs,group,picturename){
        pca<-prcomp(t(exprs),scale= T)
        pca.var<-pca$sdev
        pca.var.per<-round(pca.var/sum(pca.var)*100 ,1)
        pca.data<-data.frame(gseID=row.names(pca$x),
                             x=pca$x[,1],
                             y=pca$x[,2])
        pca.data<- merge(pca.data,group,by="gseID")
        #barplot(pca.var.per,main =paste0('pca var of',exprs,'group') ,xlab = 'Principal Comonent',ylab = 'Percent Variation')
        pp<- ggplot(data=pca.data,aes(x=x,y=y))+#,shape=group
          geom_point(aes(colour=batch,shape=group),size=5)+mytheme+ scale_color_lancet()+
          xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
          ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
        pdf(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"batch.pdf",sep = "_"),height=8,width=8)
        print(pp)
        dev.off()
        png(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"batch.png",sep = "_"),width = 800, height = 800)
        print(pp)
        dev.off()
        pp<- ggplot(data=pca.data,aes(x=x,y=y))+geom_point(aes(colour=group),size=5)+ scale_color_lancet()+mytheme+
          xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
          ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
        pdf(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"group.pdf",sep = "_"),height=8,width=8)
        print(pp)
        dev.off()
        png(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"group.png",sep = "_"),width = 800, height = 800)
        print(pp)
        dev.off()
      }
   
      ##use the function
      do_pca(exprs00_h3n2,h3n2_group,"scale_before")
      do_pca(exprs00_h3n21,h3n2_group,"scale_after")




  
#==========================[ processing the batch effects ]============================

  remove_batch<-function(exprs,group_infor,dataname){
    exprs<-as.matrix(exprs)
    group_infor<-group_infor[colnames(exprs),]
    mod0<-model.matrix(~as.factor(group),data = group_infor)
    combat_edata = ComBat(dat=exprs, batch=group_infor$batch, mod=mod0, par.prior=TRUE, prior.plots=FALSE)
    write.csv(combat_edata, paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",dataname,"_combat_data.csv",sep = ""))
  }
  
  ##apply the remove_batch function

    remove_batch(exprs00_h3n21,h3n2_group,"h3n2")
    #remove_batch(exprs00_h1n11,batch_h1n1,group_h1n1,"h1n1")
  



#[=============== check the performance use PCA again after remove batch effect ==============]

  ##load data after removing batch effects
 
    exprs00_h3n2_combat<-read.csv('h3n2_combat_data.csv',row.names = 1)

  ##apply the do_pca function to check the performance of removing batch effects
  
    do_pca(exprs00_h3n2_combat,h3n2_group,"combat_h3n2")

  
 #####################3
    exprs00_h3n2<-read.csv('exprs00_h3n2_train.csv',row.names = 1)  #87*8286
    h3n2_group <- read.csv("group_h3n2_train.csv",row.names = 1)
   
    
    #batch_group_h3n2 <- data.frame(sample=colnames(exprs00_h3n2),group_h3n2=group_h3n2,batch=batch_h3n2)
    exprs00_h3n2_1 <- as.data.frame(exprs00_h3n2)
    exprs00_h3n2_1$gene <- rownames(exprs00_h3n2_1)
    exprs00_h3n2_2 <- melt(exprs00_h3n2_1,id="gene")
    colnames(exprs00_h3n2_2)<- c("gene","gseID","expression")
    exprs00_h3n2_2 <- merge(exprs00_h3n2_2,h3n2_group,by="gseID")
    
    pdf('H3N2_before1_train.pdf',width = 15, height = 6 )
    pp<-ggplot(exprs00_h3n2_2,aes(x=gseID,y=expression))+geom_boxplot(aes(fill=group)) +#,outlier.colour = NA
      scale_fill_manual(values = c("#085A9C","#EF0808"))+mytheme+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      guides(fill=guide_legend(title=NULL))
    print(pp)
    dev.off()
    png('H3N2_before1_train.png',width = 1500, height = 800)
    print(pp)
    dev.off()
    
    exprs00_h3n21<-normalizeBetweenArrays(exprs00_h3n2)   ##use the limma to do the chip normalization
    # boxplot(exprs00_h3n21,outline=FALSE,notch=T,col=group_h3n2,las=2,main='the distribution of H3N2_group after correction')
    exprs00_h3n2_1 <- as.data.frame(exprs00_h3n21)
    exprs00_h3n2_1$gene <- rownames(exprs00_h3n2_1)
    exprs00_h3n2_2 <- melt(exprs00_h3n2_1,id="gene")
    colnames(exprs00_h3n2_2)<- c("gene","gseID","expression")
    exprs00_h3n2_2 <- merge(exprs00_h3n2_2,h3n2_group,by="gseID")
    
    pdf('H3N2_after1_train.pdf',width = 15, height = 6 )
    pp<-ggplot(exprs00_h3n2_2,aes(x=gseID,y=expression))+geom_boxplot(aes(fill=group)) +#,outlier.colour = NA
      scale_fill_manual(values = c("#085A9C","#EF0808"))+mytheme+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      guides(fill=guide_legend(title=NULL))
    print(pp)
    dev.off()
    png('H3N2_after1_train.png',width = 1500, height = 800)
    print(pp)
    dev.off()
    
    
    
    
    ##do PCA
    
    do_pca<-function(exprs,group,picturename){
      pca<-prcomp(t(exprs),scale= T)
      pca.var<-pca$sdev
      pca.var.per<-round(pca.var/sum(pca.var)*100 ,1)
      pca.data<-data.frame(gseID=row.names(pca$x),
                           x=pca$x[,1],
                           y=pca$x[,2])
      pca.data<- merge(pca.data,group,by="gseID")
      #barplot(pca.var.per,main =paste0('pca var of',exprs,'group') ,xlab = 'Principal Comonent',ylab = 'Percent Variation')
      pp<- ggplot(data=pca.data,aes(x=x,y=y))+#,shape=group
        geom_point(aes(colour=batch,shape=group),size=5)+mytheme+ scale_color_lancet()+
        xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
        ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
      pdf(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"batch.pdf",sep = "_"),height=8,width=8)
      print(pp)
      dev.off()
      png(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"batch.png",sep = "_"),width = 800, height = 800)
      print(pp)
      dev.off()
      pp<- ggplot(data=pca.data,aes(x=x,y=y))+geom_point(aes(colour=group),size=5)+ scale_color_lancet()+mytheme+
        xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
        ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
      pdf(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"group.pdf",sep = "_"),height=8,width=8)
      print(pp)
      dev.off()
      png(paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",picturename,"group.png",sep = "_"),width = 800, height = 800)
      print(pp)
      dev.off()
    }
    
    ##use the function
    do_pca(exprs00_h3n2,h3n2_group,"scale_before_train")
    do_pca(exprs00_h3n21,h3n2_group,"scale_after_train")
    
    
    
    
    
    #==========================[ processing the batch effects ]============================
    
    remove_batch<-function(exprs,group_infor,dataname){
      exprs<-as.matrix(exprs)
      group_infor<-group_infor[colnames(exprs),]
      mod0<-model.matrix(~as.factor(group),data = group_infor)
      combat_edata = ComBat(dat=exprs, batch=group_infor$batch, mod=mod0, par.prior=TRUE, prior.plots=FALSE)
      write.csv(combat_edata, paste("/home/dulab/Documents/wrok/flu_paper/data/result_fi/",dataname,"_combat_data_train.csv",sep = ""))
    }
    
    ##apply the remove_batch function
    
    remove_batch(exprs00_h3n21,h3n2_group,"h3n2")

    
    
    
    #[=============== check the performance use PCA again after remove batch effect ==============]
    
    ##load data after removing batch effects
    
    exprs00_h3n2_combat<-read.csv('h3n2_combat_data_train.csv',row.names = 1)
    
    
    ##apply the do_pca function to check the performance of removing batch effects
    
    do_pca(exprs00_h3n2_combat,h3n2_group,"combat_h3n2_train")






