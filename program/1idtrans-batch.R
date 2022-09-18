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
  
  setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  
  
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
  #setwd('~/result1/0rawdata/')
  gse30550_sym<-read.csv('gse30550_sym_matrix.csv',row.names = 1)
  gse30550_asym<-read.csv('gse30550_asym_matrix.csv',row.names = 1)
  
  gse61754_sym_shedding<-read.csv('gse61754_sym_shedding.csv',row.names = 1)
  gse61754_asym_nonshedding<-read.csv('gse61754_asym_nonshedding.csv',row.names = 1)
  
  gse17156_sym<-read.csv('gse17156_sym.csv',row.names = 1)
  gse17156_asym<-read.csv('gse17156_asym.csv',row.names = 1)
}

#=================================[ pre-processing data ]==========================
{
  ##combine the two groups
  gse30550_exprs0<-as.data.frame(t(rbind(gse30550_sym,gse30550_asym))) 
  #gse30550_group0<-c(rep(1,9),rep(0,8))
  gse30550_group <- read.csv("gse30550_group.csv",row.names = 1)
  
  gse61754_sym_shedding_exprs0<-as.data.frame(cbind(gse61754_sym_shedding,gse61754_asym_nonshedding))
 # gse61754_sym_shedding_group0<-c(rep(1,5),rep(0,3))
  gse61754_group <- read.csv("gse61754_group_asym_nonshedding.csv",row.names = 1)
  
  gse17156_exprs0<-as.data.frame(cbind(gse17156_sym,gse17156_asym))
  #gse17156_group0<-c(rep(1,8),rep(0,9))
  gse17156_group <- read.csv("gse17156_group.csv",row.names = 1)
 
}


#===============================[ transform ID ]=================================

    
    probe1symbol<-function(pak,exprs){
      ids<-toTable(pak)
      ids$probe_id<-paste0('X',ids$gene_id,"_at")
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
    
  
  #gpl9188 <- read.table("D:/wrok/flu_paper/data/GPL9188_1.txt",header = T,sep="\t")
  ##apply the function to all of data
  {
    gse30550_exprs<-probe1symbol(org.Hs.egSYMBOL,gse30550_exprs0)   
    gse61754_sym_shedding_exprs<-probe2symbol(illuminaHumanv4SYMBOL,gse61754_sym_shedding_exprs0)
    gse17156_exprs<-probe2symbol(hgu133a2SYMBOL,gse17156_exprs0)
  }
  
  ##export the trans-ID expression matrix
  {
    write.csv(gse30550_exprs,'gse30550_exprs.csv')
    write.csv(gse61754_sym_shedding_exprs,'gse61754_sym_shedding_exprs')
    write.csv(gse17156_exprs,'gse17156_exprs')
  }
  

  


#[================================================ combine all the expression matrix ===================================] 
{
  
  gse61754_sym_shedding_exprs<- read.csv("gse61754_sym_shedding_exprs",row.names = 1)
  gse17156_exprs<- read.csv("gse17156_exprs",row.names = 1)
  gse30550_exprs<- read.csv("gse30550_exprs.csv",row.names = 1)
  
  
  #####
  exprs00_h3n2<-Reduce(function(x,y) merge(x,y,by='symbol'),list(gse17156_exprs,gse30550_exprs,
                                                                 gse61754_sym_shedding_exprs),accumulate = FALSE)
  row.names(exprs00_h3n2)<-exprs00_h3n2$symbol
  exprs00_h3n2<-exprs00_h3n2[,-1]
  batch_h3n2<-data.frame(subject=colnames(exprs00_h3n2),batch=factor(c(rep("GSE17156",17),rep("GSE30550",17),rep("GSE61754",8))))
  group_h3n2 <- rbind(gse17156_group,gse30550_group,gse61754_group)
  group_h3n2 <- merge(batch_h3n2,group_h3n2,by="subject",all.x=T)
  
  ##############
  
  range_fun <- function(x){
    if(nrow(x)>=2){
      gseID<- paste(x$gseID[1],x$gseID[2],sep = ",")
    }else{
      gseID<-x$gseID
    }
   y<-data.frame(gseID=gseID,subject=x$subject[1],group=x$group[1],batch=x$batch[1]) 
  return(y)}
  
  group_h3n2_list <- split(group_h3n2,list(group_h3n2$subject))
  group_h3n21<-as.data.frame(t(sapply(group_h3n2_list,range_fun)))
  group_h3n21$gseID<-unlist(group_h3n21$gseID)
  group_h3n21$subject<-unlist(group_h3n21$subject)
  group_h3n21$group<-unlist(group_h3n21$group)
  group_h3n21$batch<-unlist(group_h3n21$batch)
  setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  
  
  
  write.csv(exprs00_h3n2,'exprs00_h3n2.csv')
  write.csv(group_h3n2,'group_h3n2_whole.csv')
  write.csv(group_h3n21,'group_h3n2.csv')

}
  


#=================[ check the overall expression and do the chip normalization and do PCA to check the quality of the data ]===================

    exprs00_h3n2<-read.csv('exprs00_h3n2.csv',row.names = 1)  #87*8286
    h3n2_group <- read.csv("group_h3n2.csv",row.names = 1)
    
    exprs00_h3n2_1 <- as.data.frame(exprs00_h3n2)
    exprs00_h3n2_1$gene <- rownames(exprs00_h3n2_1)
    exprs00_h3n2_2 <- melt(exprs00_h3n2_1,id="gene")
    colnames(exprs00_h3n2_2)<- c("gene","subject","expression")
    exprs00_h3n2_2 <- merge(exprs00_h3n2_2,h3n2_group,by="subject")
    
    pdf('H3N2_before1.pdf',width = 15, height = 10 )
    pp<-ggplot(exprs00_h3n2_2,aes(x=subject,y=expression))+geom_boxplot(aes(fill=group)) +#,outlier.colour = NA
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
    colnames(exprs00_h3n2_2)<- c("gene","subject","expression")
    exprs00_h3n2_2 <- merge(exprs00_h3n2_2,h3n2_group,by="subject")
    
    pdf('H3N2_after1.pdf',width = 15, height = 10 )
    pp<-ggplot(exprs00_h3n2_2,aes(x=subject,y=expression))+geom_boxplot(aes(fill=group)) +#,outlier.colour = NA
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
        pca.data<-data.frame(subject=row.names(pca$x),
                             x=pca$x[,1],
                             y=pca$x[,2])
        pca.data<- merge(pca.data,group,by="subject")
        #barplot(pca.var.per,main =paste0('pca var of',exprs,'group') ,xlab = 'Principal Comonent',ylab = 'Percent Variation')
        pp<- ggplot(data=pca.data,aes(x=x,y=y))+#,shape=group
          geom_point(aes(colour=batch),size=5)+mytheme+ scale_color_lancet()+
          xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
          ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
        pdf(paste("/home/dulab/Documents/wrok/flu_paper/data/result_newest/",picturename,"batch.pdf",sep = "_"),height=8,width=8)
        print(pp)
        dev.off()
        png(paste("/home/dulab/Documents/wrok/flu_paper/data/result_newest/",picturename,"batch.png",sep = "_"),width = 800, height = 800)
        print(pp)
        dev.off()
        pp<- ggplot(data=pca.data,aes(x=x,y=y))+geom_point(aes(colour=group),size=5)+ scale_color_lancet()+mytheme+
          xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
          ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
        pdf(paste("/home/dulab/Documents/wrok/flu_paper/data/result_newest/",picturename,"group.pdf",sep = "_"),height=8,width=8)
        print(pp)
        dev.off()
        png(paste("/home/dulab/Documents/wrok/flu_paper/data/result_newest/",picturename,"group.png",sep = "_"),width = 800, height = 800)
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
    write.csv(combat_edata, paste("/home/dulab/Documents/wrok/flu_paper/data/result_newest/",dataname,"_combat_data.csv",sep = ""))
  }
  
  ##apply the remove_batch function

    remove_batch(exprs00_h3n21,h3n2_group,"h3n2")

#[=============== check the performance use PCA again after remove batch effect ==============]

  ##load data after removing batch effects
 
    exprs00_h3n2_combat<-read.csv('h3n2_combat_data.csv',row.names = 1)
  
    do_pca(exprs00_h3n2_combat,h3n2_group,"combat_h3n2")






