#=============================[ extract the standardized data from GEO ]============================
rm(list=ls())
options(stringsAsFactors = F)

#==========================[ load packages ]========================
{
  library(BiocManager)
  # BiocManager::install("sva")
  #BiocManager::install("GEOquery")
  # BiocManager::install("illuminaHumanv4.db")
  #BiocManager::install("hgu133a2.db")
  # BiocManager::install("org.Hs.eg.db")
  library(hgu133a2.db)
  library(org.Hs.eg.db)
  library(illuminaHumanv4.db)
  library(GEOquery)
  library(tidyr)
  library(sva)
  library(dplyr)
  library(limma)
  library(openxlsx)
  #library(xlsx)
  library(RankProd)
  
}

#=========================[ set working directory ]==========================

setwd("/home/dulab/Documents/wrok/flu_paper/data/result2/result_new/")  

#======================[ load RMA normalized data from GEO ]========================

  ##load the data
    #GSE30550
    load("GSE30550_eSet.Rdata")#raw_data/gse30550/data/
    exprs3<-exprs(gset[[1]])
    group3<-pData(gset[[1]])

#=============================[ data pre-processing ]=======================

  #================[ extract the baseline time point,the influenza subject,group by symptoms/(symptoms and shedding) and it's expression matrix ]===============
  ##gse30550

    ###extract the pre-0 time point and the flu subject
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/") 
#######################
###grouping by the symptom and shedding
group3_sym_index<-group3[group3$"clinic_pheno:ch1"=='Symptomatic',]
group3_asym_index<-group3[group3$"clinic_pheno:ch1"=='Asymptomatic',]


subject_name_sym <- apply(data.frame(group3_sym_index$source_name_ch1),1,function(x) unlist(strsplit(paste(unlist(strsplit(x,split = " "))[4],unlist(strsplit(x,split = " "))[5],sep = ""),split = ","))[1])
subject_name_asym<- apply(data.frame(group3_asym_index$source_name_ch1),1,function(x) unlist(strsplit(paste(unlist(strsplit(x,split = " "))[4],unlist(strsplit(x,split = " "))[5],sep = ""),split = ","))[1])

gse30550_group_sym  <- data.frame(gseID=rownames(group3_sym_index),Time=group3_sym_index$`time_hpi:ch1`,subject=subject_name_sym,group='Symptomatic')
gse30550_group_asym  <- data.frame(gseID=rownames(group3_asym_index),Time=group3_asym_index$`time_hpi:ch1`,subject=subject_name_asym,group='Asymptomatic')
gse30550_group <- rbind( gse30550_group_sym, gse30550_group_asym)
rownames( gse30550_group)<- gse30550_group$gseID
write.csv(gse30550_group,"gse30550_group_alltime.csv")

###generate the mean matrix
gse30550_sym<-as.data.frame(t(exprs3[,group3_sym_index$geo_accession]))
gse30550_asym<-as.data.frame(t(exprs3[,group3_asym_index$geo_accession]))


write.csv(gse30550_sym,'gse30550_sym_matrix_alltime.csv',row.names = T)
write.csv(gse30550_asym,'gse30550_asym_matrix_alltime.csv',row.names = T)

#######################################################################        
gse30550_sym_sub<-read.csv('gse30550_sym_matrix_alltime.csv',row.names = 1)
gse30550_asym_sub<-read.csv('gse30550_asym_matrix_alltime.csv',row.names = 1)
sym_shedding_exprs0<-as.data.frame(t(rbind(gse30550_sym_sub,gse30550_asym_sub)))
        # dee2_sym_shedding_group0<-c(rep(1,9),rep(0,6))
        
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
        sym_shedding_exprs<-probe1symbol(org.Hs.egSYMBOL,sym_shedding_exprs0)
        
        write.csv(sym_shedding_exprs[,-ncol(sym_shedding_exprs)],"gse30550_sym_shedding_exprs_alltime.csv",row.names = T)
##################################deg
        exprs_h3n2<-read.csv('gse30550_sym_shedding_exprs_alltime.csv',row.names = 1)   #72*8286
        #exprs<-exprs_h3n2[1:50,]
        #exprs.gnames<-rownames(exprs_h3n2)
        gse30550_group <- read.csv("gse30550_group_alltime.csv",row.names = 1)
        gse30550_group$batch <- "gse30550"
        #h3n2_group<-h3n2_group[colnames(exprs_h3n2),]
        re_arrange_fun <- function(x,y,basetime){
          FC <-data.frame()
          degsss <- data.frame()
         tt <- setdiff(x$Time,basetime)
         for (gg in unique(x$group)){
           sub_group1 <- subset(x,group==gg)
           for(i in tt){
             sub_group <- subset(sub_group1,Time==basetime|Time==i)
             exprs.cl <- as.numeric(as.factor(sub_group$Time))
             exprs.cl[exprs.cl==1]<-1
             exprs.cl[exprs.cl==2]<-0
             sub_y <- y[,sub_group$gseID]
             RP.adv.out <- RP.advance(sub_y,exprs.cl,sub_group$batch, logged=TRUE,na.rm=FALSE,  gene.names=rownames(sub_y),
                                      plot = FALSE,rand = 123,MinNumOfValidPairs =1)     ###class2:time0 ,class1:other time
             
             degs<-topGene(RP.adv.out,cutoff = 0.05,method="pfp",logged=TRUE,logbase=2,gene.names=rownames(sub_y)) 
             if(is.null(degs$Table1)==0){
              down<-as.data.frame(degs$Table1)
             down$degtype <- "down" 
             write.csv(down,paste(i,gg,'downdeg_h3n2.csv',sep = "_"))
             if(is.null(degs$Table2)==0){
            up<-as.data.frame(degs$Table2)
             up$degtype <- "up"
             write.csv(up,paste(i,gg,'updeg_h3n2.csv',sep="_"))  
            degss <- rbind(up,down) 
             degss$type <- i
             degss$group <-gg
             degss$gene<-rownames(degss)
             degsss <- rbind(degsss,degss)
             write.csv(degss,paste(i,gg,'deg_h3n2.csv',sep = "_")) 
             }else {
               degss <- down 
               degss$type <- i
               degss$group <-gg
               degss$gene<-rownames(degss)
               degsss <- rbind(degsss,degss)
               write.csv(degss,paste(i,gg,'deg_h3n2.csv',sep = "_")) 
            }
                
             }else{
               if(is.null(degs$Table2)==0){
                 up<-as.data.frame(degs$Table2)
                 up$degtype <- "up"
                 write.csv(up,paste(i,gg,'updeg_h3n2.csv',sep="_"))  
                 degss <- up 
                 degss$type <- i
                 degss$group <-gg
                 degss$gene<-rownames(degss)
                 degsss <- rbind(degsss,degss)
                 write.csv(degss,paste(i,gg,'deg_h3n2.csv',sep = "_")) 
               }else {
                 next;
               }
             }
             geneFC <- RP.adv.out$AveFC
             colnames(geneFC)<-c("AveFC")
             geneFC <- as.data.frame(geneFC)
             geneFC$pfp <-RP.adv.out$pfp
             geneFC$pval <- RP.adv.out$pval
             geneFC$type <- i
             geneFC$group <-gg
             geneFC$gene<-rownames(geneFC)
             write.csv(geneFC,paste(i,gg,"gene_fc.csv",sep="_"),row.names = F,quote =F)
             FC <- rbind(FC,geneFC)
        print(i) }
           print(gg)
           }
         write.csv(FC,"fc_all.csv",row.names = F,quote =F)
         write.csv(degsss,"deg_all.csv",row.names = F,quote =F)
         
        return(FC)}
       
       FC<- re_arrange_fun(gse30550_group,exprs_h3n2,basetime = "Hour 00")
       
       
       re_arrange_fun <- function(x,y,basetime){
         FC <-data.frame()
         degsss <- data.frame()
         tt <- setdiff(x$Time,basetime)
         for (gg in unique(x$group)){
           sub_group1 <- subset(x,group==gg)
           for(i in tt){
             sub_group <- subset(sub_group1,Time==basetime|Time==i)
             exprs.cl <- as.numeric(as.factor(sub_group$Time))
             exprs.cl[exprs.cl==1]<-1
             exprs.cl[exprs.cl==2]<-0
             sub_y <- y[,sub_group$gseID]
             RP.adv.out <- RP.advance(sub_y,exprs.cl,sub_group$batch, logged=TRUE,na.rm=FALSE,  gene.names=rownames(sub_y),
                                      plot = FALSE,rand = 123,MinNumOfValidPairs =1)     ###class2:time0 ,class1:other time
             
             degs<-topGene(RP.adv.out,cutoff = 0.05,method="pfp",logged=TRUE,logbase=2,gene.names=rownames(sub_y)) 
             if(is.null(degs$Table1)==0){
               down<-as.data.frame(degs$Table1)
               down$degtype <- "down" 
               write.csv(down,paste(i,gg,'downdeg_h3n21.csv',sep = "_"))
               if(is.null(degs$Table2)==0){
                 up<-as.data.frame(degs$Table2)
                 up$degtype <- "up"
                 write.csv(up,paste(i,gg,'updeg_h3n21.csv',sep="_"))  
                 degss <- rbind(up,down) 
                 degss$type <- i
                 degss$group <-gg
                 degss$gene<-rownames(degss)
                 degsss <- rbind(degsss,degss)
                 write.csv(degss,paste(i,gg,'deg_h3n21.csv',sep = "_")) 
               }else {
                 degss <- down 
                 degss$type <- i
                 degss$group <-gg
                 degss$gene<-rownames(degss)
                 degsss <- rbind(degsss,degss)
                 write.csv(degss,paste(i,gg,'deg_h3n21.csv',sep = "_")) 
               }
               
             }else{
               if(is.null(degs$Table2)==0){
                 up<-as.data.frame(degs$Table2)
                 up$degtype <- "up"
                 write.csv(up,paste(i,gg,'updeg_h3n21.csv',sep="_"))  
                 degss <- up 
                 degss$type <- i
                 degss$group <-gg
                 degss$gene<-rownames(degss)
                 degsss <- rbind(degsss,degss)
                 write.csv(degss,paste(i,gg,'deg_h3n21.csv',sep = "_")) 
               }else {
                 next;
               }
             }
             geneFC <- RP.adv.out$AveFC
             colnames(geneFC)<-c("AveFC")
             geneFC <- as.data.frame(geneFC)
             geneFC$pfp <-RP.adv.out$pfp
             geneFC$pval <- RP.adv.out$pval
             geneFC$type <- i
             geneFC$group <-gg
             geneFC$gene<-rownames(geneFC)
             write.csv(geneFC,paste(i,gg,"gene_fc1.csv",sep="_"),row.names = F,quote =F)
             FC <- rbind(FC,geneFC)
             print(i) }
           print(gg)
         }
         write.csv(FC,"fc_all1.csv",row.names = F,quote =F)
         write.csv(degsss,"deg_all1.csv",row.names = F,quote =F)
         
         return(FC)}
       
       FC<- re_arrange_fun(gse30550_group,exprs_h3n2,basetime = "Baseline")