#=============================[ extract the standardized data from GEO ]============================
rm(list=ls())
options(stringsAsFactors = F)
gc()
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
}

#=========================[ set working directory ]==========================

  #setwd("/home/dulab/Documents/wrok/flu_paper/data/")  
  
#======================[ load RMA normalized data from GEO ]========================

  ###Other data sets were downloaded in the same way

   # f='GSE73072_eSet.Rdata'
  #    gset <- getGEO('GSE73072', destdir=".",
  #                   AnnotGPL = F,    
  #                   getGPL = F)      
  #    save(gset,file=f)  

  
  ##load the data
  {
    #setwd("~/result1/")
    #GSE73072
    load("GSE73072_eSet.Rdata")#raw_data/gse73072/base_data/
    exprs1<-exprs(gset[[1]])
    group1<-pData(gset[[1]])
    #GSE61754
    load("GSE61754_eSet.Rdata")#raw_data/gse61754/data/
    exprs4<-exprs(gset[[1]])
    group4<-pData(gset[[1]])

  }


#=============================[ data pre-processing ]=======================

  #================[ extract the baseline time point,the influenza subject,group by symptoms/(symptoms and shedding) and it's expression matrix ]===============
  ##GSE73072
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi/")  

    ###extract the pre-0 time point and the flu subject
   
      hour<-c('hour 0')
      virus<-c('H3N2')
      gse73072_pre0_gsm<-group1[group1$`time point:ch1` %in% hour,] 
      gse73072_pre0_flu_gsm<-gse73072_pre0_gsm[gse73072_pre0_gsm$`virus:ch1` %in% virus,]
      gse73072_exp<-exprs1[,row.names(gse73072_pre0_flu_gsm)]   ##12023*139
   
    ###grouping by symptoms and shedding and calculate the mean of gene expression value up to baseline time point, generate the mean expression matrix.
    
      ###H3N2 DEE2 have two time point:-23,0
      
        ####grouping index are from paper
        sym_index<-c(1,5,6,7,8,10,12,13,15)
        asym_index<-c(2,3,4,9,11,14,16,17)
        #shedding_index<-c(1,2,4,5,6,7,8,10,12,13,15)
        #nonshedding_index<-c(3,9,11,14,16,17)
        #sym_shedding_index<-c(1,5,6,7,8,10,12,13,15)
        #asym_nonshedding_index<-c(3,9,11,14,16,17)
      
        sym_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE2 subject ',sym_index)),]
        asym_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE2 subject ',asym_index)), ]
        #shedding_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE2 subject ',shedding_index)),]
        #nonshedding_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE2 subject ',nonshedding_index)),]
        #sym_shedding_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in%c(paste0('H3N2 DEE2 subject ',sym_shedding_index)),]
        #asym_nonshedding_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE2 subject ',asym_nonshedding_index)),]
        
        
        #subject_name_sym <- apply(data.frame(sym_shedding_gsm$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE2",sep = "_"))
        #subject_name_asym<- apply(data.frame(asym_nonshedding_gsm$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE2",sep = "_"))
        
        subject_name_sym <- apply(data.frame(sym_gsm$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE2",sep = "_"))
        subject_name_asym<- apply(data.frame(asym_gsm$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE2",sep = "_"))
        
        
        #gse73072_group_sym <- data.frame(gseID=sym_shedding_gsm$geo_accession, subject=subject_name_sym,group="Symptomatic")
        #gse73072_group_asym <- data.frame(gseID=asym_nonshedding_gsm$geo_accession, subject=subject_name_asym,group="Asymptomatic")
        
        gse73072_group_sym <- data.frame(gseID=sym_gsm$geo_accession, subject=subject_name_sym,group="Symptomatic")
        gse73072_group_asym <- data.frame(gseID=asym_gsm$geo_accession, subject=subject_name_asym,group="Asymptomatic")
        
        gse73072_group<- rbind(gse73072_group_sym,gse73072_group_asym)
        rownames(gse73072_group)<-gse73072_group$gseID
        write.csv(gse73072_group,"gse73072_group.csv")
        
        
       
        dee2_sym_mean_matrix0<-as.data.frame(exprs1[,sym_gsm$geo_accession])
        dee2_asym_mean_matrix0<-as.data.frame(exprs1[,asym_gsm$geo_accession])
        
        
        ####Export the mean matrix
        write.csv(dee2_sym_mean_matrix0,'dee2_sym_mean_matrix.csv',quote = F,row.names = T)
        write.csv(dee2_asym_mean_matrix0,'dee2_asym_mean_matrix.csv',row.names = T)
  
     
      
     
      ###H3N2 DEE5 :-30
        hour<-c('hour -30')
        virus<-c('H3N2')
        gse73072_pre0_gsm<-group1[group1$`time point:ch1` %in% hour,] 
        gse73072_pre0_flu_gsm<-gse73072_pre0_gsm[gse73072_pre0_gsm$`virus:ch1` %in% virus,]
        gse73072_exp<-exprs1[,row.names(gse73072_pre0_flu_gsm)]   ##12023*13
       
        
         sym_index<-c(1,2,3,4,6,8,9,15,16,18,19,20,21)
        asym_index<-c(5,7,10,11,12,13,14,17)
        #shedding_index<-c(1,2,4,6,7,17,18,19,20,21)
        #nonshedding_index<-c(3,5,8,9,10,11,12,13,14,15,16)
        #sym_shedding_index<-c(1,2,4,18,19,20,21)
        #asym_nonshedding_index<-c(5,10,11,12,13,14)
        
        dee5_sym_index<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE5 subject ',sym_index)),]
        dee5_asym_index<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE5 subject ',asym_index)),]
        #dee5_shedding_index<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE5 subject ',shedding_index)),]
        #dee5_nonshedding_index<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE5 subject ',nonshedding_index)),]
        #dee5_sym_shedding_index<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE5 subject ',sym_shedding_index)),]
        #dee5_asym_nonshedding_index<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE5 subject ',asym_nonshedding_index)),]
        
        #subject_name_sym <- apply(data.frame(dee5_sym_shedding_index$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE5",sep = "_"))
        #subject_name_asym<- apply(data.frame(dee5_asym_nonshedding_index$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE5",sep = "_"))
        
         #dee5_group_sym <- data.frame(gseID=dee5_sym_shedding_index$geo_accession, subject=subject_name_sym,group="Symptomatic")
        #dee5_group_asym <- data.frame(gseID=dee5_asym_nonshedding_index$geo_accession, subject=subject_name_asym,group="Asymptomatic")
        
        subject_name_sym <- apply(data.frame(dee5_sym_index$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE5",sep = "_"))
        subject_name_asym<- apply(data.frame(dee5_asym_index$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE5",sep = "_"))
        
        dee5_group_sym <- data.frame(gseID=dee5_sym_index$geo_accession, subject=subject_name_sym,group="Symptomatic")
        dee5_group_asym <- data.frame(gseID=dee5_asym_index$geo_accession, subject=subject_name_asym,group="Asymptomatic")
        
        
        dee5_group<- rbind(dee5_group_sym,dee5_group_asym)
        rownames(dee5_group)<-dee5_group$gseID
        write.csv(dee5_group,"dee5_group.csv")
        ###calculate the mean value
        dee5_sym_mean_matrix<-as.data.frame(exprs1[,dee5_sym_index$geo_accession])
        dee5_asym_mean_matrix<-as.data.frame(exprs1[,dee5_asym_index$geo_accession])
       
        
        ###export the mean matrix
        write.csv(dee5_sym_mean_matrix,'dee5_sym_mean_matrix.csv',quote = F,row.names = T)
        write.csv(dee5_asym_mean_matrix,'dee5_asym_mean_matrix.csv',row.names = T)
      
  
  ##GSE61754 
  
    ###extract the non-vaccinate subject
    group4_nonva_index<-group4[group4$"vaccination status:ch1"=='Control',]
    ###extract the baseline subject
    group4_nonva_base_index<-group4_nonva_index[group4_nonva_index$"timepoint:ch1"=='Pre-challenge',]
    
    ###group by the symptom only and generate the expression matrix
    asym_index<-group4_nonva_base_index$"symptom severity:ch1"=="None"
    group4_asym_index<-group4_nonva_base_index[group4_nonva_base_index$"symptom severity:ch1"=="None",]
    group4_sym_index<-group4_nonva_base_index[asym_index=='FALSE',]
    
    gse61754_symonly<-exprs4[,group4_asym_index$geo_accession]
    gse61754_asymonly<-exprs4[,group4_sym_index$geo_accession]
    
    gse61754_group_sym  <- data.frame(gseID=rownames(group4_sym_index),subject="subject",group='Symptomatic')
    gse61754_group_asym  <- data.frame(gseID=rownames(group4_asym_index),subject="subject",group='Asymptomatic')
    gse61754_group <- rbind( gse61754_group_sym, gse61754_group_asym)
    rownames( gse61754_group)<- gse61754_group$gseID
    write.csv(gse61754_group,"gse61754_group.csv")
    write.csv(gse61754_symonly,"gse61754_sym.csv")
    write.csv(gse61754_asymonly,"gse61754_asym.csv")
    

  
