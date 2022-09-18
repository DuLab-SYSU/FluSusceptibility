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
  #library(xlsx)
}

#=========================[ set working directory ]==========================

  setwd("/home/dulab/Documents/wrok/flu_paper/data/result2/result_new/")  
  
#======================[ load RMA normalized data from GEO ]========================
{
  ##load the data
    #GSE17156
    load("GSE17156_eSet.Rdata")
    exprs1<-exprs(gset[[1]])#raw_data/GSE17156/data/
    group1<-pData(gset[[1]])
    #GSE30550
    load("GSE30550_eSet.Rdata")#raw_data/gse30550/data/
    exprs2<-exprs(gset[[1]])
    group2<-pData(gset[[1]])
    #GSE61754
    load("GSE61754_eSet.Rdata")#raw_data/gse61754/data/
    exprs3<-exprs(gset[[1]])
    group3<-pData(gset[[1]])
}


#=============================[ data pre-processing ]=======================

  #================[ extract the baseline time point,the influenza subject,group by symptoms/(symptoms and shedding) and it's expression matrix ]===============
  ##GSE73072
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_newest/")  
  ##GSE17156 :baseline only
  
    ###extract the flu subject
    group1_flu_index<-group1[80:113,]
    ###extract the baseline time point
    group1_flu_base_index<-group1_flu_index[group1_flu_index$"timepoint:ch1"=='baseline',]
    ###grouping by the symptom and shedding
    gse17156_base_sym_index<-group1_flu_base_index[group1_flu_base_index$"symptom group:ch1"=='Symptomatic',]
    gse17156_base_asym_index<-group1_flu_base_index[group1_flu_base_index$"symptom group:ch1"=='Asymptomatic',]
    ###develop expression matrix
    gse17156_sym<-exprs1[,gse17156_base_sym_index$geo_accession]
    gse17156_asym<-exprs1[,gse17156_base_asym_index$geo_accession]
    ###export the expression matrix
    write.csv(gse17156_sym,'gse17156_sym.csv')
    write.csv(gse17156_asym,'gse17156_asym.csv')
    
    gse17156_group_sym <- data.frame(gseID=rownames(gse17156_base_sym_index),subject=rownames(gse17156_base_sym_index),group="Symptomatic")
    gse17156_group_asym <- data.frame(gseID=rownames(gse17156_base_asym_index),subject=rownames(gse17156_base_asym_index),group="Asymptomatic")
    gse17156_group <- rbind(gse17156_group_asym,gse17156_group_sym)
    rownames(gse17156_group)<-gse17156_group$gseID
    ###export the expression matrix
    write.csv(gse17156_group,'gse17156_group.csv')
  
  
  ##GSE30550:two time point: -24 ,0 
 
    ###extract the baseline(-24) and 0h
    group2_base_index<-group2[group2$"time_hpi:ch1" %in% c('Hour 00'),] #,"Baseline"
    ###grouping by the symptom and shedding
    group2_sym_index<-group2_base_index[group2_base_index$"clinic_pheno:ch1"=='Symptomatic',]
    group2_asym_index<-group2_base_index[group2_base_index$"clinic_pheno:ch1"=='Asymptomatic',]
    
    
    subject_name_sym <- apply(data.frame(group2_sym_index$source_name_ch1),1,function(x) unlist(strsplit(paste(unlist(strsplit(x,split = " "))[4],unlist(strsplit(x,split = " "))[5],sep = ""),split = ","))[1])
    subject_name_asym<- apply(data.frame(group2_asym_index$source_name_ch1),1,function(x) unlist(strsplit(paste(unlist(strsplit(x,split = " "))[4],unlist(strsplit(x,split = " "))[5],sep = ""),split = ","))[1])
    gse30550_group_sym  <- data.frame(gseID=rownames(group2_sym_index),subject=subject_name_sym,group='Symptomatic')
    gse30550_group_asym  <- data.frame(gseID=rownames(group2_asym_index),subject=subject_name_asym,group='Asymptomatic')
    gse30550_group <- rbind( gse30550_group_sym, gse30550_group_asym)
    rownames( gse30550_group)<- gse30550_group$gseID
    write.csv(gse30550_group,"gse30550_group.csv")
    
    ###generate  matrix
    gse30550_sym<-as.data.frame(t(exprs2[,group2_sym_index$geo_accession]))
    gse30550_asym<-as.data.frame(t(exprs2[,group2_asym_index$geo_accession]))
    
    #sym_sub<-c(rep(1:7,times=c(2,2,2,2,2,2,2)),8,9,9)
    #asym_sub<-c(rep(10:17,times=c(2,2,2,2,2,2,2,2)))
    gse30550_sym_sub<-cbind(subject=subject_name_sym,gse30550_sym)
    gse30550_asym_sub<-cbind(subject=subject_name_asym,gse30550_asym)
    ###export the matrix
    write.csv(gse30550_sym_sub,'gse30550_sym_matrix.csv',row.names = FALSE)
    write.csv(gse30550_asym_sub,'gse30550_asym_matrix.csv',row.names = FALSE)
  
  
  ##GSE61754 
  
    ###extract the non-vaccinate subject
    group3_nonva_index<-group3[group3$"vaccination status:ch1"=='Control',]
    ###extract the baseline subject
    group3_nonva_base_index<-group3_nonva_index[group3_nonva_index$"timepoint:ch1"=='Pre-challenge',]
    
    ###group by the symptom only and generate the expression matrix
    asym_index<-group3_nonva_base_index$"symptom severity:ch1"=="None"
    group3_asym_index<-group3_nonva_base_index[group3_nonva_base_index$"symptom severity:ch1"=="None",]
    group3_sym_index<-group3_nonva_base_index[asym_index=='FALSE',]
    
    gse61754_symonly<-exprs3[,group3_asym_index$geo_accession]
    gse61754_asymonly<-exprs3[,group3_sym_index$geo_accession]
    
    gse61754_group_sym  <- data.frame(gseID=rownames(group3_sym_index),subject="subject",group='Symptomatic')
    gse61754_group_asym  <- data.frame(gseID=rownames(group3_asym_index),subject="subject",group='Asymptomatic')
    gse61754_group <- rbind( gse61754_group_sym, gse61754_group_asym)
    rownames( gse61754_group)<- gse61754_group$gseID
    write.csv(gse61754_group,"gse61754_group.csv")
    #write.csv(gse61754_symonly,"gse61754_sym.csv")
    #write.csv(gse61754_asymonly,"gse61754_asym.csv")
    
    ###grouping by the shedding only and generate the expression matrix
    group3_shedding_index<-group3_nonva_base_index[group3_nonva_base_index$"viral shedding:ch1"=='Positive',]
    group3_nonshedding_index<-group3_nonva_base_index[group3_nonva_base_index$"viral shedding:ch1"=='Negative',]
    
    gse61754_shedding<-exprs3[,group3_shedding_index$geo_accession]
    gse61754_nonshedding<-exprs3[,group3_nonshedding_index$geo_accession]
    
    ###group by the symptom and shedgse61754_symonlyding and generate the expression matrix
    group3_sym_shedding_index<-group3_sym_index[group3_sym_index$"viral shedding:ch1"=='Positive',]
    group3_asym_nonshedding_index<-group3_asym_index[group3_asym_index$"viral shedding:ch1"=='Negative',]
    gse61754_group_sym  <- data.frame(gseID=rownames(group3_sym_shedding_index),subject=rownames(group3_sym_shedding_index),group='Symptomatic')
    gse61754_group_asym  <- data.frame(gseID=rownames(group3_asym_nonshedding_index),subject=rownames(group3_asym_nonshedding_index),group='Asymptomatic')
    gse61754_group <- rbind( gse61754_group_sym, gse61754_group_asym)
    rownames( gse61754_group)<- gse61754_group$gseID
    write.csv(gse61754_group,"gse61754_group_asym_nonshedding.csv")
    
    gse61754_sym_shedding<-exprs3[,group3_sym_shedding_index$geo_accession]
    gse61754_asym_nonshedding<-exprs3[,group3_asym_nonshedding_index$geo_accession]
    
    ###export the expression matrix together
    #write.csv(gse61754_symonly,'gse61754_symonly.csv')
    #write.csv(gse61754_asymonly,'gse61754_asymonly.csv')
    #write.csv(gse61754_shedding,'gse61754_shedding.csv')
    #write.csv(gse61754_nonshedding,'gse61754_nonshedding.csv')
    write.csv(gse61754_sym_shedding,'gse61754_sym_shedding.csv')
    write.csv(gse61754_asym_nonshedding,'gse61754_asym_nonshedding.csv')
  

