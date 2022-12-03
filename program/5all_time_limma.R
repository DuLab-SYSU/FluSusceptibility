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
  library(tidyverse)
  
  
}

#=========================[ set working directory ]==========================

setwd("/home/dulab/Documents/wrok/flu_paper/data/result2/result_new/")  

#======================[ load RMA normalized data from GEO ]========================
#############
virus<-c('H3N2')
#gse73072_pre0_gsm<-group1[group1$`time point:ch1` %in% hour,] 
#gse73072_pre0_flu_gsm<-group1[group1$`virus:ch1` %in% virus,]
hour<-c('hour 0','hour 5','hour 12','hour 22',' hour 29','hour 36','hour 46','hour 53','hour 60',
        'hour 70','hour 77','hour 84','hour 94',
        'hour 101','hour 108','hour 118','hour 125','hour 132','hour 142','hour 166')
gse73072_pre0_gsm<-group1[group1$`time point:ch1` %in% hour,] 
gse73072_pre0_flu_gsm<-gse73072_pre0_gsm[gse73072_pre0_gsm$`virus:ch1` %in% virus,]
gse73072_exp<-exprs1[,row.names(gse73072_pre0_flu_gsm)]   ##12023*139


sym_index<-c(1,5,6,7,8,10,12,13,15)
asym_index<-c(2,3,4,9,11,14,16,17)


sym_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE2 subject ',sym_index)),]
asym_gsm<-gse73072_pre0_flu_gsm[gse73072_pre0_flu_gsm$"subject:ch1" %in% c(paste0('H3N2 DEE2 subject ',asym_index)), ]

subject_name_sym <- apply(data.frame(sym_gsm$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE2",sep = "_"))
subject_name_asym<- apply(data.frame(asym_gsm$`subject:ch1`),1,function(x) paste(unlist(strsplit(x,split = " "))[3],unlist(strsplit(x,split = " "))[4],"DEE2",sep = "_"))

setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi") 


gse73072_group_sym <- data.frame(gseID=sym_gsm$geo_accession, Time=sym_gsm$`time point:ch1`,subject=subject_name_sym,group="Symptomatic")
gse73072_group_asym <- data.frame(gseID=asym_gsm$geo_accession,Time=asym_gsm$`time point:ch1`, subject=subject_name_asym,group="Asymptomatic")

gse73072_group<- rbind(gse73072_group_sym,gse73072_group_asym)
rownames(gse73072_group)<-gse73072_group$gseID
write.csv(gse73072_group,"gse73072_group_dee2_alltime.csv")



dee2_sym_mean_matrix0<-as.data.frame(exprs1[,sym_gsm$geo_accession])
dee2_asym_mean_matrix0<-as.data.frame(exprs1[,asym_gsm$geo_accession])

####
write.csv(dee2_sym_mean_matrix0,'dee2_sym_mean_matrix_dee2_alltime.csv',quote = F,row.names = T)
write.csv(dee2_asym_mean_matrix0,'dee2_asym_mean_matrix_dee2_alltime.csv',row.names = T)
#############################################################################################3
dee2_sym_sub<-read.csv('dee2_sym_mean_matrix_dee2_alltime.csv',row.names = 1)
dee2_asym_sub<-read.csv('dee2_asym_mean_matrix_dee2_alltime.csv',row.names = 1)
sym_shedding_exprs0<-as.data.frame(cbind(dee2_sym_sub,dee2_asym_sub))
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
dee2_sym_exprs<-probe1symbol(org.Hs.egSYMBOL,sym_shedding_exprs0)
write.csv(dee2_sym_exprs[,-ncol(dee2_sym_exprs)],"dee2_sym_exprs_alltime.csv",row.names = T)

################################################
setwd("/home/dulab/Documents/wrok/flu_paper/data/result_fi") 

dee2_group <- read.csv("gse73072_group_dee2_alltime.csv",row.names = 1)
dee2_group$batch <- "gse73072_dee2"

dee2_group <- dplyr::filter(dee2_group, Time != "Baseline") %>%
  mutate(type = str_replace(Time, "hour ", "h")) 




exprs_dee2<-read.csv('dee2_sym_exprs_alltime.csv',row.names = 1)   #72*8286
#####################################
dee2_group$time_lev <- paste0(dee2_group$group,".",dee2_group$type)

lev <- unique(dee2_group$time_lev)
f <- factor(dee2_group$time_lev, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev
corfit<-duplicateCorrelation(exprs_dee2,design,block = dee2_group$gseID)
corfit$consensus
fit <- lmFit(exprs_dee2, design,block = dee2_group$gseID,correlation=corfit$consensu)

con_pair <- dee2_group[dee2_group$group=="Symptomatic",]
con_pair <- setdiff(con_pair$time_lev,"Symptomatic.h0")
con_pair<-   paste0(con_pair,"-","Symptomatic.h0")

cont.sym <- makeContrasts(
  "Symptomatic.h5-Symptomatic.h0" ,  "Symptomatic.h12-Symptomatic.h0" , "Symptomatic.h22-Symptomatic.h0" ,
  "Symptomatic.h36-Symptomatic.h0" , "Symptomatic.h46-Symptomatic.h0" , "Symptomatic.h53-Symptomatic.h0" ,
   "Symptomatic.h60-Symptomatic.h0" , "Symptomatic.h70-Symptomatic.h0",  "Symptomatic.h77-Symptomatic.h0" ,
   "Symptomatic.h84-Symptomatic.h0"  ,"Symptomatic.h94-Symptomatic.h0" , "Symptomatic.h101-Symptomatic.h0",
  "Symptomatic.h108-Symptomatic.h0" ,"Symptomatic.h118-Symptomatic.h0" ,"Symptomatic.h125-Symptomatic.h0",
  "Symptomatic.h132-Symptomatic.h0" ,"Symptomatic.h142-Symptomatic.h0", "Symptomatic.h166-Symptomatic.h0",
   levels=design)
 fit2 <- contrasts.fit(fit, cont.sym)
 fit3 <- eBayes(fit2,trend = TRUE)
 DEG=topTable(fit3, adjust="BH")

 dt <- decideTests(fit3,method = "separate",adjust.method = "BH",p.value = 0.05,lfc = 0)
 head(dt)
 
 summary(dt)
 write.fit(fit3,dt,adjust="BH",file ="sym_all_deg_limma.txt")
 
 deg <- read.table("sym_all_deg_limma.txt")
 deg$gene <- rownames(deg)
 groupinfo <- deg[,c(76:94)]
 groupinfo_melt <- melt(groupinfo,id="gene")
 colnames(groupinfo_melt)<-c("gene","type","group")
 type <- gsub("Results.Symptomatic.","",groupinfo_melt$type)
 type <- gsub(".Symptomatic.h0","",type)
 groupinfo_melt$type<- type
 groupinfo_melt <- within(groupinfo_melt,{
   degtype<-NA
   degtype[group==-1]<-"DOWN"
   degtype[group==1]<-"UP"
   degtype[group==0]<-"NotSig"
   
 })
 
 pvalue <- deg[,c(38:55,94)]
 pvalue_melt <- melt(pvalue,id="gene")
 colnames(pvalue_melt)<-c("gene","type","pvalue")
 type <- gsub("P.value.Symptomatic.","",pvalue_melt$type)
 type <- gsub(".Symptomatic.h0","",type)
 pvalue_melt$type<- type
 
 padj <- deg[,c(56:73,94)]
 padj_melt <- melt(padj,id="gene")
 colnames(padj_melt)<-c("gene","type","padj_BH")
 type <- gsub("P.value.adj.Symptomatic.","",padj_melt$type)
 type <- gsub(".Symptomatic.h0","",type)
 padj_melt$type<- type
 
 
 lfc <-deg[,c(2:19,94)]
 
 lfc_melt <- melt(lfc,id="gene")
 colnames(lfc_melt)<-c("gene","type","lfc")
 type <- gsub("Coef.Symptomatic.","",lfc_melt$type)
 type <- gsub(".Symptomatic.h0","",type)
 lfc_melt$type<- type
 
 deg_sym <- cbind(groupinfo_melt[,c(1,2,4)],pvalue_melt[,3],padj_melt[,3],lfc_melt[,3])
 colnames(deg_sym) <- c("gene","type","degtype","pvalue","padj_BH","lfc")
 deg_sym$group <- "Symptomatic"
###############################################################################################
 
 ################################################################################################
 con_pair <- dee2_group[dee2_group$group=="Asymptomatic",]
 con_pair <- setdiff(con_pair$time_lev,"Asymptomatic.h0")
 con_pair<-   paste0(con_pair,"-","Asymptomatic.h0")
cont.asym<- makeContrasts(
  "Asymptomatic.h5-Asymptomatic.h0" ,  "Asymptomatic.h12-Asymptomatic.h0" , "Asymptomatic.h22-Asymptomatic.h0" ,
   "Asymptomatic.h36-Asymptomatic.h0",  "Asymptomatic.h46-Asymptomatic.h0"  ,"Asymptomatic.h53-Asymptomatic.h0" ,
  "Asymptomatic.h60-Asymptomatic.h0",  "Asymptomatic.h70-Asymptomatic.h0"  ,"Asymptomatic.h77-Asymptomatic.h0" ,
  "Asymptomatic.h84-Asymptomatic.h0",  "Asymptomatic.h94-Asymptomatic.h0"  ,"Asymptomatic.h101-Asymptomatic.h0",
  "Asymptomatic.h108-Asymptomatic.h0", "Asymptomatic.h118-Asymptomatic.h0" ,"Asymptomatic.h125-Asymptomatic.h0",
  "Asymptomatic.h132-Asymptomatic.h0" ,"Asymptomatic.h142-Asymptomatic.h0", "Asymptomatic.h166-Asymptomatic.h0",
   levels=design)
fit2 <- contrasts.fit(fit, cont.asym)
fit3 <- eBayes(fit2,trend = TRUE)
DEG=topTable(fit3, adjust="BH")

dt <- decideTests(fit3,method = "separate",adjust.method = "BH",p.value = 0.05,lfc = 0)
head(dt)

summary(dt)
write.fit(fit3,dt,adjust="BH",file = "asym_all_deg_limma.txt")
#######################################################
deg <- read.table("asym_all_deg_limma.txt")
deg$gene <- rownames(deg)
groupinfo <- deg[,c(76:94)]
groupinfo_melt <- melt(groupinfo,id="gene")
colnames(groupinfo_melt)<-c("gene","type","group")
type <- gsub("Results.Asymptomatic.","",groupinfo_melt$type)
type <- gsub(".Asymptomatic.h0","",type)
groupinfo_melt$type<- type
groupinfo_melt <- within(groupinfo_melt,{
  degtype<-NA
  degtype[group==-1]<-"DOWN"
  degtype[group==1]<-"UP"
  degtype[group==0]<-"NotSig"
  
})

pvalue <- deg[,c(38:55,94)]
pvalue_melt <- melt(pvalue,id="gene")
colnames(pvalue_melt)<-c("gene","type","pvalue")
type <- gsub("P.value.Asymptomatic.","",pvalue_melt$type)
type <- gsub(".Asymptomatic.h0","",type)
pvalue_melt$type<- type

padj <- deg[,c(56:73,94)]
padj_melt <- melt(padj,id="gene")
colnames(padj_melt)<-c("gene","type","padj_BH")
type <- gsub("P.value.adj.Asymptomatic.","",padj_melt$type)
type <- gsub(".Asymptomatic.h0","",type)
padj_melt$type<- type


lfc <-deg[,c(2:19,94)]

lfc_melt <- melt(lfc,id="gene")
colnames(lfc_melt)<-c("gene","type","lfc")
type <- gsub("Coef.Asymptomatic.","",lfc_melt$type)
type <- gsub(".Asymptomatic.h0","",type)
lfc_melt$type<- type

deg_asym <- cbind(groupinfo_melt[,c(1,2,4)],pvalue_melt[,3],padj_melt[,3],lfc_melt[,3])
colnames(deg_asym) <- c("gene","type","degtype","pvalue","padj_BH","lfc")
deg_asym$group <- "Asymptomatic"

deg_all <- rbind(deg_sym,deg_asym)
write.csv(deg_all,"deg_limma_alltime_summary.csv")
#deg_melt <- melt(deg,id="gene")
################################################################################################################
############################################################################################################

cont.sym <- makeContrasts(
  "Symptomatic.h5-Symptomatic.h0" ,  "Symptomatic.h12-Symptomatic.h5" , "Symptomatic.h22-Symptomatic.h12" ,
  "Symptomatic.h36-Symptomatic.h22" , "Symptomatic.h46-Symptomatic.h36" , "Symptomatic.h53-Symptomatic.h46" ,
  "Symptomatic.h60-Symptomatic.h53" , "Symptomatic.h70-Symptomatic.h60",  "Symptomatic.h77-Symptomatic.h70" ,
  "Symptomatic.h84-Symptomatic.h77"  ,"Symptomatic.h94-Symptomatic.h84" , "Symptomatic.h101-Symptomatic.h94",
  "Symptomatic.h108-Symptomatic.h101" ,"Symptomatic.h118-Symptomatic.h108" ,"Symptomatic.h125-Symptomatic.h118",
  "Symptomatic.h132-Symptomatic.h125" ,"Symptomatic.h142-Symptomatic.h132", "Symptomatic.h166-Symptomatic.h142",
  levels=design)
fit2 <- contrasts.fit(fit, cont.sym)
fit2 <- eBayes(fit2)
DEG=topTable(fit2, adjust="BH")

DEG = na.omit(DEG)
DEG$change = ifelse(DEG$F>0&DEG$adj.P.Val<=0.05,"UP",ifelse(DEG$F<0&DEG$adj.P.Val<=0.05,"DOWN","NOT"))
table(DEG$change)

