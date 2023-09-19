rm(list=ls())
library(multtest)
library(Seurat)
library(dplyr)
library(tidyverse)
getwd()
setwd("D:/单细胞/sc")
dir="D:/单细胞/sc"
list.files(dir)
sc<-Read10X(data.dir = dir) 
#读取说明文件
#meta_p<-read.table("D:\\单细胞\\SraRunTable.txt",sep=",",header = T)
#table(meta_p$Tissue)
#table(meta_p$Donor)
#table(meta_p$Site)
#table(meta_p$Site,meta$Donor)
#meta<-meta_p[,c(13,28)]
#meta$group<-paste0(substr((meta$Site),1,1),substr((meta$Donor),9,9))
#meta1<-str_replace(meta$group,c("P"),c("T"))
#meta$group<-meta1
sc<- CreateSeuratObject(counts = sc,min.cells = 3,min.features =200)
patient_id<-substr(rownames(sc@meta.data),18,18)
meta_sce<-data.frame(patient_id=patient_id,row.names=rownames(sc@meta.data))
sc<-Read10X(data.dir = dir) 
sc<- CreateSeuratObject(counts = sc,meta.data=meta_sce,min.cells = 3,min.features =200)
table(sc$patient_id)

ncol(sc)
nrow(sc)
sc@meta.data
