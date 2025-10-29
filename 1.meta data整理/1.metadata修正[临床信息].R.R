# TODO: 添加临床信息信息到seurat对象
# 
# Author: heshengyuan
###############################################################################

#'-Part [0]-
#' 内容：前期准备，路径设置和函数调用
#'#######
# R包加载
library(Seurat)
library(ggsci)
library(ggplot2)
library(openxlsx)
library(reshape2)

#工作路径设置
setwd("/IData/WorkSpace/heshengyuan/Bone/Results") 

#变量设置
outDir = "/IData/WorkSpace/heshengyuan/Bone/RData"
RDS.name = "clustered[Clinical].rds" # 根据注释信息分类后保存的seurat对象名

#数据加载
scRNAseq.profile.obj <- readRDS("/IData/WorkSpace/heshengyuan/Bone/RData/clustered[CellTypes].rds")


#'-Part [1]-
#' 内容：添加临床信息到Seurat对象
#'#######
# 读取临床信息
clinicalInfo = read.xlsx("/IData/WorkSpace/heshengyuan/Bone/clinical_information.xlsx")
rownames(clinicalInfo) = clinicalInfo[, "sampleID"] # 更改行名

# 添加临床信息到meta data
scRNAseq.profile.obj@meta.data = cbind(scRNAseq.profile.obj@meta.data, clinicalInfo[scRNAseq.profile.obj@meta.data$sampleID, -1])
saveRDS(scRNAseq.profile.obj, file = paste(outDir, RDS.name, sep="/")) # 保存数据
