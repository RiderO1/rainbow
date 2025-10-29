# TODO: 提取重聚类后子类细胞的top30marker基因
# 
# Author: heshengyuan
###############################################################################


#'-Part [0]-
#' 内容：前期准备，路径设置和函数调用
#'#######
# R包加载
library(Seurat)
library(tidyverse)

#工作路径设置
setwd("/IData/WorkSpace/heshengyuan/Bone/Results") 

#变量设置
outDir = "/IData/WorkSpace/heshengyuan/Bone/Results"
t.outDir = paste(outDir, "Table", "cluster再调整", "cluster_marker", sep = "/")

#数据加载
scRNAseq.profile.obj <- readRDS("/IData/WorkSpace/heshengyuan/Bone/RData/clustered[new].RDS")


#' Point [1]
#'内容：巨噬细胞
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="Macrophage")
# 识别marker
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/识别所有Cluster的markers.R") # FindAllMarkers.ClusterCell()
markers <- FindAllMarkers.ClusterCell(t.seurat.reduced, Idents.col="new.cluster")
# 提取top30 marker
tmpmarker <- dplyr::group_by(markers, cluster) %>% 
        dplyr::slice_max(order_by = avg_log2FC, n = 30)
marker.list = lapply(unique(tmpmarker$cluster),function(x){
    tmp = pull(dplyr::filter(tmpmarker, cluster == x),gene)
    res = paste(x, paste0(tmp, collapse = ","),sep = ":")
    return(res)
})
marker.df = do.call("rbind",marker.list)
write.table(marker.df, file = paste(t.outDir, "cluster marker[Macrophage].txt", sep = "/"),
            row.names = FALSE, col.names = FALSE,quote = FALSE) # 将marker输出到表格中

#' Point [2]
#'内容：T细胞
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="T cell")
# 识别marker
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/识别所有Cluster的markers.R") # FindAllMarkers.ClusterCell()
markers <- FindAllMarkers.ClusterCell(t.seurat.reduced, Idents.col="new.cluster")
# 提取top30 marker
tmpmarker <- dplyr::group_by(markers, cluster) %>% 
        dplyr::slice_max(order_by = avg_log2FC, n = 30)
marker.list = lapply(unique(tmpmarker$cluster),function(x){
    tmp = pull(dplyr::filter(tmpmarker, cluster == x),gene)
    res = paste(x, paste0(tmp, collapse = ","),sep = ":")
    return(res)
})
marker.df = do.call("rbind",marker.list)
write.table(marker.df, file = paste(t.outDir, "cluster marker[T cell].txt", sep = "/"),
            row.names = FALSE, col.names = FALSE,quote = FALSE) # 将marker输出到表格中

#' Point [3]
#'内容：成纤维细胞
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="Fibroblast")
# 识别marker
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/识别所有Cluster的markers.R") # FindAllMarkers.ClusterCell()
markers <- FindAllMarkers.ClusterCell(t.seurat.reduced, Idents.col="new.cluster")
# 提取top30 marker
tmpmarker <- dplyr::group_by(markers, cluster) %>% 
        dplyr::slice_max(order_by = avg_log2FC, n = 30)
marker.list = lapply(unique(tmpmarker$cluster),function(x){
    tmp = pull(dplyr::filter(tmpmarker, cluster == x),gene)
    res = paste(x, paste0(tmp, collapse = ","),sep = ":")
    return(res)
})
marker.df = do.call("rbind",marker.list)
write.table(marker.df, file = paste(t.outDir, "cluster marker[Fibroblast].txt", sep = "/"),
            row.names = FALSE, col.names = FALSE,quote = FALSE) # 将marker输出到表格中

#' Point [4]
#'内容：内皮细胞
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="Endothelial cell")
# 识别marker
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/识别所有Cluster的markers.R") # FindAllMarkers.ClusterCell()
markers <- FindAllMarkers.ClusterCell(t.seurat.reduced, Idents.col="new.cluster")
# 提取top30 marker
tmpmarker <- dplyr::group_by(markers, cluster) %>% 
        dplyr::slice_max(order_by = avg_log2FC, n = 30)
marker.list = lapply(unique(tmpmarker$cluster),function(x){
    tmp = pull(dplyr::filter(tmpmarker, cluster == x),gene)
    res = paste(x, paste0(tmp, collapse = ","),sep = ":")
    return(res)
})
marker.df = do.call("rbind",marker.list)
write.table(marker.df, file = paste(t.outDir, "cluster marker[Endothelial cell].txt", sep = "/"),
            row.names = FALSE, col.names = FALSE,quote = FALSE) # 将marker输出到表格中

