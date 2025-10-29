# TODO: 从整体的聚类结果中的选择特定的cluster重新聚类
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

#数据加载
scRNAseq.profile.obj <- readRDS("/IData/WorkSpace/heshengyuan/Bone/RData/clustered[Clinical].rds")
scRNAseq.profile.obj@meta.data$new.cluster = as.character(scRNAseq.profile.obj@meta.data$cell.cluster.1)

#工作路径设置
setwd("/IData/WorkSpace/heshengyuan/Bone/Results") 

#变量设置
outDir = "/IData/WorkSpace/heshengyuan/Bone/Results"
# outDir = "/IData/TMP/heshengyuan/tmpresult/SingleCell"
if(!file.exists(outDir)) 
	dir.create(outDir)
# marker表格存储路径
t.outDir = paste(outDir, "Table", "cluster再调整", "cluster_marker", sep = "/")
if(!file.exists(t.outDir)) 
	dir.create(t.outDir, recursive = TRUE)
if(!file.exists(file.path(outDir, "Figure", "cluster再调整")))
    dir.create(file.path(outDir, "Figure", "cluster再调整"), recursive = TRUE)

#'-Part [1]-
#' 内容：重新聚类Macrophage
#'#######
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="Macrophage")
#' Point [1]
#'内容：识别高变异基因并绘制手肘图
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/ClusteringSeurat.R")
t.seurat.reduced = Clustering.ndims.Elbow(t.seurat.reduced, de.novo = TRUE, out.path = file.path(outDir, "Figure", "cluster再调整", "ElbowPlot[Macrophage].pdf"))

#' Point [2]
#'内容：细胞聚类
t.seurat.reduced = ClusteringSeurat(t.seurat.reduced, de.novo = FALSE, ndim = 20, resolution = 0.4)
t.seurat.reduced@meta.data$clusterid = paste("Macrophage", as.character(t.seurat.reduced@meta.data[, "integrated_snn_res.0.4"]), sep = ".")

#' Point [3]
#'内容：将Macrophage的聚类结果整合到总得 seurat对象中
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/整合-融合seurat的meta信息.R")
scRNAseq.profile.obj = MergeCluster(seurat.object = scRNAseq.profile.obj, sub.seurat.object = t.seurat.reduced, from = "clusterid", to = "new.cluster")



#'-Part [2]-
#' 内容：重新聚类 T cell
#'#######
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="T cell")
#' Point [1]
#'内容：识别高变异基因并绘制手肘图
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/ClusteringSeurat.R")
t.seurat.reduced = Clustering.ndims.Elbow(t.seurat.reduced, de.novo = TRUE, out.path = file.path(outDir, "Figure", "cluster再调整", "ElbowPlot[T cell].pdf"))

#' Point [2]
#'内容：细胞聚类
t.seurat.reduced = ClusteringSeurat(t.seurat.reduced, de.novo = FALSE, ndim = 20, resolution = 0.3)
t.seurat.reduced@meta.data$clusterid = paste("T cell", as.character(t.seurat.reduced@meta.data[, "integrated_snn_res.0.3"]), sep = ".")

#' Point [3]
#'内容：将T cell的聚类结果整合到总得 seurat对象中
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/整合-融合seurat的meta信息.R")
scRNAseq.profile.obj = MergeCluster(seurat.object = scRNAseq.profile.obj, sub.seurat.object = t.seurat.reduced, from = "clusterid", to = "new.cluster")


#'-Part [3]-
#' 内容：重新聚类 Fibroblast
#'#######
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="Fibroblast")
#' Point [1]
#'内容：识别高变异基因并绘制手肘图
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/ClusteringSeurat.R")
t.seurat.reduced = Clustering.ndims.Elbow(t.seurat.reduced, de.novo = TRUE, out.path = file.path(outDir, "Figure", "cluster再调整", "ElbowPlot[T cell].pdf"))

#' Point [2]
#'内容：细胞聚类
t.seurat.reduced = ClusteringSeurat(t.seurat.reduced, de.novo = FALSE, ndim = 20, resolution = 0.4)
t.seurat.reduced@meta.data$clusterid = paste("Fibroblast", as.character(t.seurat.reduced@meta.data[, "integrated_snn_res.0.4"]), sep = ".")

#' Point [3]
#'内容：将T cell的聚类结果整合到总得 seurat对象中
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/整合-融合seurat的meta信息.R")
scRNAseq.profile.obj = MergeCluster(seurat.object = scRNAseq.profile.obj, sub.seurat.object = t.seurat.reduced, from = "clusterid", to = "new.cluster")


#'-Part [4]-
#' 内容：重新聚类 Endothelial cell
#'#######
t.seurat.reduced = subset(scRNAseq.profile.obj, subset = cell.cluster.1=="Endothelial cell")
#' Point [1]
#'内容：识别高变异基因并绘制手肘图
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/ClusteringSeurat.R")
t.seurat.reduced = Clustering.ndims.Elbow(t.seurat.reduced, de.novo = TRUE, out.path = file.path(outDir, "Figure", "cluster再调整", "ElbowPlot[T cell].pdf"))

#' Point [2]
#'内容：细胞聚类
t.seurat.reduced = ClusteringSeurat(t.seurat.reduced, de.novo = FALSE, ndim = 20, resolution = 0.2)
t.seurat.reduced@meta.data$clusterid = paste("Endothelial cell", as.character(t.seurat.reduced@meta.data[, "integrated_snn_res.0.2"]), sep = ".")

#' Point [3]
#'内容：将T cell的聚类结果整合到总得 seurat对象中
source("/pub5/xiaoyun/BioY/heshengyuan/somemission/0.单细胞数据处理/1.头颈癌单细胞数据处理/0.function/3.cluster再调整/整合-融合seurat的meta信息.R")
scRNAseq.profile.obj = MergeCluster(seurat.object = scRNAseq.profile.obj, sub.seurat.object = t.seurat.reduced, from = "clusterid", to = "new.cluster")

saveRDS(scRNAseq.profile.obj, file = paste("/IData/WorkSpace/heshengyuan/Bone/RData", "clustered[new].RDS", sep="/"))
