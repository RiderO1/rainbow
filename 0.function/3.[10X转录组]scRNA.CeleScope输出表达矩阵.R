# TODO: 将celescope结果转化为seurat对象
# 
# Author: heshengyuan
####修改日志 - 2023-04-13
#' 修改了获取样本结果文件夹的方式
#' 原来根据样本名匹配，会受到样本命名方式的影响
#' 由于CeleScope结果文件夹下除了shell文件夹外，其他文件夹即为样本结果文件夹，改成通过匹配非shell文件夹外的文件夹获取样本结果
#'
##
###############################################################################

#' ######
#' @description 提取CeleScope输出的barcodes.tsv、features.tsv、matrix.mtx文件信息，转化为seurat对象
#' @param File.dir，character；CeleScope的结果文件路径，包含每个样本结果文件夹
#' @Return results，list；存储每个样本对应的Seurat对象
#' ######
Extract.CeleScope.scRNAseq.10x <- function(File.dir)
{
    library(Seurat)

    #获得样本名
	# sample.names <- dir(File.dir, pattern = "S\\d+")
	sample.names <- list.dirs(File.dir, recursive = FALSE, full.names = FALSE)
	sample.names <- sample.names[grep("shell", sample.names, invert = TRUE)]
    # 读取每个样本的单细胞数据（raw count）
    results = lapply(sample.names, function(i){
        #' Point [1]
        #' 内容：获取路径信息
        tmp.path = file.path(File.dir, i, "05.count")
        # 获取barcodes.tsv、features.tsv、matrix.mtx所在的文件夹路径
        File = list.files(tmp.path, pattern = "filtered_feature_bc_matrix", full.names = TRUE)

        #' Point [2]
        #' 内容：读取数据
        cat(format(Sys.time(), "[%b-%d,%H:%M:%S] "))
        cat("正在读取样本", i, "...\n")
        temp <- Read10X(File)	#读取重要的三个文件(barcodes.tsv、features.tsv、matrix.mtx)

        #' Point [3]
        #' 内容：转化成seurat对象
        # 建立样本名字meta
        t.meta = data.frame(sampleID = rep(i, ncol(temp)), stringsAsFactors = FALSE)
        rownames(t.meta) = colnames(temp)

        # 统一使用seurat对象存储原始count数据
        temp.seurat <- CreateSeuratObject(counts = temp, meta.data = t.meta)

        # 产生线粒体基因count占比信息，方便未来过滤使用
        temp.seurat[["percent.mt"]] = PercentageFeatureSet(temp.seurat, pattern = "^MT-")
        return(temp.seurat)
    })
    names(results) = sample.names
    return(results)
}


#######--2022年08月05日--#######
#' -----------测试案例-----------
# inputDir = "/home/xteam/Data/Projects/Bone/Results/Celescope"
# tmp = Extract.CeleScope.scRNAseq.10x(inputDir)
#' -----------结束测试-----------
#####################
