

###############
####heatmap热图
#用行名提取数据
rm(list = ls())
## 加载表达数据
setwd("D://SomeProje/Significant expression genes  Project/0401/")
load(file = "exprSet31240_rmdup.Rdata")
exprSet<-exprSet2
## 加载差异列表
load(file = "allDiff.Rda")
library(dplyr)
library(tibble)
diffLab <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 1e-4) %>% 
  filter(abs(logFC) > 4) %>% 
  column_to_rownames()

#write.table(diffLab,"/Users/louis/Documents/2022/差异表达基因/part1/GSE131240_p0.05FC1.txt",sep = '\t',quote=F,row.names = T)
## 加载热图的R包
library(pheatmap)
##用名称提取部分数据用作热图绘制
heatdata <- exprSet[rownames(diffLab),]
heatdata <- heatdata[,c(1,3,5,7,2,4,6,8)]
##制作一个分组信息用于注释
group <- c(rep("DMSO",4),rep("CHIR",4)) 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)
#annotation_col$group <- factor(annotation_col$group,levels = c("con","treat"),ordered = F)

## 直接作图
pdf("GSE131240_heapmap.pdf",10,14)
pheatmap(heatdata,cluster_rows= T,annotation_col =annotation_col,cellwidth = 12, cellheight = 12)
dev.off()
#如果注释出界，可以通过调整格子比例和字体修正
# 
# pheatmap(heatdata, #热图的数据
#          cluster_rows = TRUE,#行聚类
#          cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
#          annotation_col =annotation_col, #标注样本分类
#          annotation_legend=TRUE, # 显示注释
#          show_rownames = F,# 显示行名
#          scale = "row", #以行来标准化，这个功能很不错
#          color =colorRampPalette(c("blue", "white","red"))(100),#调色
#          #filename = "heatmap_F.pdf",#是否保存
#          cellwidth = 12, cellheight = 12,# 格子比例
#          fontsize = 10)
# 
# ### 可以重新筛选，阈值设大一点
# diffLab <- allDiff %>% 
#   rownames_to_column() %>% 
#   filter(adj.P.Val < 0.2) %>% 
#   filter(abs(logFC) >1) %>% 
#   column_to_rownames()
# heatdata <- exprSet[rownames(diffLab),]
##制作一个分组信息用于注释
group <- c(rep("con",4),rep("treat",4)) 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pdf("new_GSE131240_heapmap.pdf",15,15)
 pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 20, cellheight = 12,# 格子比例
         fontsize = 10)
dev.off()


## 加载R包
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="heatmap.pptx")
## 记住保存
