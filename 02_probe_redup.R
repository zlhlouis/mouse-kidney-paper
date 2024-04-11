

rm(list = ls())
load(file = "exprSet31240_readGSE.Rdata")
boxplot(exprSet,outline=FALSE, notch=T, las=2)

## 自动log化，解释
## 注意此时数据的格式是矩阵
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { 
  #ex[which(ex <= 0)] <- NaN
  ## 取log2
  exprSet <- log2(ex)
  print("log2 transform finished")
  }else{
    print("log2 transform not needed")
  }

library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
## 这步把矩阵转换为很重要
exprSet <- as.data.frame(exprSet)

#################################################################
## 探针基因名转换

## 安装R包

options(BioC_mirror="http://bioconductor.org")
if(!require("mouse430a2.db")) BiocManager::install("mouse430a2.db",update = F,ask = F)
library(mouse430a2.db)
#获取探针
probe2symbol_df <- toTable(get("mouse430a2SYMBOL"))

#save(probe2symbol_df,file = "probe2symbol_df.Rdata")
#load(file = "probe2symbol_df.Rdata")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("clusterProfiler")
library("org.Mm.eg.db")
library("clusterProfiler")
rownames(exprSet)<-gsub("\\..*","",rownames(exprSet))
probe2symbol_df <- bitr(geneID = rownames(exprSet),fromType="ENSEMBL",toType=c("SYMBOL"),OrgDb=org.Mm.eg.db)
colnames(probe2symbol_df)<-c("probe_id","symbol")
id<-intersect(ENTREZID$ENSEMBL,rownames(exprSet))

exprSet1<-exprSet[pmatch(id,rownames(exprSet)),]
ENTREZID1 <- ENTREZID[pmatch(id,ENTREZID[,1]),]
exprSet2<-cbind(exprSet1,ENTREZID1)
exprSet3<-exprSet2[,c(-9,-10)]

#看一下symbol有没有重复，发现只有18834个，
length(unique(probe2symbol_df$symbol))
## 而探针和基因的对应关系要更多
nrow(probe2symbol_df)
# 所以需要去重，多个探针对应一个基因

###探针转换以及去重，获得最终的表达矩阵
library(dplyr)
library(tibble)
exprSet <- exprSet %>% 
  ## 行名转列名
  rownames_to_column("probe_id") %>% 
  #合并探针的信息
  inner_join(probe2symbol_df,by="probe_id") %>% 
  #去掉多余信息
  select(-probe_id) %>%  
  #重新排列
  select(symbol,everything()) %>%  
  #求出平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  # 去重，symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  #反向选择去除rowMean这一列
  select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("symbol")
 exprSet[exprSet == 0] = 0
cf_DFinf2NA <- function(x)
{
  for (i in 1:ncol(x)){
    x[,i][is.infinite(x[,i])] = NA
  }
  return(x)
}
exprSet1 <- cf_DFinf2NA(exprSet)
exprSet2<-na.omit(exprSet1)
write.table(exprSet2,"GSE31240_expression.profile",sep='\t',quote=F)
save(exprSet2,file = "exprSet31240_rmdup.Rdata")
