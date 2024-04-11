

##############
##############
####差异分析
#############
#加载limma包，用于校正和比较差异
##############
##############
####差异分析
#############
#加载limma包，用于校正和比较差异
rm(list = ls())
setwd("E://SomeProje/Significant expression genes  Project/part1")
library(limma)
#data <- read.table("par1exprSet_rmdup.Rdata",sep='\t',head=T)
load(file = "exprSet31240_rmdup.Rdata")
exprSet<-exprSet2
exprSet<-exprSet2[,c(1,3,5,7,2,4,6,8)]
#rownames(allDiff)
#Apcdd1<-allDiff[which(rownames(allDiff)=="Apcdd1"),]
#write.table(exprSet,"exprSet.txt",sep='\t',quote=F)
#differential差异分析
##1.构建分组矩阵,有两种方式
##https://dwz.cn/wR2qU7s9
#这一段有两种方法，但是初学的时候特别容易误解，这一步分完全就是独立的
## 创建分组
### 这一步根据样本来就行
group <- c(rep("DMSO",4),rep("CHIR",4)) 
## 分组变成向量，并且限定leves的顺序
## levels里面，把对照组放在前面
group <- factor(group,levels = c("DMSO","CHIR"),ordered = F)
## 构建比较矩阵
design <- model.matrix(~group)
## 比较矩阵命名
colnames(design) <- levels(group)
design

#2.线性模型拟合
fit <- lmFit(exprSet,design)
#3.贝叶斯检验
fit2 <- eBayes(fit)
#4.输出差异分析结果,其中coef的数目不能超过design的列数
# 此处的2代表的是第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
save(allDiff,file = "../0401/allDiff.Rda")
write.table(allDiff,"../0401/GSE131240.diff.txt",sep='\t',quote=F)
#找出差异两倍以上，pvalue小于0.05，
diffLab <- subset(allDiff,abs(logFC) > 5 & adj.P.Val < 1e-5)

###################################################################################
#library(dplyr)
# diffLab <- allDiff %>% 
#   filter(adj.P.Val < 0.05) %>% 
#   filter(abs(logFC) >1)

library(magrittr)
library(tibble)
diffLab <- allDiff %>%
  rownames_to_column() %>%
  dplyr::filter(adj.P.Val < 1e-5) %>% 
  dplyr::filter(abs(logFC) > 5) %>%
  column_to_rownames()
save(diffLab,group,file = "../0401/diffLab.Rda")

## 作图环节
####################################################################################
####################################################################################
## 1.把数据调整成可以作图的格式

exprSet <- t(exprSet)
## 矩阵转数据框
exprSet <- as.data.frame(exprSet)
## 增加一列为第一列
dd <- cbind(group=group,exprSet)

## 2.作图展示
library(ggplot2)
ggplot(data = dd,aes(x=group,y=CD36,fill=group))+
  geom_boxplot()+
  geom_point()+
  theme_bw()

## 3.steal plot
my_comparisons <- list(
  c("DMSO","CHIR")
)
library(ggpubr)
ggboxplot(
  dd, x = "group", y = "CD36",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("AggreDMSO24h", "AggreBIO24h")
  )
  library(ggpubr)
  ggboxplot(
    dd, x = "group", y = gene,
    color = "group", palette = c("#00AFBB", "#E7B800"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

diffplot("DUSP6")
diffplot("MOXD1")


## 4.多个基因作图查看
## 先把基因提取出来
genelist <- rownames(diffLab)
## 再提取表达量
data <- dd[,c("group",genelist)]
## 用gather调整数据
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
## 多基因作图
#################################
pdf("../0401/GSE131240_alldiff.pdf",10,6)
ggplot(data = data,aes(x=gene,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  stat_compare_means(aes(group=group), label = "p.signif", method = "t.test")
dev.off()
## 尝试更清晰的展示
pdf("../0401/formatGSE131240.pdf",30,25)
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~gene)+
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "t.test")
dev.off()
