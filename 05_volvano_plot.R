

###############
###############
##volcano火山图
###############
###############
rm(list = ls())
setwd("D://SomeProje/Significant expression genes  Project/GSE131240")
##用ggplot2
library(ggplot2)
library(ggrepel)
library(dplyr)
load(file = "allDiff.Rda")
data <- allDiff
data$gene <- rownames(data)
## 仔细观察data数据
## logFC，P.Value，gene
pdf("volcano.pdf",20,10)
ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val))) +
  ## 三个部分分别画点
  geom_point(data=subset(data,abs(data$logFC) >= 1),aes(size=abs(logFC)),color="black",alpha=0.3) +
  geom_point(data=subset(data,data$adj.P.Val < 1e-3 & data$logFC > 4),aes(size=abs(logFC)),color="red",alpha=0.4) +
  geom_point(data=subset(data,data$adj.P.Val < 1e-3 & data$logFC < -4),aes(size=abs(logFC)),color="green",alpha=0.4) +
  ## 画线
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  ## 主题
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+
  ## 标签
  geom_text_repel(data=subset(data, abs(logFC) > 5 ), aes(label=gene),col="black",alpha = 0.8)
dev.off()


##另一种风格
# library(ggplot2)
# library(ggrepel)
# data <- allDiff
# data$gene <- rownames(data)
# data$significant <- as.factor(data$adj.P.Val<1e-4 & abs(data$logFC) > 5)
# data$gene <- rownames(data)
# pdf("volcano.pdf",20,10)
# ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
#   geom_point(alpha=0.8, size=1.2,col="black")+
#   geom_point(data=subset(data, logFC > 1),alpha=0.8, size=1.2,col="red")+
#   geom_point(data=subset(data, logFC < -1),alpha=0.6, size=1.2,col="blue")+
#   labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
#   theme(plot.title = element_text(hjust = 0.4))+
#   geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
#   geom_vline(xintercept = c(0.5,-0.5),lty=4,lwd=0.6,alpha=0.8)+
#   theme_bw()+
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),   
#         axis.line = element_line(colour = "black")) +
#   geom_point(data=subset(data, abs(logFC) >= 6),alpha=0.8, size=3,col="green")+
#   geom_text_repel(data=subset(data, abs(logFC) > 6), 
#                   aes(label=gene),col="black",alpha = 0.8)
# 
# dev.off()



## 加载R
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="volcano.pptx")

## GEO教程长期更新的链接:
## https://codingsoeasy.com/archives/geo