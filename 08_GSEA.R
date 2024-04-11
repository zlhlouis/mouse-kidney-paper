
#############
#############
rm(list = ls())
setwd("D://SomeProje/Significant expression genes  Project/GSE131240")
library(clusterProfiler)
#load(file = "par2allDiff.Rda")
load(file = "allDiff.Rda")
#获得基因列表
gene <- rownames(allDiff)
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
#save(gene,file="gene_GSEA.Rdata")
## load(file="gene_GSEA.Rdata")
gene_df <- data.frame(logFC=allDiff$logFC,
                      SYMBOL = rownames(allDiff))
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$SYMBOL
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
#library(clusterProfiler)

################################################################
#################################################################
### GSEA 变化在于gene set
#################################################################
### 1.hallmarks gene set
## 读入hallmarks gene set？
library(msigdbr)
hallmarks <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_t2g = hallmarks %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
#enricher(gene = gene_symbols_vector, TERM2GENE = msigdbr_t2g, ...)

### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =msigdbr_t2g,nPerm = 10000)

yd <- data.frame(y)
write.table(yd,"gsea.pro",sep='\t',quote=F)
### 看整体分布
#dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)

######################### 点图
pdf("gseadotplot.pdf",25,12)
library(stringi)
library(ggplot2)
dotplot(y,showCategory=24,split=".sign")+
  facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) stri_sub(x,10))
dev.off()
#########################折线图
library(enrichplot)
pdf("gsea.pdf",20,15)
gseaplot2(y, geneSetID = 1:10)
dev.off()
### cutting edge作图
pdf("cnerplot.pdf",20,20)
cnetplot(y,showCategory = 10,foldChange = geneList,colorEdge = T)
dev.off()
####################################################################
### 筛选特定的通路来画图
### install.packages("clusterProfiler.dplyr-master/",repos=NULL,type="source")
### https://yulab-smu.github.io/clusterProfiler-book/chapter13.html
library(clusterProfiler)
### 
yd <- data.frame(y)
### 提取NES大于0的，也就是激活的
y1 <- filter(y,NES > 0)
dotplot(y1,showCategory=30,split=".sign")+facet_grid(~.sign)
y2 <- filter(y,NES < 0)
dotplot(y2,showCategory=30,split=".sign")+facet_grid(~.sign)
### 自定义画图
pdf("gseadotplot_direction.pdf",25,12)
ggplot(y, showCategory = 30, aes(NES, forcats::fct_reorder(Description, NES))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  ##scale_color_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL)
dev.off()
#################################################################
### 2.kegg 通路
## 读入kegg gene set
hallmarks <- read.gmt("mh.all.v2023.1.Mm.symbols.gmt")
kegg <- GSEA(geneList,TERM2GENE =hallmarks,nPerm = 10000)
keggd <- data.frame(kegg)
### 看整体分布
pdf("gseadotplot_KEGG.pdf",25,12)
dotplot(kegg,showCategory=12,split=".sign")+facet_grid(~.sign)
dev.off()
#################################################################
### 3.转录因子
## 读入转录因子
hallmarks <- read.gmt("ENCODE_TF_ChIP-seq_2015.txt")
tfbs <- GSEA(geneList,TERM2GENE =hallmarks,nPerm = 10000)
tfbsd <- data.frame(tfbs)
### 看整体分布
dotplot(tfbs,showCategory=12,split=".sign")+facet_grid(~.sign)


