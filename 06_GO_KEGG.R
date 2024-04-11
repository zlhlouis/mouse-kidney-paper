
#########
#########
# GO分析，KEGG分析
###########################
###########################
rm(list = ls())
setwd("E://SomeProje/Significant expression genes  Project/GSE131240//")
#setwd("F://差异表达基因/8-21/GSE9629/")
library(clusterProfiler)
load(file = "allDiff.Rda")
diffLab <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 1e-2) %>% 
  filter(abs(logFC) > 2) %>% 
  column_to_rownames()
#获得基因列表
gene1 <- rownames(subset(diffLab,diffLab$logFC>0))
gene2 <- rownames(subset(diffLab,diffLab$logFC<0))
#基因名称转换，返回的是数据框
gene1 = bitr(gene1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
gene2 = bitr(gene2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
#save(gene,file = "gene.Rdata")
#load(file = "file = "gene.Rdata")
head(gene)

#**GO分析**
#细胞组分CC
if(F)
  {
  ego_CC <- enrichGO(gene = gene$ENTREZID,
                     OrgDb= "org.Mm.eg.db",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
}
save(ego_CC,file = "ego_CC.Rdata")
#write.table(ego_CC,file="./part1/ego_CC.pro",sep='\t',quote=F)
library(tidyr)
load(file = "ego_CC.Rdata")
eGo <- separate(data=ego_CC, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=ego_CC, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(ego_CC, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor 


#**作图**
#条带图
pdf("131240_CC_barplo.pdf",10,12)
barplot(ego_CC,showCategory=30)
dev.off()
# 点图
pdf("part1_overlap_CC_dotplo.pdf",10,6)
dotplot(ego_CC)
dev.off()
# GO 作图
pdf("GSE9629_CC_goplo.pdf",10,6)
goplot(ego_CC)
dev.off()

#GO分析的生物过程BP，需要网络
if(F){
  ego_BP <- enrichGO(gene = gene$ENTREZID,
                     OrgDb= org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
}
write.table(ego_BP,"ego_BP.pro",sep='\t',quote=F)
#load(file = "ego_BP_20190316.Rdata")
#**作图**
#条带图
library(RColorBrewer)

pdf("BP_barplot.pdf",10,6)
barplot(ego_BP)
dev.off()
# 点图
pdf(".BP_dotplot.pdf",10,6)
dotplot(ego_BP)
dev.off()
# GO 作图
pdf("./part1/BP_goplot.pdf",20,12)
goplot(ego_BP)
dev.off()
###########################
pdf("./part1/BP_plot.pdf",20,15)
cnetplot(go1, categorySize="pvalue", foldChange=gene1)
dev.off()

#GO分析分子功能MF：
if(F){
  ego_MF <- enrichGO(gene = gene$ENTREZID,
                     OrgDb= org.Mm.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
}
write.table(ego_MF,file = "./part1/ego_MF.pro",sep='\t',quote=F)

load(file = "ego_MF.Rdata")
#**作图**
#条带图
pdf("./part1/MF_barplot.pdf",10,6)
barplot(ego_MF)
dev.off()
# 点图
pdf("./part1/MF_dotplot.pdf",10,6)
dotplot(ego_MF)
dev.off()
# GO 作图
pdf("./part1/MF_goplot.pdf",20,6)
ego_MF1=na.omit(ego_MF)
goplot(ego_MF1)
dev.off()
###另外一种保存图片的方式
if(F){
  pdf(file = "dotplot_20190428.pdf")# 尝试随便命名
  dotplot(ego_BP,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
  dev.off()
}
# 推荐保存方法，plots，export，device

############################################
############################################
############################### zhaohui 

eG <- enrichGO(gene = gene$ENTREZID, #需要分析的基因的EntrezID
               OrgDb = "org.Mm.eg.db", #人基因数据库
               pvalueCutoff =0.01, #设置pvalue界值
               qvalueCutoff = 0.01, #设置qvalue界值(FDR校正后的p值）
               ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
               readable =T)
write.table(eG,file="eG.txt", sep="\t", quote=F, row.names = F)
######################  barplot
pdf(file="eGO_barplot.pdf",width = 8,height = 10) 
barplot(eG, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
        showCategory =10, #只显示前10
        split="ONTOLOGY") + #以ONTOLOGY类型分开
  facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
dev.off()

################  dotplot
dotplot(eG,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =5,#只显示前5
        split="ONTOLOGY") + #以ONTOLOGY类型分开
  facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分屏绘图

###########################
eGo <- read.table("eG.txt",header=TRUE,sep="\t",quote = "") 
eGo <- separate(data=eGo, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
 #计算Enrichment Factor 
eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)

pdf("0626eGO_barplot1.pdf", 11, 10) 
ggplot(eGo10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()+
  facet_wrap( ~ ONTOLOGY)
dev.off()

#######################  GO bubble
library(GOplot)
GOterms = data.frame(category = eGo10$ONTOLOGY,
                     ID = eGo10$ID,
                     term = eGo10$Description, 
                     genes = gsub("/", ",", eGo10$geneID), 
                     adj_pval = eGo10$p.adjust)
diffLab1 <- diffLab[pmatch(gene$SYMBOL,rownames(diffLab)),]
genelist <- data.frame(ID = gene$SYMBOL, logFC = diffLab1$logFC)
circ <- circle_dat(GOterms, genelist)
pdf("0626eGO_bubble_plot2.pdf", 12, 10)
GOBubble(circ, labels = 3, # 标注的界值：-log(adjusted p-value) (默认5)
         table.legend =T, #不显示右侧legend
         ID=T, # 标注term名称
         display='single')

reduced_circ <- reduce_overlap(circ,overlap=0.75) 
GOBubble(reduced_circ, labels = 2.8,display = 'single', bg.col = T)

GOBubble(circ, title = 'Bubble plot with background colour',
         display = 'multiple', bg.col = T, labels = 3)
dev.off()

GOCircle(circ)
pdf("eGO_circle_plot.pdf", 16, 14)
GOCircle(circ,rad1=20, #内环半径
         rad2=30, #外环半径
         label.size= 5, #标签大小
         label.fontface = 'bold', #标签字体
         nsub=20, #显示的terms数，前10个。（画图前需要先把数据整理好，想要哪些term）
         zsc.col = c('red', 'steelblue'), # z-score颜色
         lfc.col = c('red', 'steelblue')) # 基因up或down的颜色
dev.off()

###########################################  guozi
if(F){
  go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all")
}
##save(go,file="go.Rdata")
load(file = "go.Rdata")
library(ggplot2)
pdf("BP_varplotGO.pdf",10,6)
barplot(ego_BP, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
dev.off()
pdf("BP_GO_cnetplot.pdf",25,18)
cnetplot(ego_BP,categorySize="pvalue",circular=T,colorEdge=T，)
dev.off()

go<-as.data.frame(ego_BP)
GO<-ego_BP[1:9,c(1,2,8,6)]
GO$geneID<-str_replace_all(GO$geneID,"/",",")
names(GO)<-c("ID","Term","Genes","adj_pval")
GO$Category="BP"
genedata<-as.data.frame(cbind(rownames(dat),dat[,1]))
colnames(genedata)<-c("ID","logFC")
circ<-circle_dat(GO,genedata)
############hexiantu



############tiaoxingtu
pdf("BP_GObar.pdf",10,6)
GOBar(circ,display="multiple")
dev.off()
###################
pdf("BP_GOBubble.pdf",20,10)
GOBubble(circ,title="Bubble plot",colour = c("orange","darkred","gold"),display="multiple",label=3)
dev.off()
pdf("BP_GOCluster.pdf",10,6)
GOCluster(circ,EC$process,clust.by = "logFC")
dev.off()


####################################################
#**KEGG分析**

library(KEGG.db)
EGG <- enrichKEGG(gene = gene2$ENTREZID,
                 organism = 'mmu',
                 pvalueCutoff = 0.05,
                 use_internal_data =F)
write.table(EGG,file = "down_EGG.pro",sep='\t',quote=F)
## 画图
test <- as.data.frame(EGG)
pdf("KEGG_barplot.pdf",10,6)
barplot(EGG)
dev.off()
pdf("down_KEGG_dotplot.pdf",10,6)
dotplot(EGG)
dev.off()
## rich factor
## 定制一个图片
if(T){
  x = EGG
  ## 内置的函数可以转换为数据框
  df = data.frame(x)
  dd =x@result
  ## 计算富集分数
  dd$richFactor =dd$Count / as.numeric(sub("/\\d+", "", dd$BgRatio))
  ## 提取p值小于0.05 的
  dd <- dd[dd$p.adjust < 0.05,]
  
  library(ggplot2)
  ## 正式画图
  pdf("KEGG_plot",10,6)
  ggplot(dd,aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
    ## 画横线
    geom_segment(aes(xend=0, yend = Description)) +
    ## 画点
    geom_point(aes(color=p.adjust, size = Count)) +
    ## 调整颜色的区间,begin越大，整体颜色越明艳
    scale_color_viridis_c(begin = 0.3, end = 1) +
    ## 调整泡泡的大小
    scale_size_continuous(range=c(2, 10)) +
    theme_bw() + 
    xlab("Rich factor") +
    ylab(NULL) + 
    ggtitle("")
}
dev.off()
### 使用clusterProfiler.dplyr的mutate功能
library(clusterProfiler.dplyr)
### richFactor
y=mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
test <- as.data.frame(y)
pdf("KEGG.pdf",20,10)
ggplot(y,showCategory = 30,aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
  ## 画横线
  geom_segment(aes(xend=0, yend = Description)) +
  ## 画点
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 调整颜色的区间,begin越大，整体颜色越明艳
  scale_color_viridis_c(begin = 0.3, end = 1) +
  ## 调整泡泡的大小
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Rich factor") +
  ylab(NULL) + 
  ggtitle("")

dev.off()
###################################################################################################################
## 图谱展示
KEGG_df <- data.frame(EGG)
browseKEGG(EGG, 'hsa04110')


