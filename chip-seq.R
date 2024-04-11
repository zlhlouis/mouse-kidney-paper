
rm(list=ls())
setwd("E://SomeProje/Significant expression genes  Project/0706Manuscript")
library(readxl)
data1 <- read_excel("CHiP_data.xlsx",sheet = 1)
data2 <- read_excel("CHiP_data.xlsx",sheet = 2)
id <- intersect(data1$GeneSymbol,data2$GeneSymbol)
id <- as.data.frame(id)
write.table(data21,"chip_new_overlap.txt",sep='\t',quote=F)
data11 <- data1[pmatch(id$id,data1$GeneSymbol),]
data21 <- data2[pmatch(id$id,data2$GeneSymbol),]
colnames(data11) <- c("GSM980186_Chromosome","GSM980186_Start","GSM980186_End","GeneSymbol")
colnames(data21) <- c("GSM980187_Chromosome","GSM980187_Start","GSM980187_End","GeneSymbol")
data_new <- cbind(data11,data21)

deg1 <- read.table("../0424_GSE39583_20325/diffLab_part1_chayi.txt",sep='\t',head=T,row.names = 1)
deg2 <- read.table("../0424_GSE39583_20325/diffLab_part2_chayi.txt",sep='\t',head=T,row.names = 1)
deg3 <- read.table("../0424_GSE39583_20325/diffLab_part3_chayi.txt",sep='\t',head=T,row.names = 1)

deg1_chip <- intersect(data_new$GeneSymbol,rownames(deg1))
deg2_chip <- intersect(data_new$GeneSymbol,rownames(deg2))
deg3_chip <- intersect(data_new$GeneSymbol,rownames(deg3))
deg1_chip1 <- data_new[pmatch(deg1_chip,data_new$GeneSymbol),]
deg2_chip1 <- data_new[pmatch(deg2_chip,data_new$GeneSymbol),]
deg3_chip1 <- data_new[pmatch(deg3_chip,data_new$GeneSymbol),]

write.table(deg1_chip1,"chip1_overlap.csv",sep=',',quote=F)
write.table(deg2_chip1,"chip2_overlap.csv",sep=',',quote=F)
write.table(deg3_chip1,"chip3_overlap.csv",sep=',',quote=F)

#####################################  DEG1 plot
deg3_chip1<- deg3_chip1[,-8]
data<- deg3_chip1[,c(4,1)]
data[,2] <- str_replace_all(data[,2],"chr","")
data$GSM980186_Chromosome <- as.numeric(data$GSM980186_Chromosome)
# Set a number of 'empty bar'
empty_bar <- 0
# Add lines to the initial dataset
to_add <- matrix(NA, empty_bar, ncol(data))
colnames(to_add) <- colnames(data)
data <- rbind(data, to_add)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Make the plot
pdf("deg3_circle.pdf",8,6)
ggplot(data, aes(x=as.factor(id), y=GSM980186_Chromosome)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", fill=alpha("#20B2AA", 0.3)) +
  ylim(-1.5,20) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar(start = 0) + 
  geom_text(data=label_data, aes(x=id, y=GSM980186_Chromosome+0.2, label=GeneSymbol, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) 
dev.off()

####################### venn plot
venn.diagram(x=list(data1$GeneSymbol,data2$GeneSymbol),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#58A4C3','#E4C755'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#58A4C3','#E4C755'), # 填充色 配色https://www.58pic.com/
             
             category.names = c("GSM980186", "GSM980187") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 1, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='chip_veen.pdf'
             
)


neen_data <- cbind(data1$GeneSymbol,data2$GeneSymbol)
venn_list <- as.list(data1$GeneSymbol,data2$GeneSymbol)              # 制作韦恩图搜所需要的列表文件
venn_list <- purrr::map(venn_list, na.omit) # 删除列表中每个向量中的NA

# 绘图
library(ggvenn)
ggvenn(
  data = venn_list,         # 数据列表
  columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c('#58A4C3','#E4C755'), # 填充颜色
  fill_alpha = 0.5,         # 填充透明度
  stroke_color = "white",   # 边缘颜色
  stroke_alpha = 0.5,       # 边缘透明度
  stroke_size = 0.5,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 4             # 交集个数文字大小
)

