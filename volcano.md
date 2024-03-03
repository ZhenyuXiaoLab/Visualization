# 火山图

## 导入DEG数据

```R
#导出DEG
#cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)  #find cluster0 的差异基因
#write.table(cluster0.markers, './cluster0.markers.txt',sep='\t')    #导出 cluster0 的差异基因

#导入DEG
c0 <- read.table('./cluster0.markers.txt')
c1 <- read.table('./cluster1.markers.txt')
c <- rbind(c0, c1)
```



```R
# 自定义log2FC和p-val阈值
FC = 1.5
log2(FC)
P.Value = 0.05

# 根据阈值给DEG添加change列存放差异结果信息
k1 <- (c$p_val < P.Value) & (c$avg_log2FC < -log2(FC))
k2 <- (c$p_val < P.Value) & (c$avg_log2FC > log2(FC))
c <- mutate(c, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
head(c)
```



```R
# 添加一列GeneName
c$GeneName <- rownames(c)
# 筛选显著性top10标签：
top10 <- filter(c, change != "none") %>%
  distinct(GeneName, .keep_all = T) %>%
  top_n(10, abs(avg_log2FC))

# top100
top100 <- filter(c, change != "none") %>%
  distinct(GeneName, .keep_all = T) %>%
  top_n(100, abs(avg_log2FC))
```



```R
# 加载R包
library(dplyr)
library(ggplot2)
library(ggVolcano)
library(RColorBrewer)
library(ggrepel)
```



## volcano1

```R
huoshan <- ggplot(data = c, 
            aes(x = avg_log2FC, 
                y = -log10(p_val))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-log2(FC), log2(FC)), lty = 4, col = "black", lwd = 0.8) +  # log2FC阈值线
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +         # p-val阈值线
  theme_bw()
ggsave("./figure/huoshan.png",plot=huoshan,device="png",width=8)


# 添加top10基因标签
huoshanlabel <- huoshan +
  geom_text_repel(data = top10,
                  aes(x = avg_log2FC, y = -log10(p_val), label = GeneName),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 0.5, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.3,
                  nudge_x = 0,
                  nudge_y = 1)
ggsave("./figure/huoshanlabel.png",plot=huoshanlabel,device="png",width=8)
```

![1](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/1.png)



## volcano2

```R
# 换种主题
huoshan <- ggplot(data=c, 
            aes(x = avg_log2FC, y = -log10(p_val),color=change)) + 
  geom_point(alpha=0.4, size=1.75) + 
  theme_set(theme_set(theme_bw(base_size=20))) + 
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  ggtitle("this_tile") + 
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('blue','black','red'))
ggsave("./figure/huoshan.png",plot=huoshan,device="png",width=8)


# 添加top10基因标签
huoshanlabel <- huoshan + geom_text_repel(data = top10, aes(x = avg_log2FC, 
                                    y = -log10(p_val), 
                                    label = GeneName),
                    size = 3,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE)
ggsave("./figure/huoshanlabel.png",plot=huoshanlabel,device="png",width=8)


# 修改label样式
huoshanlabel2 <- huoshan +
  geom_point(size = 3, shape = 1, data = top10) +
  ggrepel::geom_label_repel(
    aes(label = GeneName),
    data = top10,
    color="black"
  )
ggsave("./figure/huoshanlabel2.png",plot=huoshanlabel2,device="png",width=8)
```

![微信截图_20231224182716](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/微信截图_20231224182716.png)



## volcano3

```R
#根据反比例函数 y = 1/X和自己设定的阈值自定义双曲线函数
f <- function(x){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(P.Value))
  dff <- rbind(data.frame(x = inputx + log2(FC), y = y),
               data.frame(x = -(inputx + log2(FC)), y = y))
  return(dff)
}


#根据函数生成所需的曲线数组坐标：
dff_curve <- f(4)
head(dff_curve)
#需要新增曲线数值列：
##每列log2FoldChange值在曲线上对应的y轴坐标；
c$curve_y <- case_when(
  c$avg_log2FC > 0 ~ 1/(c$avg_log2FC-log2(FC)) + (-log10(P.Value)),
  c$avg_log2FC <= 0 ~ 1/(-c$avg_log2FC-log2(FC)) + (-log10(P.Value))
)

#根据双曲线的阈值范围新增上下调分组标签，存放在change2列：
c$change2 <- case_when(
  -log10(c$p_val) > c$curve_y & c$avg_log2FC >= log2(FC) ~ 'up',
  -log10(c$p_val) > c$curve_y & c$avg_log2FC <= -log2(FC) ~ 'down',
  TRUE ~ 'none'
)

#根据双曲线筛选显著性top10标签：
c$GeneName <- rownames(c)
top10_shuangquxian <- filter(c, change2 != "none") %>%
  distinct(GeneName, .keep_all = T) %>%
  top_n(10, abs(avg_log2FC))


mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# 制作双曲线火山图
huoshan <- ggplot(data = c,
             aes(x = avg_log2FC, y = -log10(p_val))) +  #使用新的分组
  geom_point(size = 2.2, aes(color = change2)) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2.5)) +  #限制x坐标轴范围
  scale_y_continuous(expand = expansion(add = c(2, 0)),limits = c(0, 300), breaks = seq(0, 300, by = 100)) +  #限制y坐标轴范围
  scale_colour_manual(values = c("#4A1985","#d8d8d8","#F8B606")) +  # 修改颜色
  geom_line(data = dff_curve,
            aes(x = x, y = y), #曲线坐标
            color = "black",lty = "dashed", size = 0.7) +  # 阈值线使用做好的双曲线
  mytheme
  #geom_vline(xintercept = c(-log2(FC), log2(FC)), lty = 4, col = "black", lwd = 0.8) +  # log2FC阈值线
  #geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8)           # p-val阈值线

ggsave("./figure/huoshan.png",plot=huoshan,device="png",width=8)
```

![image-20231224183621223](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/image-20231224183621223.png)



## volcano4

```R
# ggVolcano包作图
# ggVolcano包会将原数据列名改掉
# 不需要其他处理，导入数据直接进行
library(ggVolcano)

# 加一列GeneName
c$GeneName <- rownames(c)
# 加一列regulate存放上下调信息
c <- add_regulate(c, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val",log2FC = 1, fdr = 0.05)

# 作火山图
huoshan <- ggvolcano(c, x = "log2FoldChange", y = "padj",
                    fills = c("#e94234","#b4b4d8","#269846"), #修改散点颜色
                    colors = c("#e94234","#b4b4d8","#269846"), #修改散点颜色
          label = "GeneName", label_number = 10, output = FALSE)
ggsave("./figure/huoshan4.png",plot=huoshan,device="png",width=8)


# 作渐变火山图
gradual_huoshan <- gradual_volcano(c, x = "log2FoldChange", y = "padj",
                                   fills = brewer.pal(5, "RdYlBu"), #利用RColorBrewer包修改渐变颜色
                                   colors = brewer.pal(8, "RdYlBu"), #利用RColorBrewer包修改渐变颜色
                label = "GeneName", label_number = 10, output = FALSE)
ggsave("./figure/gradual_huoshan.png",plot=gradual_huoshan,device="png",width=8)
```

![微信截图_20231224185007](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/微信截图_20231224185007.png)



## volcano5

```R
# EnhancedVolcano包作图
library(EnhancedVolcano)

huoshan <- EnhancedVolcano(
  c,
  x = "avg_log2FC",
  y = "p_val",
  lab = c$GeneName,    # 基因名
  selectLab = rownames(top10),
  pCutoff = 0.001,      # p值阈值：水平线
  FCcutoff = 0.8,       # FC阈值：垂直线
  cutoffLineWidth = 0.7,      # 虚线粗细
  cutoffLineType = "twodash", # 线的类型
  xlim = c(-5.0, 5.0),  # x轴起始
  ylim = c(0, 300),     # y轴起始
  pointSize = 3,        # 点的大小
  labSize = 2.5,        # 标签大小
  xlab = bquote(~Log[2] ~ "fold change"),     # 此行为默认x轴名
  ylab = bquote(~-Log[10] ~ italic(p-value)), # 修改y轴名为斜体
  axisLabSize = 12,     # 坐标轴字体大小
  title = "Gene differential expression",     # 修改标题
  titleLabSize = 16,    # 标题大小
  subtitle = bquote(italic("Volcano plot")),  # 修改子标题：斜体
  subtitleLabSize = 14, # 子标题大小
  legendLabSize = 11,   # legend大小
  col = c("#a5a5a5", "#8fc490","#9eb5eb","#f87f7f"), # legend颜色
  colAlpha = 0.4,       # 透明度
  gridlines.major = TRUE,     # 背景网格
  gridlines.minor = TRUE,
  #drawConnectors = TRUE,  # 标记基因是否显示标记点
  #widthConnectors = 1.0   # 标记点参数
  #colConnectors = 'black' # 标记点颜色
)
ggsave("./figure/huoshan.png",plot=huoshan,device="png",width=8)
```

![image-20231224203846944](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/image-20231224203846944.png)



## volcano6

```R
# 可以导入多组DEG
c0$group <- 'cluster0'
c1$group <- 'cluster1'
c2$group <- 'cluster2'
c3$group <- 'cluster3'
c4$group <- 'cluster4'
c5$group <- 'cluster5'
c6$group <- 'cluster6'
c7$group <- 'cluster7'
c <-rbind(c0, c1,c2,c3,c4,c5,c6,c7)
head(c)

# 自定义log2FC和p-val阈值
FC = 1.5
log2(FC)
P.Value = 0.05

# 根据阈值给DEG添加change列存放差异结果信息
k1 <- (c$p_val < P.Value) & (c$avg_log2FC < -log2(FC))
k2 <- (c$p_val < P.Value) & (c$avg_log2FC > log2(FC))
c <- mutate(c, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
head(c)

# 添加一列GeneName
c$GeneName <- rownames(c)
# 筛选显著性top10标签：
top10 <- filter(c, change != "none") %>%
  distinct(GeneName, .keep_all = T) %>%
  top_n(10, abs(avg_log2FC))

# 如果数据中没有size信息，引入top100基因作散点个大小差异渐变
# top100
top100 <- filter(c, change != "none") %>%
  distinct(GeneName, .keep_all = T) %>%
  top_n(100, abs(avg_log2FC))

####################### 以上涉及可能重复部分，自定义修改
# 开始作图
# 提出x轴和y轴数据
fc <-c$avg_log2FC
names(fc) <-rownames(c)

p <- -log10(c$p_val)
names(p) <-rownames(c)

# 给每个cluster的泡泡一种颜色
# 先自定义足够多的颜色
mycol <- c("#B2DF8A","#FB9A99","#33A02C","#E75B5B","#B15928","#8858BC","#CAB2D6","#A6CEE3","#1F78B4","#FDBF6F","#999999","#FF7F00")
cols.names <- unique(c$group)
cols.code <- mycol[1:length(cols.names)]
names(cols.code) <- cols.names
col <- paste(cols.code[as.character(c$group)],"BB", sep="")
# Highlight 感兴趣的部分，这里选择cluster0/1/2三个类，以及top10/top100基因
i <- c$group %in% c("cluster0","cluster1") # 可以设置感兴趣的pathway
ii <- c$group %in% c("cluster2")
top <- c$GeneName %in% top10$GeneName
topp <- c$GeneName %in% top100$GeneName

### size列
# sizes <- c$size
# names(sizes) <- rownames(c )

### pval列
# pp <- c$p_val
# names(pp) <- rownames(c)

pdf('huoshan.pdf')
par(xpd = FALSE,
    mar = par()$mar + c(0,0,0,6)) # 在右侧留出画图例的地方
plot(fc, p, #log = 'y',           # x轴，y轴，对y轴进行log变换
     col = paste(cols.code[as.character(c$group)], "BB", sep = ""),
     pch = 16,
     ylab = bquote(~Log[10]~"P value"), xlab = "log2FC",
     cex = ifelse(top,2,ifelse(topp,1,0.5)), # 大泡泡画top10基因，中泡泡画top10-100基因，小泡泡画不感兴趣的基因，对通路pathway同理
     xlim = range(fc * 1.2))

# 添加横线
abline(h=-log10(0.05), lty=2, lwd=1)
# 添加竖线
abline(v=-log2(1.5), col="blue", lty=2, lwd=1)
abline(v=log2(1.5), col="red", lty=2, lwd=1)

# 添加size的图例
par(xpd = TRUE) #all plotting is clipped to the figure region
f <- c(0.5,1,2) # 泡泡大小数值，于前文一致
s <- sqrt(f*3)  # 设置图例泡泡大小，按等比例修改
legend("topright",
       inset=c(-0.2,0), #把图例画到图外
       legend=f, pch=16, pt.cex=s, bty='n', col=paste("#88888888"))
# 添加cluster颜色的图例
legend("bottomright", 
       inset=c(-0.25,0), #把图例画到图外
       pch=16, col=cols.code, legend=cols.names, bty="n")
dev.off()
```

![image-20231224212547806](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/image-20231224212547806.png)



## volcano7

```R
#自定义颜色与主题：
mycol <- c("#ffa500","#d8d8d8","#B3A9EB")

mytheme <- theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        plot.margin = margin(15,5.5,5.5,5.5)) #图边缘空白调整，顺时针方向(上右下左)

#常规火山图绘制：
p <- ggplot(data = c,
            aes(x = avg_log2FC, y = -log10(p_val))) + #建立映射
  geom_point(alpha = 0.4, size = 2.6, 
             aes(color = change)) + #添加散点
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + #自定义散点颜色
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2)) + #x轴范围限制
  scale_y_continuous(expand = expansion(add = c(2, 0)), limits = c(0, 200), breaks = seq(0, 200, by = 20)) + #y轴范围限制
  geom_hline(yintercept = c(-log10(P.Value)),size = 0.7,color = "black",lty = "dashed") + #水平阈值线
  geom_vline(xintercept = c(-log2(FC), log2(FC)),size = 0.7,color = "black",lty = "dashed") + #垂直阈值线
  mytheme

#自定义主题(将原本y轴隐藏掉)：
mytheme2 <- theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(15,5.5,5.5,5.5),
        axis.line.y = element_blank(), #y轴线隐藏
        axis.text.y = element_blank(), #y轴文本标签隐藏
        axis.ticks.y = element_blank()) #y轴刻度线隐藏

p1 <- ggplot() +
  geom_vline(xintercept = 0, size = 0.5, color = "black") + #先在原点（0,0）添加一条垂直线，作为y轴的替代
  geom_point(data = c,
             aes(x = avg_log2FC, y = -log10(p_val)),
             size = 2.6, color = 'grey', alpha = 0.6) +
  scale_x_continuous(limits = c(-5, 5),
                     breaks = seq(-5, 5, by = 2)) +
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 210),
                     breaks = seq(0, 210, by = 20))
p2 <- p1 + mytheme2

# y轴刻度值添加：
## 要根据前文y轴范围限制进行对应创建
### 先创建刻度值数据集；
dty <- data.frame(x = rep(0,10),
                  y = seq(20,200,20),
                  label = rep('—',10)) #"伪造"刻度线
dty


p3 <- p2 +
  geom_text(data = dty,
            aes(x = x, y = y, label = label), fontface = 'bold') +
  geom_text(data = dty,
            aes(x = x, y = y, label = y), hjust = 1.5)

#筛选显著性top10标签：
##c$GeneName <- rownames(c)
#top10 <- filter(c, change != "none") %>%
#  distinct(GeneName, .keep_all = T) %>%
#  top_n(10, abs(avg_log2FC))

#筛选显著性top20标签：
top20 <- filter(c, change != "none") %>%
  distinct(GeneName, .keep_all = T) %>%
  top_n(20, abs(avg_log2FC))

#筛选显著性top10-20标签，此为手动比对筛选
#top1020 <-top20[c(5,6,7,9,10,11,12,15,16,20),]

#top10标注：
p4 <- p3 +
  geom_point(data = top20,
             aes(x = avg_log2FC, y = -log10(p_val)),
             color = '#FE88B1', size = 3.5, alpha = 1) +
  geom_text_repel(data = top20,
                  aes(x = avg_log2FC, y = -log10(p_val), label = GeneName),
                  size = 3, fontface = 'bold.italic')
#top1020标注：
p5 <- p4 +
  geom_point(data = top1020,
             aes(x = avg_log2FC, y = -log10(p_val)),
             color = '#9EB9F3', size = 3.5, alpha = 1) +
  geom_text_repel(data = top1020,
                  aes(x = avg_log2FC, y = -log10(p_val), label = GeneName),
                  size = 3, fontface = 'bold.italic')
huoshan <-p5
ggsave("./figure/huoshan.png",plot=huoshan,device="png",width=8)
```

![image-20231224213557244](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/image-20231224213557244.png)

# multi-volcano

输入文件格式

![微信截图_20231224213931](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/微信截图_20231224213931.png)

```R
# 示例数据
url <- "https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/Volcano/mutiVolcano.txt"
  
destination<- "./mutiVolcano.txt"

# 下载示例数据  
download.file(url, destination)


#加载包
library(tidyverse)
library(ggrepel)   # 用于标记

# 自定义一个函数，之后我们直接调用它
mutiVolcano = function(df,         # 绘图数据
                       P = 0.05,   # P值卡值
                       FC = 1.5,   # FC卡值
                       GroupName = c("Sig","Not Sig"),      # 分组标签
                       pointColor = c("#CC3333","#0099CC"), # 分组散点的颜色
                       barFill = "#efefef",  # 柱子的颜色
                       pointSize = 0.9,      # 散点的大小
                       labeltype = "1",      # 标记差异基因的选项，标记类型有"1"和"2"两种选项
                       labelNum = 5,         # 当标记类型为1时，待标记的散点个数
                       labelName =NULL,      # 当标记类型为2时，待标记的散点名称
                       tileLabel =  "Label", # 标记比较对的选项，选项有“Label”和“Num”，Label时显示分组名称，Num时显示数字，防止因为标签太长导致的不美观
                       tileColor = NULL      # 比较对的颜色
                       ){
  # 数据分组 根据p的卡值分组
  dfSig = df %>% 
    mutate(log2FC = log2(FC)) %>%
    filter(FC>{{FC}} | FC <(1/{{FC}})) %>%
    mutate(Group = ifelse(PValue<0.05,GroupName[[1]],GroupName[[2]])) %>%
    mutate(Group = factor(Group,levels=GroupName)) %>%
    mutate(Cluster = factor(Cluster,levels=unique(Cluster)))   # Cluster的顺序是文件中出现的顺序
  
  # 柱形图数据整理
  dfBar = dfSig %>%
    group_by(Cluster) %>%
    summarise(min = min(log2FC,na.rm = T),
              max = max(log2FC,na.rm = T)
              )
  # 散点图数据整理
  dfJitter = dfSig %>%
    mutate(jitter = jitter(as.numeric(Cluster),factor = 2))
  dfJitter
  
  # 整理标记差异基因的数据
  if(labeltype == "1"){
    # 标记一
    # 每组P值最小的几个
    dfLabel = dfJitter %>%
      group_by(Cluster) %>%
      slice_min(PValue,n=labelNum,with_ties = F) %>%
      ungroup()
  }else if(labeltype == "2"){
    # 标记二
    # 指定标记
    dfLabel = dfJitter %>%
      filter(Name %in% labelName)
  }else{
    dfLabel = dfJitter %>% slice()
  }
    
  # 绘图
  p = ggplot()+
    # 绘制柱形图
    geom_col(data = dfBar,aes(x=Cluster,y=max),fill = barFill)+
    geom_col(data = dfBar,aes(x=Cluster,y=min),fill = barFill)+
    # 绘制散点图
    geom_point(data = dfJitter,
               aes(x = jitter, y = log2FC, color = Group),
               size = pointSize,
               show.legend = NA
               )+
    # 绘制中间的标签方块
    ggplot2::geom_tile(data = dfSig,
                       ggplot2::aes(x = Cluster, y = 0, fill = Cluster), 
                       color = "black",
                       height = log2(FC) * 1.5,
                       # alpha = 0.3,
                       show.legend = NA
                       ) + 
    # 标记差异基因
    ggrepel::geom_text_repel(
      data = dfLabel,
      aes(x = jitter,                   # geom_text_repel 标记函数
          y = log2FC,          
          label=Name),        
      min.segment.length = 0.1,
      max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
      size=3,                                  # 字体大小
      box.padding=unit(0.5,'lines'),           # 标记的边距
      point.padding=unit(0.1, 'lines'), 
      segment.color='black',                   # 标记线条的颜色
      show.legend=F)#+
  
  if(tileLabel=="Label"){
     p =
      p +
      geom_text(data = dfSig,aes(x = Cluster,y = 0,label = Cluster))+
      ggplot2::scale_fill_manual(values = tileColor,
                                 guide = NULL # 不显示该图例
                                 )
  }else if(tileLabel=="Num"){
    # 如果比较对的名字太长，可以改成数字标签
    p =
      p +
      geom_text(data = dfSig,aes(x = Cluster,y = 0,label = as.numeric(Cluster)),show.legend = NA)+
      ggplot2::scale_fill_manual(values = tileColor,
                                 labels = c(paste0(1:length(unique(dfSig$Cluster)),": ",unique(dfSig$Cluster))))
  }

  
    
  # 修改主题
  p = p+ggplot2::scale_color_manual(values = pointColor)+
    theme_classic()+
    ggplot2::scale_y_continuous(n.breaks = 5) + 
    ggplot2::theme(
                   legend.position = "right", 
                   legend.title = ggplot2::element_blank(), 
                   legend.background = ggplot2::element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   # axis.title.x = element_blank(),
                   axis.line.x = element_blank()
                   ) + 
    ggplot2::xlab("Clusters") + ggplot2::ylab("log2FC") + 
    # ggplot2::guides(fill = ggplot2::guide_legend())
    guides(color = guide_legend(override.aes = list(size = 3)))
    
    return(p)
}

df = read.delim("./mutiVolcano.txt") %>%  
  as_tibble() %>%
  set_names(c("Name","FC","PValue","Cluster"))

# 调用函数画图
mutiVolcano(
  df = df,    # 绘图数据
  P = 0.05,   # P值卡值
  FC = 1.5,   # FC卡值
  GroupName = c("Sig","Not Sig"),      # 分组标签
  pointColor = c("#CC3333","#0099CC"), # 分组散点的颜色
  barFill = "#efefef",   # 柱子的颜色
  pointSize = 0.9,       # 散点的大小
  labeltype = "1",       # 标记差异基因的选项，标记类型有"1"和"2"两种选项
  labelNum = 5,          # 当标记类型为1时，待标记的散点个数
  labelName =c("ID1","ID2029"),   # 当标记类型为2时，待标记的散点名称
  tileLabel =  "Label",           # 标记比较对的选项，选项有“Label”和“Num”，Label时显示分组名称，Num时显示数字，防止因为标签太长导致的不美观
  tileColor = RColorBrewer::brewer.pal(length(unique(df$Cluster)),"Set3")   # 比较对的颜色
)

# 保存图片
ggsave("./figure/multiVolcano.png",width = 8,height = 6,dpi=600)
```

![image-20231224214040933](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/image-20231224214040933.png)

# G0-volcano

go输入文件格式（列名也要相同）

![微信截图_20231224220031](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/微信截图_20231224220031.png)

测试结果得出，两个文件（DEG文件和GO文件）基因数目可以不匹配

```R
# 使用ggVolcano包
# ggVolcano包会将原数据列名改掉
library(readxl)
library(ggVolcano)
library(ggplot2)
library(RColorBrewer)

# 读取DEG数据
c1 <- read.table('./cluster1.markers.txt')
c2 <- read.table('./cluster2.markers.txt')
c <- rbind(c1, c2)
# 添加一列上下调信息
c <- add_regulate(c, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val",log2FC = 1, fdr = 0.05)
# 添加一列GeneName
c$GeneName <- rownames(c)


# 读取go信息
g1 <- read_excel('./1_GO_enrichment.txt.xlsx')
g2 <- read_excel('./2_GO_enrichment.txt.xlsx')

# 根据自己的go信息文件制作符合的文件格式，最后存放在变量ggg中
g1 <- g1[,c(3,9)]
g2 <- g2[,c(3,9)]
g <- rbind(g1, g2)
g[,3] <- g[,1]
g <- g[,-1]

gg <- data.frame()
genes <- rownames(c)

for (gene in genes){
    for (i in 1:nrow(g)){
        if(grepl(gene, g[i,1])){
            gg <- rbind(gg, data.frame(GeneName = gene, GO_terms = g[i,2]))
            }
        }
    
    
}
colnames(gg) <-c('Gene.names','term')

# 使用duplicated标识重复项的位置
duplicate_rows <- duplicated(gg$term)

# 使用逻辑索引删除重复项及其整一行
filtered_data <- gg[!duplicate_rows, ]
# 使用table函数计算每个字符的出现次数
char_counts <- table(gg$term)

# 使用sort函数按照出现次数降序排序
sorted_counts <- sort(char_counts, decreasing = TRUE)
sorted_counts[1:5]

tt <-c('cell-substrate junction','focal adhesion ','mononuclear cell differentiation','positive regulation of cell activation','regulation of peptidase activity')
ggg <-data.frame()

for (t in tt){
    for (i in 1:nrow(gg)){
        if(gg[i,2]==t){
            ggg <-rbind(ggg,data.frame(Gene.names=gg[i,1], term=gg[i,2]))
            }
        }
}

# 作go火山图
GOhs <-term_volcano(c, ggg, # 输入deg文件和go文件
             x = "log2FoldChange", y = "padj",
             label = "GeneName", label_number = 10, output = FALSE)
# 选择RColorBrewer颜色
deg_point_fill <- brewer.pal(5, "RdYlBu")
# 
GOhuoshan <-term_volcano(c, ggg,
             x = "log2FoldChange", y = "padj",
             normal_point_color = "#75aadb",  # 非显著基因散点的颜色
             deg_point_fill = deg_point_fill, # 显著基因填充色
             deg_point_color = "grey",        # 显著基因描边色
            deg_point_size =2.5,
             legend_background_fill = "#deeffc", # 图例框背景色
             label = "GeneName", label_number = 10, output = FALSE) # 显著基因标签
ggsave("./GOhuoshan.png",plot=GOhuoshan,device="png",width=8)
```

![image-20231224221400955](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/image-20231224221400955.png)
