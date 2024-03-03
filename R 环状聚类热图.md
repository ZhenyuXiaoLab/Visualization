# R 环状聚类热图

## 1、加载R包

```R
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(magrittr)
library(ggnewscale)
library(ggtreeExtra)
library(RColorBrewer)
```



## 2、生成一份可供使用的data

```R
data_frame2 <- read.csv("D:/Work/Plotplot/Circular clustering Heatmap/data.csv")
write.table(data_frame2,file = "data2.tsv",sep = "\t",row.names = FALSE)
```

> #创建示例数据框架
>
> data_frame <- data.frame(
>   Group = c("zz","hh","hh","zz"),
>   id = c(1, 2, 3, 4),  # 生物门分类信息？
>   Phylum = c("Phylum_A", "Phylum_B", "Phylum_A", "Phylum_C"),  # 样本编号？
>   Abundance = c(10, 15, 20, 12)  # 对应生物门在样本中的丰度信息)
>
> #将数据框架写入文件，以供您的代码使用
>
> write.table(data_frame, file = "your_data_file.tsv", sep = "\t", row.names = FALSE)
>



## 3、数据清洗

```R
data <- read_tsv("D:/Work/Plotplot/Circular clustering Heatmap/data2.tsv") %>% 
  select(-Group) %>%          #从数据中移除名为 "Group" 的列
  group_by(id,Phylum) %>%     #将数据按id、Phylum列进行分组
  summarise(across(where(is.numeric), ~ sum(.x, na.rm=TRUE))) %>%    
                              #对每个数值型列进行汇总操作：across 用于选择多个列，
                              #where(is.numeric) 用于筛选数值型列，sum总和计算，
                              #na.rm=TRUE 表示在计算总和时忽略缺失值。
  pivot_wider(names_from = "Phylum",values_from = "Abundance")
             #将数据从长格式转换为宽格式，其中 "Phylum" 列的唯一值将成为新的列，
             #并且 "Abundance" 列的值将填充到相应的单元格中
```



## 4、构建表达信息、树文件、分组文件

```R
##构建表达信息
exp <- data %>% pivot_longer(-id)

##构建树文件
tree <- data %>% column_to_rownames(var="id") %>% 
  dist() %>% ape::bionjs()

##构建分组文件
group <- read_tsv("D:/Work/Plotplot/Circular clustering Heatmap/data2.tsv") %>% 
  select(id,Group) %>% 
  mutate(group="group")
```



## 5、数据可视化

```R
ggtree(tree,branch.length = "none", layout = "circular",
       linetype = 1,size = 0.5, ladderize = T)+
                           #tree：树形图的输入数据。
                           #branch.length = "none"：指定树形图中不显示分支长度。
                           #layout = "circular"：设置树形图的布局为圆形。
                           #linetype = 1：指定线条类型为1。
                           #size = 0.5：设置绘图元素的尺寸为0.5。
                           #ladderize = T：指定在绘制树时对其进行阶梯化（ladderize）处理，
                                          #即对树进行重新排列，使其分支更清晰。

  layout_fan(angle =180)+      #调整开口角度

  theme(plot.margin=margin(0,1,-7,0,"cm"))+     #设置绘图的主题，包括边距的设置

  geom_tiplab(offset=7,show.legend=FALSE,size=2.5,
              color = "black",starstroke = 0)+
                           #geom_tiplab()：在树形图的叶节点处显示标签。
                           #offset = 10.5：设置标签的偏移量。🔺.
                           #show.legend = FALSE：不显示图例。
                           #size = 2.8：设置标签的大小。
                           #starstroke = 0：设置叶节点标签的边框粗细为0。

  geom_fruit(data=exp,geom=geom_tile,                 ##⭐value图例标签
             mapping=aes(y=id,x=name,fill=value),
             pwidth=0.6,offset=0.02,
             axis.params=list(axis="x",text.angle=-90,text.size=2,hjust=0))+
                           #geom_fruit()：在树形图上添加“果实”（即其他数据）。
                           #data = exp：要添加到树形图上的数据。
                           #geom = geom_tile：使用geom_tile几何图形添加数据。
                           #mapping =……设置数据映射到图形的方式，y轴对应id，x轴对应name，
                                     #填充颜色根据value的值。
                           #pwidth = 0.6：指定每个“瓷砖”（tile）的宽度。
                           #offset = 0.02：设置偏移量。
                           #axis.params = ……指定轴参数，包括轴的方向、文本角度、文本大小和水平对齐方式。

  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdBu")))+
                      #设置填充颜色的渐变色条。
                      #使用RColorBrewer包中的RdBu调色板，并翻转颜色顺序。

  new_scale_fill()+   #创建一个新的填充颜色比例尺

  geom_fruit(data=group,geom=geom_tile,                 ##⭐Group图例标签
             mapping=aes(y=id,x=group,fill=Group),color="white",
             pwidth=1,offset=0.4)+

  scale_fill_manual(values = c("#EDB749","#3CB2EC","#9C8D58"))  #手动设置填充颜色。
```



![环状聚类热图-效果示例](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/环状聚类热图-效果示例.png)

————————————————
From  https://mp.weixin.qq.com/s/p26BqddVYBwJSvf55P_QGg