# R 相关性热图

##From: https://blog.csdn.net/dege857/article/details/132097240

```R
setwd("D:/Work/Plotplot/相关性热图")

library(tidyverse)
library(linkET)       ##devtools::install_github("Hy4m/linkET", force = TRUE)
library(RColorBrewer)
library(ggtext)
library(magrittr)
library(psych)
library(reshape)

##⭐导入数据
table1 <- read.csv("D:/Work/Plotplot/相关性热图/数据文件/varespec.csv",sep=',',header=TRUE)
table2 <- read.csv("D:/Work/Plotplot/相关性热图/数据文件/varechem.csv",sep=',',header=TRUE)
   #varespec数据框：24行，44列（44个物种的估计覆盖值）
   #varechem数据框：24行，14列（与varespec中相同地点的土壤特性）
```



## linkET函数介绍

```R
##correlate函数可以计算数据的相关性
correlate(table1)
#也可以计算不同数据的相关性系数
correlate(table1[1:30], table2)
#算出了相关性就能生成图形了，先来单个数据的
correlate(table2) %>% 
  as_md_tbl() %>% 
  qcorrplot() +
  geom_square()
#再来个双数据的
correlate(table1[1:30], table2) %>% 
  qcorrplot() +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
```

```R
##Qcorrplot函数能把系数化的矩阵图示化
#Type系数可以控制我们取局部图形，比如我只想取下半截
qcorrplot(correlate(table2), type = "lower") +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))


##使用set_corrplot_style()函数和geom_square()函数还可以定制自己的风格
#比如我想改方框为圆形，定义颜色为红色、蓝色和白色

set_corrplot_style(colours = c("red", "white", "blue"))
qcorrplot(correlate(table2), type = "lower") +
  geom_shaping(marker = "circle")

##想要从新回复成系统自定义颜色可以使用
set_default_style()
```



## 相关性热图绘制

#### 1、初步绘制

```R
#绘制前要进行一个曼特尔试验，R包作者是这样说的，相异矩阵的Mantel和偏Mantel检验。
                   #（介绍见后）
#注意：spec_select选择分类的是列的引索值，varespec数据刚好44列
mantel <- mantel_test(table1,table2,
                      spec_select = list(Spec01 = 1:7,
                                         Spec02 = 8:18,
                                         Spec03 = 19:37,
                                         Spec04 = 38:44))

#得出每个类别的R值和P值后我们对他们进行分段表示
mantel <- mantel %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

#计算好后就可以进一步绘图了
qcorrplot(correlate(table2), type = "lower", diag = FALSE) +
              #diag = FALSE：将对角线上的元素设为FALSE → 不显示对角线上的相关性
  geom_square() +                 #添加正方形到图中，每个正方形代表一个相关系数
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature())     #指定连接线的曲率
```

> **曼特尔试验**
>
> ​        ① 相关系数（Mantel's r）： 衡量两个矩阵之间的相关性程度。取值范围在-1到1之间，其中1表示完全正相关，-1表示完全负相关，0表示无相关性。
> ​        ② p值（Mantel's p）： 是用来检验Mantel's r值是否在零假设下显著的概率。p值小于显著性水平（通常是0.05）时，我们可以拒绝零假设，认为两个矩阵之间存在显著的相关性。
>
> ​        总体来说，Mantel's r值告诉我们两个矩阵之间的关联程度，而p值告诉我们这种关联是否是由于随机因素引起的。如果p值很小，我们通常会认为两个矩阵之间存在显著的关联。



#### 2、修改参数美化

① 版本1

```R
qcorrplot(correlate(table2), type = "lower", diag = FALSE) +
  geom_square() +    
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +  

  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +

  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),             #order指图例出现的顺序
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3),
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 
```

![图1](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/图1.png)



② 版本2：

（From： https://mp.weixin.qq.com/s/SWh5rSRPy-eANw69x9hAdg）

```R
qcorrplot(correlate(table2,method = "spearman"),diag=F,type="upper")+
  geom_tile()+                     ##每一格都是填满的
  geom_mark(size=2.5,sep="\n")+    ##内部标记
  geom_couple(aes(colour=pd,size=rd),data=mantel,label.colour = "black",
              curvature=nice_curvature(0.15),
              nudge_x=0.2,      #在x轴上微调标签的位置(水平向右偏移了0.2个单位)
              label.fontface=2, #定义文字标签的字体粗细，设置为2，可能表示粗体
              label.size =4,
              drop = T)+        #是否删除数据中的缺失值。这里设置为 T，表示删除

  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdBu")))+  ##颜色翻转
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values =c("#D95F02","#1B9E77","#A2A2A288")) +
  guides(size = guide_legend(title = "cor",override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "P_value",override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "spearman's r",order = 3))+

  theme(plot.margin = unit(c(0,0,0,-1),units="cm"),   ##绘图区域边距（上右下左）
               panel.background = element_blank(),           ##绘图区域背景（无）
               axis.text=element_markdown(color="black",size=10),   #坐标轴刻度文本样式
               legend.background = element_blank(),          ##图例背景
               legend.key = element_blank())                 ##设置图例的键（legend key）

##🔺.导出高清图片？？
ggsave("D:/mantel-linkET.tiff",p1,width=8,height=6)
```

![图3](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/图3.jpg)



其他教程：https://blog.csdn.net/qq_35294674/article/details/130950109?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-5-130950109-blog-132097240.235