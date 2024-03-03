# 第一步：预处理

### 1、所需包的安装和载入

```R
library("ComplexHeatmap")  ##BiocManager::install("ComplexHeatmap")
library("circlize")        ##install.packages("circlize")
```



### 2、设置工作目录和载入本地数据

```R
setwd("D:/Work/Plotplot/Circular Heatmap/")
data <- read.csv("D:/Work/Plotplot/Circular Heatmap/example_data.CSV",
                 header=T,
                 row.names = 1)   ##将CSV文件中的第一列作为数据框的行名
head(data)
```


> head(data)
>       Sample1   Sample2  Sample3    Sample4   Sample5
> TSPAN6   3.6652797 4.8182373 5.900849 2.46545665 4.7428646
> TNMD     0.2568738 1.8581361 1.110583 0.02039193 0.3249350
> DPM1     5.6631558 6.1038026 6.441289 5.32127896 4.6557647
> SCYL3    1.7349495 1.8291205 1.476910 2.06575159 2.0488829
> C1orf112 1.2511585 1.8115178 1.539492 1.91868000 0.9876195
> FGR      0.9622511 0.9052452 1.362060 2.71399279 1.7190058



### 3、转化为矩阵并对其进行标准化

```R
madt <- as.matrix(data)
madt2 <- t(scale(t(madt)))      
          ## t():转置矩阵（行列互换）  用了两次，转置再还原
          ## scale():标准化数据→确保每一列的值具有标准正态分布(均值为0，标准差为1)

Heatmap(madt2)    # 默认参数绘制普通热图

```

![微信图片_20231123095428](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/微信图片_20231123095428.png)



### 4、计算数据大小范围，重新定义热图颜色梯度

```R
range(madt2)
mycol <- colorRamp2(c(-1.7, 0, 1.7),c("blue", "white", "red"))
```

> range(madt2)
> [1] -1.667386  1.666624



# 第二步：circos.heatmap() 绘制环形热图

### 1、绘制基础环形热图

```R
circos.heatmap(madt2,col=mycol)
circos.clear()   ##🔺.绘制完成后需要使用此函数完全清除布局
```



### 2、添加参数进行美化

#### ① 添加聚类树

```R
circos.par(gap.after=c(50))   #c(50) 意味着在每个轨道之后添加一个宽度为50的间隔
                              #调整圆环首尾间的距离，数值越大，距离越宽

circos.heatmap(madt2,col=mycol,
               dend.side="inside",        #dend.side：控制行聚类树的方向
               rownames.side="outside",
               rownames.col="black",
               rownames.cex=0.9,
               rownames.font=1,
               cluster=TRUE)              #cluster=TRUE为对行聚类
circos.clear()

```

#### ② 聚类树的美化

```R
library(dendextend)  ##install.packages("dendextend")  改颜色
library(dendsort)    ##install.packages("dendsort")    聚类树回调

circos.par(gap.after=c(50))
circos.heatmap(madt2,col=mycol,
               dend.side="inside",rownames.side="outside",
               track.height = 0.25,       #轨道的高度，数值越大圆环越粗
               rownames.col="black",
               rownames.cex=0.9,
               rownames.font=1,
               cluster=TRUE,
               dend.track.height=0.15,    #调整行聚类树的高度
               dend.callback=function(dend,m,si) {
                 color_branches(dend,k=15,col=1:15)  #修改聚类树颜色
               })
```

> 参数解释
>     #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用
>
> ​    三个参数：dend--当前扇区的树状图(是一个树状图对象，表示要进行颜色处理的树状图)；
>
> ​                       m--当前扇区对应的子矩阵(这可能是与树状图有关的其他参数或对象)；
>
> ​                       si--当前扇区的名称同样，可能是与树状图有关的其他参数或对象).
>
> ​    #color_branches 函数用于给树状图的分支上色，其中的参数包括：
>
> ​                       dend: 要处理的树状图对象。
>
> ​                       k: 表示颜色的数量或分组数量。在这里是15。
>
> ​                       col: 指定用于上色的颜色向量。这里使用了1到15的颜色。

#### ③ 添加图例标签

```R
lg <- Legend(title="Exp",col_fun = mycol,direction = c("vertical"))
grid.draw(lg)
   #"horizontal": 条目水平排列。"vertical": 条目垂直排列
   #"auto": 自动选择水平或垂直排列，具体取决于图例的大小和内容
```

#### ④ 添加列名

```R
circos.track(track.index=get.current.track.index(),      # 获取当前轨道的索引
             bg.border=NA,                   # 设置背景边界颜色为透明
             panel.fun=function(x,y){        # panel.fun:用于在轨道上绘制内容
                if(CELL_META$sector.numeric.index==1){
                   cn=colnames(madt2)
                   n=length(cn)
                   circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.8,"mm"),#x坐标
                   4.5+(1:n)*1.25,   #y坐标
                   cn,cex=0.8,adj=c(0,1),facing="inside")
                }
             })   
circos.clear()
```

> #error🔺.Note: 5 points are out of plotting region in sector 'group', track '3'.
>
> ##circos.track参数解释
>
> ​          panel.fun: 在轨道上绘制内容。这里它被用于在第一个扇区的轨道上添加文本标签。
>
> ​          CELL_META$sector.numeric.index: 表示当前扇区的索引。
>
> ​              #在这里，代码检查是否是第一个扇区。……👇
>
> ​          colnames(madt2): 获取数据矩阵 madt2 的列名。
>
> ​          circos.text: 在环状图上添加文本标签。在这里，它用于添加列名。
>
> ​          convert_x(0.8, "mm"): 将0.8毫米的值从相对坐标转换为绝对坐标。
>
> ​          5是一个基础的 y 坐标，(1:n) * 1.1 创建一个包含从 1 到 n 的整数的序列，
>
> ​                                             并将每个整数乘以 1.1。这样可以在垂直方向上平均分布文本标签。
>
> #cn: 包含列名的向量，这些列名将作为文本标签添加到环状图上。
>
>   cex=0.8: 控制文本的缩放比例，这里设置为 0.8 表示文本的大小为默认大小的 80%。
>
>   adj=c(0, 1): 控制文本标签的对齐方式。在这里表示文本标签的水平对齐方式为左对齐（0），垂直对齐方式为顶部对齐（1）。
>
>   facing="inside": 指定文本标签朝向的方向。"inside"表示文本标签朝向环状图的内部
>
> #🔺.CELL_META$cell.xlim[2]表示从CELL_META对象中提取名为cell.xlim的列表的第二个元素，
>                       #这个元素可能是一个包含有关当前单元格 x 坐标限制的信息的列表。


![微信图片_20231123101040](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/微信图片_20231123101040.png)

### 3、修改配色完整走一遍

```R
#这里代码和上文相同，仅改变了颜色和circos.par（圆环首位的距离）
mycol2 <- colorRamp2(c(-1.7, 0, 1.7),c("#57ab81", "white", "#ff9600"))

circos.par(gap.after=c(50))
circos.heatmap(madt2,col=mycol2,dend.side="inside",rownames.side="outside",
               track.height = 0.25,
               rownames.col="black",
               rownames.cex=0.9,
               rownames.font=1,
               cluster=TRUE,
               dend.track.height=0.15,      ##🔺.聚类树最里面有点重叠？？
               dend.callback=function(dend,m,si) {
                 color_branches(dend,k=15,col=1:15)
               }
)

lg=Legend(title="Exp",col_fun=mycol2,direction = c("vertical"))
grid.draw(lg)

circos.par(points.overflow.warning=FALSE)

circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  if(CELL_META$sector.numeric.index==1){
    cn=colnames(madt2)
    n=length(cn)
    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.8,"mm"),#x坐标
                4.5+(1:n)*1.255,#y坐标
                cn,cex=0.8,adj=c(0,1),facing="inside")
  }
},bg.border=NA)

  ##问题Note: 5 points are out of plotting region in sector 'group', track '3'.
  ##解决？circos.par(points.overflow.warning=FALSE) 屏蔽错误？
circos.clear()
```

![1700705672203](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/1700705672203.png)

————————————————
From  https://mp.weixin.qq.com/s/PZgLoUtpVj9e7iBRZkregQ