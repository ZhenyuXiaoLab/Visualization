# R表达分布气泡图

From： https://mp.weixin.qq.com/s/w2uhNWR3M2JXItg7yrRVUw

```R
setwd("D:/Work/Plotplot/表达分布气泡图")

#相关R包载入：
library(Seurat)
library(dplyr)
library(cols4all)
library(ggplot2)

#单细胞测试数据读入：
seurat_object <- readRDS("……")
```



## ⭐ DotPlot表达分布气泡图绘制

#### 1、数据处理

```R
#上调基因提取：avg_log2FC > 0?  
up_diff <- FindAllMarkers(seurat_object,
                          logfc.threshold = 0.6,    #logfc > 0.6
                          min.pct = 0.25,
                          only.pos = T)
head(up_diff)

up_top5 <- my_data %>%        #提取上调Top5
  group_by(cluster)%>%
  top_n(n = 5, wt = avg_log2FC)
head(up_top5)

#提取up_top5数据框gene列的唯一基因，并展示其数量以及前几个基因的值
genes <- unique(up_top5$gene)
length(genes);head(genes)
```



#### 2、DotPlot函数绘图（Seurat自带函数）

```R
DotPlot(my_data, 
               features = genes) #或自行准备目标基因列表


#图表美化：
##ggplot2对象，相关主题函数可直接叠加

mycolor <- c("darkblue","purple","orange")   #c("963c59","white","#5597ad")

options(repr.plot.width=8, repr.plot.height=6)
p0 <- DotPlot(seurat_object,
              features = genes) +
  scale_color_gradientn(colors = mycolor)+                        #绘制渐变色
  #scale_color_continuous('linear_yl_mg_bu',reverse = F) +        #配色自定义
  theme(axis.text.x = element_text(angle = 60, hjust = 1))        #x轴标签旋转60°
p0
```



![图3](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/图3.png)



# ⭐gglot2达分布气泡图绘制

#### 1、数据处理

##### ①

```R
#提取目标基因表达量数据：
dt <- seurat_object@assays$RNA@data[genes,] #这里以Top5上调为例
dtt <- expm1(dt)   #因为RNA-seq数据中的表达值通常是进行对数变换的,
                   #使用 expm1可以将其还原为非对数形式。

#根据celltype拆分数据：
dt_list <- lapply(unique(seurat_object@meta.data$seurat_clusters), FUN = function(x){
  cellid <- colnames(seurat_object)[which(seurat_object@meta.data$seurat_clusters == x)]
                   #找到属于当前细胞类型 x 的所有细胞的列名（cellid）
  data_cluster <- dtt[,cellid]    #新的矩阵：表达量--细胞类型（这里是cluster）
  return(data_cluster)   #形成一个列表，其中每个元素对应于一个细胞类型的表达量数据
})
names(dt_list) <- unique(seurat_object@meta.data$seurat_clusters)


if(F){
dt_list <- lapply(unique(seurat_object$celltype), FUN = function(x){
  cellid <- colnames(seurat_object)[which(seurat_object$celltype == x)]
  data_cluster <- dtt[,cellid]
  return(data_cluster)
})
names(dt_list) <- unique(seurat_object$celltype)
} ##这段是教程里的
```



##### ②

```R
##计算cluster每个基因的平均表达量和表达基因的细胞比例
#计算均值和表达细胞比例:
dt_sum <- list()

for(i in 1:length(dt_list)){
  means <- log1p(rowMeans(dt_list[[i]]))   #均值
                  #log1p处理，确保数据不会受到取对数时小于等于零的问题的干扰

  pro <- apply(dt_list[[i]], 1, FUN = function(x){
    pro <- sum(x > 0.5)/length(x) #表达量大于0.5的细胞视为表达基因
  })

  gene <- factor(rownames(dt_list[[i]]), levels = genes)
  celltype <- rep(names(dt_list)[i], nrow(dt_list[[i]]))
  stat <- data.frame(gene, celltype, means, pro)
  dt_sum[[i]] <- stat
}

dot_dt <- do.call("rbind", dt_sum)
head(dot_dt)        #将列表中的数据框按行绑定成一个大的数据框 dot_dt，
                    #其中包含了所有细胞类型的每个基因的平均表达量和表达基因的细胞比例
```



##### ③

```R
##归一化：    
#lapply函数对基因列表genes进行循环，dot_scale最终是包含了每个基因的归一化数据的列表
dot_scale <- lapply(genes, FUN = function(x){     
  dt1 <- dot_dt[which(dot_dt$gene == x),]        
                     #从整个数据框 dot_dt中筛选出当前基因的行
  dt1$scale <- as.vector(scale(dt1$means))        
                     #scale函数对当前基因的均值归一化，结果转变为向量储存在scale列
  return(dt1)        #返回处理后的数据框dt1
})

#得到ggplot2绘图所需长数据
dot_dt <- do.call("rbind", dot_scale)
head(dot_dt)     
      #将归一化后的数据框列表 dot_scale 合并为一个大的数据框 dot_dt。
      #这样做是为了方便后续的绘图，因为ggplot2通常需要长格式的数据

#celltype转换为因子，固定顺序
dot_dt$celltype <- factor(dot_dt$celltype,
                          levels = c('Naive CD4 T','Memory CD4 T',.....))
```



#### 2、ggplot2绘图

```R
p1 <- ggplot(dot_dt, aes(x = celltype, y = gene)) + #建立映射
  geom_point(aes(size = pro, color = scale)) #绘制散点
p1

#图表个性化调整:
options(repr.plot.width = 10,repr.plot.height = 10)
p2 <- p1 +
  scale_size_continuous(range = c(0,6)) + #气泡大小范围调整
  scale_color_gradient2(low='#5E3C99', high='#E66101', mid = 'white') + #配色自定义
  theme_classic() +   #使用经典的主题
  theme(legend.text = element_text(size = 14), #图例文本字号
        legend.title = element_text(size = 16), #图例标题字号
        axis.text = element_text(size = 14), #坐标轴标签字号
        axis.title = element_text(size = 16), #坐标轴标题字号
        axis.text.x = element_text(angle = 60, hjust = 1)) +#x轴标签旋转60°
  labs(x = 'Identity',y = 'Features', #xy轴标题修改
       color = 'Average Expression', size = 'Percent Expressed') + #图例标题修改
  coord_flip() #坐标轴翻转
p2
```

![气泡图2](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/气泡图2.png)





## 注释色块  Adobe Illustrator/Photoshop ?

![图5](https://grp-share-code.obs.cn-north-4.myhuaweicloud.com/图5.png)