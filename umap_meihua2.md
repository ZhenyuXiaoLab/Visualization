# umap图美化--示例数据是pbmc
## 常规pipeline
```
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
getwd()
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc=CreateSeuratObject(pbmc.data,min.cells = 3,min.features = 200,project = 'pbmc')
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#通过设置only.pos为true来只找上调的基因
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```
![1AEE18E8-71FA-4E4F-ACAC-2D486E7BF5FD](https://user-images.githubusercontent.com/80447904/201342923-97933883-c128-4174-87ec-5b522f947016.png)
## 准备工作
```
#配色
col=unique(c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0", "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE", "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887"))#40个颜色
#legend信息
legend_txt=c('')
for (i in 1:length(levels(pbmc@active.ident)))
{
    legend_txt[i]=paste0(i,'.',levels(pbmc@active.ident)[i])
 }
names(legend_txt)=levels(pbmc)
pbmc=RenameIdents(pbmc,legend_txt) 
pbmc$legend_txt=pbmc@active.ident
#label信息
label_txt=c('')
for (i in 1:length(levels(pbmc@active.ident)))
{
    label_txt[i]=as.character(i)
 }
names(label_txt)=levels(pbmc)
pbmc=RenameIdents(pbmc,label_txt) 
pbmc$label_txt=pbmc@active.ident
#umap坐标
umap_df=Embeddings(pbmc,reduction = 'umap')
umap = umap_df %>% 
  as.data.frame() %>% cbind(legend_txt= pbmc$legend_txt)
umap$label_txt=pbmc$label_txt
#label的umap坐标信息
cell_type_med <- umap %>%
 group_by(label_txt) %>%
 summarise(
   UMAP_1 = median(UMAP_1),
   UMAP_2 = median(UMAP_2)
 )
```
## 图1
```
p1=ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = legend_txt))+geom_point(size = 1.5 , alpha =1 ,stroke=0)  +  scale_color_manual(values = col)
p1
```
<img width="412" alt="image" src="https://user-images.githubusercontent.com/80447904/209803132-d41cbf9d-9ee7-4434-8f85-b976d2ce9096.png">

## 去掉多余的东西，比如坐标轴
```
p2=p1+theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
 p2       
```
<img width="396" alt="image" src="https://user-images.githubusercontent.com/80447904/209803246-deab9fbf-aa14-4c3c-9a53-8b68bfa449c7.png">


## 设置图例legend
```
p3=p2+theme(
          legend.title = element_blank(), #去掉legend.title 
          legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=12), #设置legend标签的大小
        legend.key.size=unit(.5,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=3)))#调整legend的可以的大小也就是图例的小圆点
p3
```
<img width="395" alt="image" src="https://user-images.githubusercontent.com/80447904/209803325-bc049f6f-730a-4aca-b3fd-735396c401fc.png">

## 设置label
```
p4=p3+geom_text(aes(x=UMAP_1,y=UMAP_2,label=label_txt),data = cell_type_med,size=5,color='black')
p4
```
<img width="394" alt="image" src="https://user-images.githubusercontent.com/80447904/209803653-57362e47-a1db-4423-bf6f-20af24f5b30f.png">

## 添加坐标轴
```
p5=p4+geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +6, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 6),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +3, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 6, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 3, label = "UMAP_2",
           color="black",size = 6, fontface="bold" ,angle=90)
p5
```
<img width="402" alt="image" src="https://user-images.githubusercontent.com/80447904/209803491-a30acc87-5f40-4ba9-a44f-ef2ec9cb134b.png">

## 旋转
```
pbmc1=pbmc
angle <- 180  # 逆时针旋转180
rotation_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2)  # 构造旋转矩阵
pbmc1@reductions$umap@cell.embeddings <- pbmc1@reductions$umap@cell.embeddings %*% rotation_matrix

colnames(pbmc1@reductions$umap@cell.embeddings)=c('UMAP_1','UMAP_2')

p1=DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5,cols = allcolour,pt.size = 0.5) + NoLegend()
p2=DimPlot(pbmc1, reduction = "umap", label = TRUE, label.size = 5,cols = allcolour,pt.size = 0.5) + NoLegend()

options(repr.plot.width=12,repr.plot.height=8)
p1|p2
```
![image](https://github.com/zrz-echo/scRNA-note/assets/80447904/87bc948f-2954-4f35-8cf5-0f8017120b51)





