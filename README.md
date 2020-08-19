# Browsing Sox2+ BTICs in Math1-Cre;SmoM2 and Math1-Cre;SmoM2;Piezo2-fl/fl SHH MB mice
#### **Welcome to Xi Huang Lab Piezo2 scRNAseq project!**

Link to Xi Huang Lab:https://lab.research.sickkids.ca/huang/

This code repository is supportive for single cell RNAseq analysis for Piezo2 project from Xi Huang Lab, Sickkids, CA. This code could be used to reproduce the scRNAseq results of Sox2+ BTICs in the paper titled as **Piezo2 governs blood-tumor-barrier and tumor quiescence depth to mask therapeutic vulnerability** that is under review. 

![Titile Figure.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Titile%20Figure.png?raw=true)



## Introduction

Medullobalstoma (MB) BTICs are mechanosensing that depends on the force-activated cation channel Piezo2. While BTICs are widely recognized as quiescent, our single cell RNAseq reveals a cryptic, two-branched BTIC trajectory that progresses from a quiescent state to two distinct cycling states. We provided the Seurat and Monocle objects to browse Sox2+ BTICs in Math1-Cre;SmoM2 and Math1-Cre;SmoM2;Piezo2-fl/fl mice. 

The raw-data of this project could be downloaded from GEO with BioProject accession ID PRJNA588007 (to be released) .  



## Tools

The pipelines for our single cell RNAseq analysis are based on Cellranger, Velocyto, SeuratWrapers, Seurat 3.0, Monocle2 and SCENIC. The links of these major pinelines are listed as below: 

​     Cellranger, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger; 

​     Velocyto, http://velocyto.org/;

​     Seurat 3.0, https://satijalab.org/seurat/; 

​     SeuratWrappers, https://github.com/satijalab/seurat-wrappers;

​     Monocle2, http://cole-trapnell-lab.github.io/monocle-release/docs/;​ 

​     SCENIC, https://github.com/aertslab/SCENIC.



## Objects

The Seurat objects of Sox2+ BTICs from Math1-Cre;SmoM2 and Math1-Cre;SmoM2;Piezo2-fl/fl could be obtained by contacting the authors.

The combined Sox2+ BTICs monocle object from both genotype could also be obtained by contacting the authors.

## Code

```R
#Browsering Sox2+ BTIC in Math1-Cre;SmoM2 or Math1-Cre;SmoM2:Piezo2-fl/fl
##Load packages
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(RColorBrewer)

##Set path to established Seurat objects of Sox2+ BTIC cells.
setwd("/mnt/hgfs/share/Xin/Revision_Upload_Objects")
load("MB_cells_Sox2_WT.Robj")
load("MB_cells_Sox2_KO.Robj")

##Switch Idents to cell groups.
Idents(MB_cells_Sox2_WT) <- "Ident"
Idents(MB_cells_Sox2_KO) <- "Ident"

###Generate UMAP plot for Sox2+ BTICs in Math1-Cre;SmoM2 or Math1-Cre;SmoM2:Piezo2-fl/fl.
p_WT <- DimPlot(MB_cells_Sox2_WT, dims = c(1,2), reduction = "umap",
                pt.size = 2.5, split.by = NULL, group.by = "Ident", 
                shape.by = NULL, order = NULL,
                label = TRUE, label.size = 4)
p_KO <- DimPlot(MB_cells_Sox2_KO, dims = c(1,2), reduction = "umap",
                pt.size = 2.5, split.by = NULL, group.by = "Ident", 
                shape.by = NULL, order = NULL,
                label = TRUE, label.size = 4)

###Plot UMAP for both genotypes
Cairo(file="UMAP_clusters.png",type="png",units="in",bg="white",width=17,height=6,pointsize=114,dpi=300)
plot_grid(p_WT,p_KO,labels = c("Math1-Cre;SmoM2","Math1-Cre;SmoM2;Piezo2-fl/fl"), label_size = 14)
dev.off()
```

![UMAP_clusters.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/UMAP_clusters.png?raw=true)

```R
##Identify Feature genes in each cluster 
###Find marker features for each cluster in Sox2+ BTIC.
MB_cells_Sox2_WT.markers <- FindAllMarkers(MB_cells_Sox2_WT,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.5,test.use = "roc",pseudocount.use = 1,verbose = TRUE)
MB_cells_Sox2_KO.markers <- FindAllMarkers(MB_cells_Sox2_KO,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.5,test.use = "roc",pseudocount.use = 1,verbose = TRUE)

###Extract Top10 feature genes of Sox2 BTIC
top10_WT <- MB_cells_Sox2_WT.markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_diff)
top10_KO <- MB_cells_Sox2_KO.markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_diff)

###Plot marker genes for each clusters in Sox2 BTIC.Top 10 genes for each clusters
p_WT <- DoHeatmap(MB_cells_Sox2_WT,features = top10_WT$gene,group.by = "Ident") + NoLegend()
p_KO <- DoHeatmap(MB_cells_Sox2_KO,features = top10_KO$gene,group.by = "Ident") + NoLegend()

###Plot heatmap of biomarkers
Cairo(file="Heatmap_of_Biomarkers_for_MB_cells_Sox2.png",type="png",units="in",bg="white",width=14,height=8,pointsize=114,dpi=300)
plot_grid(p_WT,p_KO,labels = c("Math1-Cre;SmoM2","Math1-Cre;SmoM2;Piezo2-fl/fl"))
dev.off()
```

![Heatmap_of_Biomarkers_for_MB_cells_Sox2.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Heatmap_of_Biomarkers_for_MB_cells_Sox2.png?raw=true)

```R
###Dot plot of marker genes in Sox2 BTIC.Top 10 genes for each clusters
Gene_list_WT <- unique(top10_WT$gene)
Gene_list_KO <- unique(top10_KO$gene)

#### Math1-Cre;SmoM2
p_WT <- DotPlot(MB_cells_Sox2_WT,features = as.character(Gene_list_WT),scale.by = 'radius',
        col.min = 0) +
  RotatedAxis() +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradient2(midpoint=0.8,low="black",mid="orange",high="red",space="grey")

#### Math1-Cre;SmoM2;Piezo2-fl/fl
p_KO <- DotPlot(MB_cells_Sox2_KO,features = as.character(Gene_list_KO),scale.by = 'radius',
                col.min = 0) +
  RotatedAxis() +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradient2(midpoint=0.8,low="black",mid="orange",high="red",space="grey")

#### Figure plotting
Cairo(file="Dotplot_of_each_cluster_in_MB_cells_Sox2.png",type="png",units="in",bg="white",width=12,height=10,pointsize=14,dpi=300)
plot_grid(p_WT,p_KO,labels = c("Math1-Cre;SmoM2","Math1-Cre;SmoM2;Piezo2-fl/fl"),nrow = 2)
dev.off()
```

![Dotplot_of_each_cluster_in_MB_cells_Sox2.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Dotplot_of_each_cluster_in_MB_cells_Sox2.png?raw=true)

```R
###Plot feature gene expression on UMAP. Customize your target genes in gene_list.
gene_list <-c("Nhlh1","Sstr2","Tubb3","Top2a","Mki67","Cdk1","Cenpa","Ung","Mcm2","Pcna","Olig2","Cdc20")

#### Math1-Cre;SmoM2
p_WT <- FeaturePlot(MB_cells_Sox2_WT, features = gene_list,pt.size = 1,cols = c("grey","Red"),order = TRUE, min.cutoff = NA, max.cutoff = NA, reduction = NULL,split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, blend.threshold = 1, label = FALSE, label.size = 2.5,repel = FALSE, ncol = NULL, combine = TRUE, coord.fixed = FALSE,
by.col = TRUE, sort.cell = FALSE) 

#### Math1-Cre;SmoM2;Piezo2-fl/fl
p_KO <- FeaturePlot(MB_cells_Sox2_KO, features = gene_list,pt.size = 1,cols = c("grey","Red"),order = TRUE, min.cutoff = NA, max.cutoff = NA, reduction = NULL,split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, blend.threshold = 1, label = FALSE, label.size = 2.5,repel = FALSE, ncol = NULL, combine = TRUE, coord.fixed = FALSE, by.col = TRUE, sort.cell = FALSE) 

#### Math1-Cre;SmoM2 plotting
Cairo(file="Feature_plot_of_MB_cells_Sox2_WT.png",type="png",units="in",bg="white",width=16,height=10,pointsize=14,dpi=300)
p_WT
dev.off()
```

![Feature_plot_of_MB_cells_Sox2_WT.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Feature_plot_of_MB_cells_Sox2_WT.png?raw=true)

```R
#### Math1-Cre;SmoM2;Piezo2-fl/fl plotting
Cairo(file="Feature_plot_of_MB_cells_Sox2_KO.png",type="png",units="in",bg="white",width=16,height=10,pointsize=14,dpi=300)
p_KO
dev.off()
```

![Feature_plot_of_MB_cells_Sox2_KO.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Feature_plot_of_MB_cells_Sox2_KO.png?raw=true)



```R
##RNAvelocity of Sox2+ BTICs
###Run RNA velocity to fit the gamma value.
MB_cells_Sox2_WT <- RunVelocity(object = MB_cells_Sox2_WT,deltaT = 1, kCells = 7, 
                             fit.quantile = 0.02, ncores = 14) 
MB_cells_Sox2_KO <- RunVelocity(object = MB_cells_Sox2_KO,deltaT = 1, kCells = 7, 
                                fit.quantile = 0.02, ncores = 14) 
###set color for clusters
ident.colors_WT <- c("Red","#A5B041","#16FF32","#941A41","#3283FE","#FE00FA")
names(x = ident.colors_WT) <- levels(x = MB_cells_Sox2_WT)
cell.colors_WT <- ident.colors_WT[Idents(object = MB_cells_Sox2_WT)]
names(x = cell.colors_WT) <- colnames(x = MB_cells_Sox2_WT)
ident.colors_KO <- c("#941A41","Red","#FE00FA","#16FF32","#3283FE")
names(x = ident.colors_KO) <- levels(x = MB_cells_Sox2_KO)
cell.colors_KO <- ident.colors_KO[Idents(object = MB_cells_Sox2_KO)]
names(x = cell.colors_KO) <- colnames(x = MB_cells_Sox2_KO)

###Embedding UMAP to RNA velocity.Regenerate of Figure 3e.
#### Math1-Cre;SmoM2
Cairo(file="Velocity_UMAP_MB_Sox2_WT.png",type="png",units="in",bg="white",width=5,height=5,pointsize=5,dpi=300)
show.velocity.on.embedding.cor(emb = Embeddings(object = MB_cells_Sox2_WT,reduction = "umap"),vel = Tool(object = MB_cells_Sox2_WT, slot = "RunVelocity"),n = 200, scale = "sqrt",cell.colors = ac(x = cell.colors_WT, alpha = 1.0), cex = 3, arrow.scale = 1.5, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 50, arrow.lwd = 1.25,do.par = FALSE, cell.border.alpha = 0,n.cores = 14)
dev.off()
```

![Velocity_UMAP_MB_Sox2_WT.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Velocity_UMAP_MB_Sox2_WT.png?raw=true)



```R
#### Math1-Cre;SmoM2;Piezo2-fl/fl
Cairo(file="Velocity_UMAP_MB_Sox2_KO.png",type="png",units="in",bg="white",width=5,height=5,pointsize=5,dpi=300)
show.velocity.on.embedding.cor(emb = Embeddings(object = MB_cells_Sox2_KO,reduction = "umap"),vel = Tool(object = MB_cells_Sox2_KO, slot = "RunVelocity"),n = 200,scale = "sqrt",cell.colors = ac(x = cell.colors_KO, alpha = 1.0),cex = 3, arrow.scale = 1.5, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 50,arrow.lwd = 1.25,do.par = FALSE, cell.border.alpha = 0,n.cores = 14)
dev.off()
```

![Velocity_UMAP_MB_Sox2_KO.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Velocity_UMAP_MB_Sox2_KO.png?raw=true)

```R
###Expression, phase portrait, u residual (Velocity) of Mki67, Regenerate of Figure 3f and Ext.Data.Fig 6b.
gene <- "Mki67"
####Math1-Cre;SmoM2
Cairo(file="Mki67_Velocity_Sox2_MB_WT.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(MB_cells_Sox2_WT,slot = "data", assay = "spliced"),GetAssayData(MB_cells_Sox2_WT,slot = "data", assay = "unspliced"),deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,cell.emb = Embeddings(MB_cells_Sox2_WT, "umap"), cell.colors = cell.colors_WT,old.fit = Tool(MB_cells_Sox2_WT, slot = "RunVelocity"),diagonal.quantiles = FALSE,show.gene = gene, expression.gradient = brewer.pal(9,"Greens"),residual.gradient = NULL)
```

![Mki67_Velocity_Sox2_MB_WT.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Mki67_Velocity_Sox2_MB_WT.png?raw=true)

```R
####Math1-Cre;SmoM2;Piezo2-fl/fl
Cairo(file="Mki67_Velocity_Sox2_MB_KO.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(MB_cells_Sox2_KO,slot = "data", assay = "spliced"),GetAssayData(MB_cells_Sox2_KO,slot = "data", assay = "unspliced"),deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,cell.emb = Embeddings(MB_cells_Sox2_KO, "umap"), cell.colors = cell.colors_KO,old.fit = Tool(MB_cells_Sox2_KO, slot = "RunVelocity"),diagonal.quantiles = FALSE,show.gene = gene, expression.gradient = brewer.pal(9,"Greens"),residual.gradient = NULL)
```

![Mki67_Velocity_Sox2_MB_KO.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Mki67_Velocity_Sox2_MB_KO.png?raw=true)

```R
###Comparison of Mki67 velocity pattern.(recommed reading to https://jef.works/blog/2020/01/14/rna_velocity_analysis_tutorial_tips/)
####Math1-Cre;SmoM2
xi <- median(Embeddings(object = MB_cells_Sox2_WT,reduction = "umap")[,1])
yi <- median(Embeddings(object = MB_cells_Sox2_WT,reduction = "umap")[,2])
pseudotime = atan2(Embeddings(object = MB_cells_Sox2_WT,reduction = "umap")[,2]-yi, 
                   Embeddings(object = MB_cells_Sox2_WT,reduction = "umap")[,1]-xi)/pi*180
pseudotime <- pseudotime-min(pseudotime)
pseudotime <- pseudotime/max(pseudotime)
gs <- c("Mki67")
RNAvelocity <- t(MB_cells_Sox2_WT@tools$RunVelocity$deltaE)
RNAvelocity <- RNAvelocity[,gs]
RNAvelocity_combined <- cbind(RNAvelocity,pseudotime,MB_cells_Sox2_WT@meta.data)
p_WT <- ggplot(data=RNAvelocity_combined,aes(x=pseudotime,y=RNAvelocity)) +
  geom_point(aes(colour = factor(Ident)),size = 2.5) + 
  geom_hline(yintercept = 0,color = "Grey") +
  scale_color_manual(values = ident.colors_WT) +
  geom_smooth(method = 'loess',formula = y~x,color = "black",span = 0.2,size = 0.5,linetype = "twodash") +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Streched trajectory",y="Mki67 u residuals")

####Math1-Cre;SmoM2;Piezo2-fl/fl
xi <- median(Embeddings(object = MB_cells_Sox2_KO,reduction = "umap")[,1])
yi <- median(Embeddings(object = MB_cells_Sox2_KO,reduction = "umap")[,2])
pseudotime = -atan2(Embeddings(object = MB_cells_Sox2_KO,reduction = "umap")[,2]-yi, 
                   Embeddings(object = MB_cells_Sox2_KO,reduction = "umap")[,1]-xi)/pi*180
pseudotime <- pseudotime-min(pseudotime)
pseudotime <- pseudotime/max(pseudotime)
gs <- c("Mki67")
RNAvelocity <- t(MB_cells_Sox2_KO@tools$RunVelocity$deltaE)
RNAvelocity <- RNAvelocity[,gs]
RNAvelocity_combined <- cbind(RNAvelocity,pseudotime,MB_cells_Sox2_KO@meta.data)
p_KO <- ggplot(data=RNAvelocity_combined,aes(x=pseudotime,y=RNAvelocity)) +
  geom_point(aes(colour = factor(Ident)),size = 2.5) + 
  geom_hline(yintercept = 0,color = "Grey") +
  scale_color_manual(values = ident.colors_KO) +
  geom_smooth(method = 'loess',formula = y~x,color = "black",span = 0.2,size = 0.5,linetype = "twodash") +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Streched trajectory",y="Mki67 u residuals")

####Plot MKi67 velocity along trajectory
Cairo(file="Mki67_velocity_comparison_in_MB_cells_Sox2.png",type="png",units="in",bg="white",width=12,height=8,pointsize=14,dpi=300)
plot_grid(p_WT,p_KO,labels = c("Math1-Cre;SmoM2","Math1-Cre;SmoM2;Piezo2-fl/fl"),nrow = 2)
dev.off()
```

![Mki67_velocity_comparison_in_MB_cells_Sox2.png](https://github.com/SiyiWanggou/Single-Cell-RNAseq-Code-for-Piezo2-Project/blob/master/results/Mki67_velocity_comparison_in_MB_cells_Sox2.png?raw=true)
