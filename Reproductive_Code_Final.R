#-----Monocle3-alpha and Linux------
#----Charpter 1. Loading the scRNAseq dataset----
setwd("/mnt/hgfs/share18/share1/Xin/analysis/Reproduce")
#loading the data of MS(Wt) and MSP2(Knockout)
library(Matrix)
matrix_dir = "/mnt/hgfs/share18/share1/Xin/analysis/Reproduce/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "features.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
matrix <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(matrix) = barcode.names$V1
rownames(matrix) = feature.names$V1
#Loading the dataset
setwd("/mnt/hgfs/share18/share1/Xin/analysis/Reproduce")
library(monocle)
library("ggthemes")
library("RColorBrewer")
library(VGAM)
library(Cairo)
library(dplyr)
HSMM_sample_sheet_all<-read.delim("barcodes.txt",header=T,row.names = 1) 
HSMM_gene_annotation_all<-read.delim("features.txt",header=T,row.names = 1)
pd<-new("AnnotatedDataFrame",data = HSMM_sample_sheet_all)
fd<-new("AnnotatedDataFrame",data = HSMM_gene_annotation_all)
HSMM_all<-newCellDataSet(matrix,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.1, 
                         expressionFamily = negbinomial.size()) 
#Estimate Size Factors and Dispersions
HSMM_all <- estimateSizeFactors(HSMM_all)
HSMM_all <- estimateDispersions(HSMM_all)
#-----------End. Charpter 1.--------------

#-----------Charpter 2. Identify the clusters of cell population------
HSMM_all <- detectGenes(HSMM_all,min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(HSMM_all),num_cells_expressed >=10)) 
HSMM_all <- setOrderingFilter(HSMM_all,expressed_genes) 
HSMM_all <- preprocessCDS(HSMM_all, method = 'PCA', norm_method = 'log', num_dim = 50, verbose = T, pseudo_expr = 1, residualModelFormulaStr = "~num_genes_expressed")
HSMM_all <- reduceDimension(HSMM_all, max_components = 2,
                            reduction_method = 'UMAP',
                            min_dist = 0.5,
                            n_neighbors = 30,
                            verbose = T)
HSMM_all <- clusterCells(HSMM_all, method = 'densityPeak', res = 1e-4, verbose = T, num_clusters = 17)
#UMAP plot for all 17 clusters (Extended Data Figure 1a.)
Cairo(file="Plot_UMAP_17_cluster_facet.png",type="png",units="in",bg="white",width=6,height=5,pointsize=5,dpi=300)
plot_cell_clusters(HSMM_all, 1, 2, color_by = "Cluster", cell_size = 0.5,show_group_id=TRUE) + 
  #facet_wrap(~Group) +
  theme(legend.text=element_text(size=1)) + #set the size of the text
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))  #put the color legend on the right
#geom_rug()
dev.off()
#--------- End. Charpter 2.--------

#--------- Charpter 3. Identify Marker Genes for each Cluster.-------
#Generate feature genes of each cluster
spatial_res <- principalGraphTest(HSMM_all, relative_expr = TRUE, method = c("Moran_I"), k = 25, cores = detectCores() - 2, verbose = T) 
cluster_marker_res <- find_cluster_markers(HSMM_all,
                                           spatial_res,
                                           group_by = 'Cluster',
                                           morans_I_threshold = 0.25)
write.table(cluster_marker_res,"markers_of_each_sub_cluster_17_clusters.txt",sep='\t',quote=FALSE,col.names = TRUE,row.names = TRUE)
cluster_marker_res <- read.delim("markers_of_each_sub_cluster_17_clusters.txt",header=T,row.names = 1)

#Primarily screening for feature genes
genes <- cluster_marker_res %>%
  dplyr::group_by(Group) %>% dplyr::slice(which(specificity > 0.5))
write.table(genes,"markers_of_each_sub_cluster_17_05.txt",sep='\t',quote=FALSE,col.names = TRUE,row.names = TRUE)
genes <- read.delim("markers_of_each_sub_cluster_17_05.txt",header = TRUE,row.names = 1)

Cairo(file="Plot_Marker_genes_by_group_17_05.png",type="png",units="in",bg="white",width=25,height=8,pointsize=10,dpi=300)
plot_markers_by_group(HSMM_all, genes$gene_short_name, group_by = "Cluster",ordering_type = 'maximal_on_diag',
                      axis_order = 'marker_group',
                      flip_percentage_mean = TRUE) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#Plot Heatmap for feature markers. According to the marker gene expression pattern, the MB cells could be regrouped
#into three major groups.
genes <- cluster_marker_res %>%
  dplyr::filter(specificity > 0.5) %>%
  dplyr::group_by(Group) %>%
  dplyr::arrange(Group, dplyr::desc(specificity))
Cairo(file="Plot_Marker_cluster_by_gene_17_05.png",type="png",units="in",bg="white",width=12,height=3,pointsize=10,dpi=300)
plot_markers_cluster(HSMM_all,
                     as.character(genes$gene_short_name),
                     minimal_cluster_fraction = 0.01,
                     show_rownames = T,
                     group_by = "Cluster",
                     scale_max = 3,
                     scale_min = -3,
                     return_heatmap = TRUE)
dev.off()
#Annotate Cell types of each cluster
pData(HSMM_all)$assigned_cell_type <- pData(HSMM_all)$Cluster
pData(HSMM_all)$assigned_cell_type <- dplyr::recode(pData(HSMM_all)$assigned_cell_type,
                                                    "1"="SHH-MB B1",
                                                    "5"="SHH-MB B2",
                                                    "10"="SHH-MB B3",
                                                    "2"="SHH-MB A1",
                                                    "9"="SHH-MB A2",
                                                    "8"="SHH-MB C2",
                                                    "12"="SHH-MB C1",
                                                    "3"="Sprial-Ganglion-Neurons/Brain-Stem",
                                                    "6"="Neural-stem-cell/Neural-progenitors",
                                                    "7"="Macrophagy/Monocytes",
                                                    "16"="Microgial",
                                                    "4"="Non-specific",
                                                    "11"="Neruons",
                                                    "14"="Astrocytes",
                                                    "13"="Oligodenrocytes",
                                                    "15"="Erythroblast",
                                                    "17"="Endothelial cells")
#Plot Top 3 Marker Genes for each clusters (Extended Data Figure 1b.)
testgenes <-c("Atoh1","Smo","Sox2","Dxc","Golga7b","Dyncli","Phyhip","Tmem255b","Ascl1","Pcdh15","Cytip","Ccr2","Themis2",
              "Pcp2","Gng13","Tmem88b","Mog","Fa2h","Aqp4","Acsbg1","Agt","Alas2","Hba-a1","Hba-a2","C1qb","Aif1","Ctss","Sox17","Flt4","Tie1")
Cairo(file="Plot_Marker_genes_1.png",type="png",units="in",bg="white",width=12,height=6,pointsize=10,dpi=300)
plot_markers_by_group(HSMM_all, testgenes, group_by = "assigned_cell_type",ordering_type = 'cluster_row_col',
                      axis_order = 'marker_group',
                      flip_percentage_mean = FALSE) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradient2(midpoint=1.5,low="black",mid="orange",high="red",space="grey")
dev.off()

#Identify Sox2 and Atoh1 double positive cells.
Sox2_id <- row.names(subset(fData(HSMM_all),gene_short_name == "Sox2"))
Atoh1_id <- row.names(subset(fData(HSMM_all),gene_short_name == "Atoh1"))
Sox2_Atoh1_exp <- exprs(HSMM_all)[c(Sox2_id,Atoh1_id),]
Sox2_Atoh1_exp <- as.matrix(Sox2_Atoh1_exp)
Sox2_Atoh1_exp <- t(Sox2_Atoh1_exp)
colnames(Sox2_Atoh1_exp) <- c("Sox2","Atoh1")
write.table(Sox2_Atoh1_exp,"Sox2_Atoh1_exp.txt",row.names = TRUE,col.names = TRUE,quote=FALSE,sep='\t')
Sox2_Atoh1_exp <- read.delim("Sox2_Atoh1_exp.txt",header=TRUE,row.names = 1)
Sox2_Atoh1_coexpression <- subset(Sox2_Atoh1_exp,Sox2 >= 1 & Atoh1 >= 1)
pData(HSMM_all)$Sox2_Atoh1_coexpression <- rownames(pData(HSMM_all)) %in% rownames(Sox2_Atoh1_coexpression)
Cairo(file="Plot_Sox2_Atoh1_facet.png",type="png",units="in",bg="white",width=4.9,height=3.5,pointsize=5,dpi=300)
plot_cell_clusters(HSMM_all, 1, 2, color_by = "Sox2_Atoh1_coexpression", cell_size = 0.5,show_group_id=FALSE) + 
  facet_wrap(~Group) +
  scale_color_manual(values = c("grey","red")) +
  theme(legend.text=element_text(size=1)) + #set the size of the text
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))  #put the color legend on the right
#geom_rug()
dev.off()

#Identify Sox2 and Gfap double positive cells.
Sox2_id <- row.names(subset(fData(HSMM_all),gene_short_name == "Sox2"))
Gfap_id <- row.names(subset(fData(HSMM_all),gene_short_name == "Gfap"))
Sox2_Gfap_exp <- exprs(HSMM_all)[c(Sox2_id,Gfap_id),]
Sox2_Gfap_exp <- as.matrix(Sox2_Gfap_exp)
Sox2_Gfap_exp <- t(Sox2_Gfap_exp)
colnames(Sox2_Gfap_exp) <- c("Sox2","Gfap")
write.table(Sox2_Gfap_exp,"Sox2_Gfap_exp.txt",row.names = TRUE,col.names = TRUE,quote=FALSE,sep='\t')
Sox2_Gfap_exp <- read.delim("Sox2_Gfap_exp.txt",header=TRUE,row.names = 1)
Sox2_Gfap_coexpression <- subset(Sox2_Gfap_exp,Sox2 >= 1 & Gfap >= 1)
pData(HSMM_all)$Sox2_Gfap_coexpression <- rownames(pData(HSMM_all)) %in% rownames(Sox2_Gfap_coexpression)

Cairo(file="Plot_Sox2_Gfap_facet.png",type="png",units="in",bg="white",width=4.9,height=3.5,pointsize=5,dpi=300)
plot_cell_clusters(HSMM_all, 1, 2, color_by = "Sox2_Gfap_coexpression", cell_size = 0.5,show_group_id=FALSE) + 
  facet_wrap(~Group) +
  scale_color_manual(values = c("grey","red")) +
  theme(legend.text=element_text(size=1)) + #set the size of the text
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))  #put the color legend on the right
#geom_rug()
dev.off()

#Calculate proportion of Sox2-Atoh1 positive cells in Sox2 MB and Sox2 astrocytes
Sox2_expression <- subset(Sox2_Gfap_exp,Sox2 >= 1)
pData(HSMM_all)$Sox2_expression <- rownames(pData(HSMM_all)) %in% rownames(Sox2_expression)
MS <- row.names(subset(pData(HSMM_all),Group == "MS"))
HSMM_MS <- HSMM_all[,MS]
MS_MB_cells <- row.names(subset(pData(HSMM_MS),Cluster == 1 | Cluster == 2 | Cluster == 5 | Cluster == 9 | 
                               Cluster == 10 | Cluster == 8 | Cluster == 12 ))
HSMM_MS_MB_Cells <- HSMM_MS[,MS_MB_cells]
MS_astrocytes <- row.names(subset(pData(HSMM_MS),Cluster == 14))
HSMM_MS_astrocytes <- HSMM_MS[,MS_astrocytes]
table(pData(HSMM_MS_MB_Cells)$Sox2_Atoh1,pData(HSMM_MS_MB_Cells)$Sox2_expression)
table(pData(HSMM_MS_astrocytes)$Sox2_expression)
#------- End. Charpter 3------

#------- Charpter 4. Enrichment analysis of SHH-A,SHH-B and SHH-C gene signature.------
#Identification of MB cell sub-population AUCell enrichment analysis
#Identify the subgroup of MB cell population in SHH-A, SHH-B, SHH-C signature.
#According to the Nature Article --Resolving medulloblastoma cellular architecture by single-cell genomics
BiocManager::install(c("doMC", "doRNG"))
install.packages("doMC", repos="http://R-Forge.R-project.org")
BiocManager::install(c("mixtools", "GEOquery", "SummarizedExperiment"))
BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh"))
BiocManager::install("AUCell")
BiocManager::install("GSEABase")
BiocManager::install("mixtools")
install.packages("NMF")
library(AUCell)
library(GSEABase)
#Generate MB cells population from whole dataset
MB_cells <- row.names(subset(pData(HSMM_all),Cluster == 1 | Cluster == 2 | Cluster == 5 | Cluster == 9 | 
                               Cluster == 10 | Cluster == 8 | Cluster == 12 ))
HSMM_MB <- HSMM_all[,MB_cells]
write.table(pData(HSMM_MB),"HSMM_MB_col_annotation.txt",col.names = TRUE,row.names = TRUE,quote = FALSE,sep="\t") #?
#Transfer mouse gene symbol into human gene symbol
fData(HSMM_MB)$gene_id <- row.names(fData(HSMM_MB))
target_genes <- fData(HSMM_MB)[,-2:-3]
target_genes <- unique(target_genes)
Mouse_gene_id <- target_genes$gene_short_name
library(biomaRt)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
class(human)
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")
human_gene_id <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                        values = Mouse_gene_id, mart = mouse,
                        attributesL = c("hgnc_symbol"), martL = human,
                        uniqueRows = T)
colnames(human_gene_id) <- c("gene_short_name","human_symbol")
target_genes <- merge(target_genes,human_gene_id,by="gene_short_name",all=FALSE)
HSMM_MB <- HSMM_MB[target_genes$gene_id,]
#Extract gene expression matrix
df <- exprs(HSMM_MB)
df <- as.matrix(df)
rownames(df) <- target_genes$human_symbol
geneSets <- getGmt("MB_signature.gmt.txt")
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))
cells_rankings <- AUCell_buildRankings(df, nCores=12, plotStats=TRUE,verbose =TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),verbose = TRUE,
                            nCores = 10)
#Generate the histogram of AUCell enrichment analysis
Cairo(file="Plot_AUC_of_MB_cells.png",type="png",units="in",bg="white",width=3,height=3,pointsize=3,dpi=300)
par(mfrow=c(3,1)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, thrP = 0.01, nCores = 12, 
                                             smallestPopPercent = 0.25, plotHist=TRUE, assign=TRUE,
                                             nBreaks = 500) 
dev.off()
#Thresholds selection of each geneset
warningMsg <- sapply(cells_assignment,function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]
cells_assignment$'SHHA (99g)'$aucThr$thresholds
cells_assignment$'SHHB (99g)'$aucThr$thresholds
cells_assignment$'SHHC (99g)'$aucThr$thresholds
cells_assignment$'SHHA (99g)'$aucThr$selected
cells_assignment$'SHHB (99g)'$aucThr$selected
cells_assignment$'SHHC (99g)'$aucThr$selected
SHHA_assigned<-cells_assignment$'SHHA (99g)'$assignment
SHHB_assigned<-cells_assignment$'SHHB (99g)'$assignment
SHHC_assigned<-cells_assignment$'SHHC (99g)'$assignment
name_SHHA <- "SHHA (99g)"
Cairo(file="Plot_AUC_of_SHHA.png",type="png",units="in",bg="white",width=2.5,height=2,pointsize=3,dpi=300)
AUCell_plotHist(cells_AUC[name_SHHA,],aucThr = 0.147,cex.axis = 2.0,cex.lab = 2.0)
abline(v=0.147)
dev.off()
name_SHHB <- "SHHB (99g)"
Cairo(file="Plot_AUC_of_SHHB.png",type="png",units="in",bg="white",width=2.5,height=2,pointsize=3,dpi=300)
AUCell_plotHist(cells_AUC[name_SHHB,],aucThr = 0.56,cex.axis = 2.0,cex.lab = 2.0)
abline(v=0.56)
dev.off()
name_SHHC <- "SHHC (99g)"
Cairo(file="Plot_AUC_of_SHHC.png",type="png",units="in",bg="white",width=2.5,height=2,pointsize=3,dpi=300)
AUCell_plotHist(cells_AUC[name_SHHC,],aucThr = 0.27,cex.axis = 2.0,cex.lab = 2.0)
abline(v=0.27)
dev.off()
#Exploring the cell assignment, Clustering based on SHH-A,SHH-B and SHH-C gene signature enrichment score to
#reclustered the MB cells into three groups which are consistent with UMAP clustering based 3 major MB groups.
cellsAssigned <- lapply(cells_assignment,function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned,value.name="cell")
colnames(assignmentTable)[2]<-"geneSet"
head(assignmentTable)
assignmentMat <- table(assignmentTable[,"geneSet"],assignmentTable[,"cell"])
assignmentMat[,1:3]
set.seed(123)
miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),250)]
library(NMF)
Cairo(file="Plot_AUC_of_HEATMAP.png",type="png",units="in",bg="white",width=5,height=3,pointsize=3,dpi=300)
aheatmap(miniAssigMat,scale="none",color="black",legend=FALSE)
dev.off()
library(pheatmap)
cell_names <- colnames(miniAssigMat)
cell_annotation <- pData(HSMM_MB)[cell_names,]
annotation_col <- data.frame(cell_annotation$Cluster,cell_annotation$Group)
rownames(annotation_col) <- rownames(cell_annotation)
ann_colors = list(
  cell_annotation.Cluster = c("1"="#FF0000","5"="#990000","10"="#FF9999","2"="#66CC00","9"="#33FF99","8"="#CC00CC","12"="#FF00FF"),
  cell_annotation.Group = c("MS" = "#FF0000", "MSP2" = "#004C99"))
Cairo(file="Plot_AUC_of_HEATMAP_refined.png",type="png",units="in",bg="white",width=15,height=4,pointsize=3,dpi=300)
pheatmap(miniAssigMat,annotation_col = annotation_col,annotation_colors = ann_colors,
         border_color = "black", cellwidth = 1.5, cellheight = 20,
         cutree_rows = 3, cutree_cols = 3,
         legend = FALSE)
dev.off()
#Plot marker genes of SHH-A,SHH-B and SHH-C in MB cells.Repalce "Sox2" with "Top2a","Cdk1","Mki67","Eif3e","Eef1a1","Boc","Stmn2","Map1b","Tubb2b","Atoh1" one by one.(Figure 3a)
marker_genes <- c("Sox2")
Cairo(file="Plot_Sox2_facet.png",type="png",units="in",bg="white",width=4.9,height=3.75,pointsize=5,dpi=300)
plot_cell_clusters(HSMM_all,
                   markers = as.character(marker_genes),
                   show_group_id = T, cell_size = 0.75) +
  facet_wrap(~Group) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#------- End. Charpter 4.--------

#------- Charpter 5. Enrichment analysis of Cell Cycle Phase gene signature and mTOR signaling pathway in all the cells.-------
#Clean the enrivoment. Repeat Charpter 1.2.3 to reconstruct the HSMM_all object, skip charpter 4 and resume from the following row.
fData(HSMM_all)$gene_id <- row.names(fData(HSMM_all))
target_genes <- fData(HSMM_all)[,-2:-3]
target_genes <- unique(target_genes)
Mouse_gene_id <- target_genes$gene_short_name
library(biomaRt)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
class(human)
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")
human_gene_id <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                        values = Mouse_gene_id, mart = mouse,
                        attributesL = c("hgnc_symbol"), martL = human,
                        uniqueRows = T)
colnames(human_gene_id) <- c("gene_short_name","human_symbol")
target_genes <- merge(target_genes,human_gene_id,by="gene_short_name",all=FALSE)
#load scRNAseq gene matrix and rename the gene with human gene symbol
HSMM_all <- HSMM_all[target_genes$gene_id,]
df <- exprs(HSMM_all)
df <- as.matrix(df)
rownames(df) <- target_genes$human_symbol
#load geneset file
geneSets <- getGmt("MB_signature_mTOR_Cell_Cycle.gmt.txt")
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))
cells_rankings <- AUCell_buildRankings(df, nCores=12, plotStats=TRUE,verbose =TRUE)
#Calculating the AUC score of each gene set in every cells
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),verbose = TRUE,
                            nCores = 10)
#plot the histogram
Cairo(file="Plot_AUC_of_all_cells.png",type="png",units="in",bg="white",width=6,height=6,pointsize=3,dpi=300)
set.seed(123)
par(mfrow=c(4,4))
cells_assignment <- AUCell_exploreThresholds(cells_AUC, thrP = 0.01, nCores = 12, 
                                             smallestPopPercent = 0.25, plotHist=TRUE, assign=TRUE,
                                             nBreaks = 500) 
dev.off()
#--Generate AUC score matrix and combine them to pData(HSMM_MB)
AUC_matrix <- getAUC(cells_AUC)[1:25,1:8617]
AUC_matrix <- t(AUC_matrix)
AUC_matrix <- as.data.frame(AUC_matrix)
pData(HSMM_all) <- cbind(pData(HSMM_all),AUC_matrix)
write.table(pData(HSMM_all),"pData_all_with_Genesets.txt",quote=FALSE,sep='\t',row.names = TRUE,col.names = TRUE)
#change geneset Id in text files. Change SHHA(100g) to SHHA-100g
#------- End. Charpter 5.--------

#------- Charpter 6. Identify Cell Cycle Phase--------
#Repeat Charpter 1.2.3 to reconstruct HSMM_all object and skip charpter 4 and charpter 5, resume from the following row.
#A two step method for G0 cell identification.
#Install Seurat
install.packages('Seurat')
library(Seurat)
library(dplyr)
#Transform the mouse gene symbol to human gene symbol
fData(HSMM_all)$gene_id <- row.names(fData(HSMM_all))
target_genes <- fData(HSMM_all)[,-2:-3]
target_genes <- unique(target_genes)
Mouse_gene_id <- target_genes$gene_short_name
library(biomaRt)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
class(human)
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")
human_gene_id <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                        values = Mouse_gene_id, mart = mouse,
                        attributesL = c("hgnc_symbol"), martL = human,
                        uniqueRows = T)
colnames(human_gene_id) <- c("gene_short_name","human_symbol")
target_genes <- merge(target_genes,human_gene_id,by="gene_short_name",all=FALSE)
#load scRNAseq gene matrix and rename the gene with human gene symbol
HSMM_all <- HSMM_all[target_genes$gene_id,]
df <- exprs(HSMM_all)
df <- as.matrix(df)
rownames(df) <- target_genes$human_symbol
#Load cell cycle markers from Tirosh et al, 2015.
#Segregate list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#Creat Seurat Object
MB_Seurat <- CreateSeuratObject(counts = df)
MB_Seurat <- NormalizeData(MB_Seurat)
MB_Seurat <- FindVariableFeatures(MB_Seurat,selection.method = "vst")
MB_Seurat <- ScaleData(MB_Seurat,features = rownames(MB_Seurat))
#Run PCA in Seurat
MB_Seurat <- RunPCA(MB_Seurat,features = VariableFeatures(MB_Seurat),ndims.print = 1:4,nfeatures.print = 10)
#Plot Heatmap of Principle Components
DimHeatmap(MB_Seurat,dims = c(1,2,3))
#Assign MB cells into cell cycle stage
MB_Seurat <- CellCycleScoring(MB_Seurat,s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
#View cell cycle scores and phase assignments
pData(HSMM_all) <- read.delim("pData_all_with_Genesets.txt",header = TRUE, row.names = 1)
pData(HSMM_all) <- cbind(pData(HSMM_all),MB_Seurat@meta.data)
#Save pData(HSMM_all) as pData_all_with_Cell_Cycle
write.table(pData(HSMM_all),"pData_all_with_Cell_Cycle.txt",quote=FALSE,sep='\t',row.names = TRUE,col.names = TRUE)
#-------- End. Charpter 6.----------

#-------- Charpter 7.--------
#Reconstruction of trajectory and pseudotime in Sox2+ MB cells.
#Trajectory analysis
#---Monocle 2---
setwd("E:/share18/share1/Xin/analysis/Reproduce/")
#loading the data of MS and MSP2
library(Matrix)
matrix_dir = "E:/share18/share1/Xin/analysis/Reproduce/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "features.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
matrix <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(matrix) = barcode.names$V1
rownames(matrix) = feature.names$V1

setwd("E:/share18/share1/Xin/analysis/Reproduce/")
library(monocle)
library("ggthemes")
library("RColorBrewer")
library(VGAM)
library(Cairo)
library(dplyr)
HSMM_sample_sheet_all<-read.delim("pData_all_with_Cell_Cycle.txt",header=T,row.names = 1) 
HSMM_gene_annotation_all<-read.delim("features.txt",header=T,row.names = 1)
pd<-new("AnnotatedDataFrame",data = HSMM_sample_sheet_all)
fd<-new("AnnotatedDataFrame",data = HSMM_gene_annotation_all)
HSMM_all<-newCellDataSet(matrix,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.1, 
                         expressionFamily = negbinomial.size()) 
HSMM_all <- estimateSizeFactors(HSMM_all)
HSMM_all <- estimateDispersions(HSMM_all)
HSMM_all <- detectGenes(HSMM_all,min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(HSMM_all),num_cells_expressed >=10)) 
HSMM_all <- setOrderingFilter(HSMM_all,expressed_genes)
#pData(HSMM_all) <- read.delim("pData_all_with_Cell_Cycle.txt",header=TRUE,row.names = 1)
#Generate MB cell population
MB_cells <- row.names(subset(pData(HSMM_all),Cluster == 1 | Cluster == 2 | Cluster == 5 | Cluster == 8 | Cluster == 9
                             | Cluster == 10 | Cluster == 12))
HSMM_SHH_MB <- HSMM_all[,MB_cells]
#Generate Sox2 MB cell population
Sox2_id <- row.names(subset(fData(HSMM_SHH_MB),gene_short_name == "Sox2"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth,"Stem_cell",classify_func = function(x){x[Sox2_id,] >= 1})
HSMM_SHH_MB <- classifyCells(HSMM_SHH_MB,cth,0.1)
pData(HSMM_SHH_MB)$regroup <- dplyr::recode(pData(HSMM_SHH_MB)$Cluster,
                                            "1"="SHH_B",
                                            "5"="SHH_B",
                                            "10"="SHH_B",
                                            "2"="SHH_A",
                                            "9"="SHH_A",
                                            "8"="SHH_C",
                                            "12"="SHH_C")
Stem_cells <- row.names(subset(pData(HSMM_SHH_MB),CellType == "Stem_cell"))
HSMM_SHH_Stem_cells <- HSMM_SHH_MB[,Stem_cells]
clustering_DEG_genes <- differentialGeneTest(HSMM_SHH_Stem_cells[expressed_genes,],
                                             fullModelFormulaStr = '~ regroup', cores = 12, verbose = TRUE)
write.table(clustering_DEG_genes,"DE_Gene_based_on_Subgroups.txt",sep='\t',quote=FALSE,col.names = TRUE,row.names = TRUE)
clustering_DEG_genes <- read.delim("DE_Gene_based_on_Subgroups.txt",header=T,row.names = 1)
uplimits<-0.05*(roundDown(nrow(subset(clustering_DEG_genes,qval<0.05)),dig=-2))
Gene_index <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:uplimits]
HSMM_SHH_Stem_cells <- setOrderingFilter(HSMM_SHH_Stem_cells,Gene_index) 
HSMM_SHH_Stem_cells <- reduceDimension(HSMM_SHH_Stem_cells, max_components = 2, num_dim = 6, norm_method = c("log"), pseudo_expr = 1,
                                       reduction_method = 'DDRTree', verbose = T, residualModelFormulaStr = "~num_genes_expressed", cores = detectCores() - 2)
HSMM_SHH_Stem_cells <- orderCells(HSMM_SHH_Stem_cells)
pData(HSMM_SHH_Stem_cells)$Cluster <- as.factor(pData(HSMM_SHH_Stem_cells)$Cluster)
#Plot trajectory of Sox2 MBs according to MB molecular subgroup
Cairo(file="Plot_Trajectory_by_Subgroup_facet.png",type="png",units="in",bg="white",width=4.95,height=3.75,pointsize=5,dpi=300)
plot_cell_trajectory(HSMM_SHH_Stem_cells, 1, 2, color_by = "regroup", cell_size = 0.5,state_number_size = 0.1,show_branch_points = TRUE) + 
  facet_wrap(~Group) +
  #stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme(legend.text=element_text(size=3)) + 
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#Plot expression of Cdk1,Top2a,Mki67 in trajectory plot (Figure 3b,small part, cell cycle gene expression)
#Change gene_short_name to Mki67, Top2a and Cdk1. Adjust midpoint and limits according to the gene expression.
target<-rownames(fData(HSMM_SHH_Stem_cells)[which(fData(HSMM_SHH_Stem_cells)$gene_short_name=="Cdk1"),])
pData(HSMM_SHH_Stem_cells)$Cdk1<-exprs(HSMM_SHH_Stem_cells)[target,]
Cairo(file="Plot_Trajectory_Cdk1_facet.png",type="png",units="in",bg="white",width=4.95,height=3.75,pointsize=10,dpi=300)
plot_cell_trajectory(HSMM_SHH_Stem_cells, color_by = "Cdk1", cell_size = 1) + 
  scale_colour_gradient2(low="white",mid="Magenta2",high="Red",midpoint = 10,limits = c(0,20)) +
  facet_wrap(~Group) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

#Choose the root of trajectory for pseudotime calculation
Cdk1_id <- row.names(subset(fData(HSMM_SHH_Stem_cells),gene_short_name == "Cdk1"))
Top2a_id <- row.names(subset(fData(HSMM_SHH_Stem_cells),gene_short_name == "Top2a"))
Mki67_id <- row.names(subset(fData(HSMM_SHH_Stem_cells),gene_short_name == "Mki67"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Low-Cycling", classify_func = function(x) { x[Cdk1_id,] < 1 & x[Top2a_id,] < 1 & x[Mki67_id,] < 1})
cth <- addCellType(cth, "High-Cycling", classify_func = function(x) { x[Cdk1_id,] > 1 | x[Top2a_id,] >1 | x[Mki67_id,] > 1})
HSMM_SHH_Stem_cells <- classifyCells(HSMM_SHH_Stem_cells, cth, 0.1)
#Pick up the branch with most low-cycling cells as the root of trajectory
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$CellType)[,"Low-Cycling"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
#recaculate the pseudotime according to the newly picked up trajectory root.
HSMM_SHH_Stem_cells <- orderCells(HSMM_SHH_Stem_cells, root_state = GM_state(HSMM_SHH_Stem_cells))
#plot the pseudotime (figure 3d)
Average <- (max(pData(HSMM_SHH_Stem_cells)$Pseudotime)+min(pData(HSMM_SHH_Stem_cells)$Pseudotime))/2
Cairo(file="Plot_Trajectory_Pseudotime_facet.png",type="png",units="in",bg="white",width=4.95,height=3.95,pointsize=10,dpi=300)
plot_cell_trajectory(HSMM_SHH_Stem_cells, 1, 2, color_by = "Pseudotime", cell_size = 1,nrow = 2)+ 
  scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint = Average) +
  facet_wrap(~Group) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#Plot genes branched pseudotime of Top2a,Cdk1,Mki67
#Comfirm the pseudotime is consistent with quiescence/cell cylcing alteration.
Markers_Cell_Cycle <- rownames(fData(HSMM_SHH_Stem_cells)
                               [which(fData(HSMM_SHH_Stem_cells)$gene_short_name %in% 
                                        c("Top2a","Mki67","Cdk1")),])
col_vector <- c("#00D4F5","#18F500","#EC5700")
Cairo(file="Plot_genes_branches_pseudotime_Top2a_Mki67_Cdk1_gene.png",type="png",units="in",bg="white",width=6,height=2.5,pointsize=10,dpi=300)
plot_genes_branched_pseudotime(HSMM_SHH_Stem_cells[Markers_Cell_Cycle,],color_by = "State",
                               branch_point = 1,branch_labels = c("Branch 1","Branch 2"),
                               ncol = 3) + 
  scale_color_manual(values = col_vector) +
  #scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint = Average) +
  #facet_wrap(~Group) +
  #stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

#Identify DE genes between Branch1 and Branch2 divided by branch-point 1
BEAM_res <- BEAM(HSMM_SHH_Stem_cells[expressed_genes,], branch_point = 1,cores = 10)
Sig_BEAM_res <- subset(BEAM_res,qval < 0.1)
write.table(Sig_BEAM_res,"DE_genes_based_on_Branchpoint_1_2_branches.txt",sep='\t',quote=FALSE,col.names = TRUE,row.names = TRUE)
Sig_BEAM_res <- read.delim("DE_genes_based_on_Branchpoint_1_2_branches.txt",header=TRUE, row.names = 1)
Sig_Genes <- rownames(Sig_BEAM_res)
a<-plot_genes_branched_heatmap(HSMM_SHH_Stem_cells[Sig_Genes,],branch_point = 1,  
                               num_clusters = 12, cores = 10, use_gene_short_name = T,
                               show_rownames = TRUE,return_heatmap = TRUE)
Cairo(file="Plot_heatmap_BEAM_Gene.png",
      type="png",units="in",bg="white",width=5,height=8,pointsize=10,dpi=300)
a$ph_res
dev.off()
#checked the heatmap of plot_heatmap_BEAM_Gene.png and refine the identifical genes that are enriched in pre-branch, 
#branch-1 and branch-2. Generate the refined and tidy gene list file, DE_genes_based_on_Branchpoint_1_2_branches_target.txt.
#Reproduce the plot_genes_branched_heatmap according to the refined list.(Figure 3c)
Sig_BEAM_res_target <- read.delim("DE_genes_based_on_Branchpoint_1_2_branches_target.txt",header=TRUE, row.names = 1)
Sig_Genes_target <- rownames(Sig_BEAM_res_target)
a<-plot_genes_branched_heatmap(HSMM_SHH_Stem_cells[Sig_Genes_target,],branch_point = 1,  
                               num_clusters = 3, cores = 10, use_gene_short_name = T,
                               show_rownames = TRUE,return_heatmap = TRUE)
Cairo(file="Plot_heatmap_BEAM_Gene_target.png",
      type="png",units="in",bg="white",width=5,height=8,pointsize=10,dpi=300)
a$ph_res
dev.off()
DE_genes<-a$annotation_row
write.table(DE_genes,"clusters_of_DE_genes_target.txt",col.names = TRUE,row.names = TRUE,sep="\t",quote=FALSE)
#Rename the clusters with the "Branch1-up", "Branch2-up" and "Pre-branch-up" in clusters_of_DE_genes_target.txt file with the colname as "Genes" and "Feature".

#Enrichment analysis of Branch1 and Branch2. Reproduce the results of enriched GO terms in Figure 3c.
DE_genes <- read.delim("clusters_of_DE_genes_target.txt",header=TRUE)
#Enrichment analysis of DE genes
Branch1_genes <- subset(DE_genes,Feature == "Branch1-up")
Branch2_genes <- subset(DE_genes,Feature == "Branch2-up")
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(biomaRt)
#------Branch1 genes----
SYMBOL<-Branch1_genes$Genes
human=useMart("ensembl",dataset="hsapiens_gene_ensembl")
class(human)
mouse=useMart("ensembl",dataset="mmusculus_gene_ensembl")
class(mouse)
SYMBOLmouse<-getLDS(attributes=c("mgi_symbol"),filters="mgi_symbol",
                    values=SYMBOL,mart=mouse,
                    attributesL=c("hgnc_symbol","chromosome_name","start_position","end_position"),
                    martL=human,
                    uniqueRows = T)
symbol<-SYMBOLmouse$HGNC.symbol
eg1<-bitr(symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
gmtfileGO<-system.file("extdata","c5.bp.v7.0.symbols.gmt",package="clusterProfiler")
c5<-read.gmt(gmtfileGO)
#GO enrichment ananlysis
c5<-enricher(eg1$SYMBOL,TERM2GENE = c5)
write.table(c5,"enrichment_c5_branch1.txt",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#Dot plot
Cairo(file="Results_C5_branch1_enrichment_all.png",type="png",units="in",bg="white",width=20,height=60,pointsize=12,dpi=300)
dotplot(c5,color="p.adjust",showCategory=355,font.size=15)
dev.off()
Cairo(file="Results_C5_branch1_enrichment_concise.png",type="png",units="in",bg="white",width=15,height=5,pointsize=12,dpi=300)
dotplot(c5,color="p.adjust",showCategory=15,font.size=15)
dev.off()
#------Branch2 genes----
SYMBOL<-Branch2_genes$Genes
human=useMart("ensembl",dataset="hsapiens_gene_ensembl")
class(human)
mouse=useMart("ensembl",dataset="mmusculus_gene_ensembl")
class(mouse)
SYMBOLmouse<-getLDS(attributes=c("mgi_symbol"),filters="mgi_symbol",
                    values=SYMBOL,mart=mouse,
                    attributesL=c("hgnc_symbol","chromosome_name","start_position","end_position"),
                    martL=human,
                    uniqueRows = T)
symbol<-SYMBOLmouse$HGNC.symbol
eg1<-bitr(symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
gmtfileGO<-system.file("extdata","c5.bp.v7.0.symbols.gmt",package="clusterProfiler")
c5<-read.gmt(gmtfileGO)
#GO enrichment ananlysis
c5<-enricher(eg1$SYMBOL,TERM2GENE = c5)
write.table(c5,"enrichment_c5_branch2.txt",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#Dot plot
Cairo(file="Results_C5_branch2_enrichment_all.png",type="png",units="in",bg="white",width=15,height=30,pointsize=12,dpi=300)
dotplot(c5,color="p.adjust",showCategory=355,font.size=15)
dev.off()
Cairo(file="Results_C5_branch2_enrichment_concise.png",type="png",units="in",bg="white",width=12.6,height=5,pointsize=12,dpi=300)
dotplot(c5,color="p.adjust",showCategory=15,font.size=15)
dev.off()

#Plot marker gene heatmap for Cell-Cycle-DNA-replication (branch1),Cell-Cycle-Spindle-Organization (Branch2)
#This example shows how to plot individualized BEAM heatmap based on targting pathways.
marker <- rownames(fData(HSMM_SHH_Stem_cells)
                   [which(fData(HSMM_SHH_Stem_cells)$gene_short_name %in% 
                            c("Rfc1","Rad51","Prim2","Tipin","Rpa1","Fbxo5","Cdc45","Rfc2","E2f8",
                              "Atad5","Fen1","Brca2","Lig1","Prim1","Slbp","E2f7","Pcna","Rfc4",
                              "Cdt1","Rfc5","Rfc3","Pola1","Smc3","Rpa2","Gmnn","Ccnb1","Tpx2","Aspm","Cenpe","Cdc20","Knstrn","Psrc1")),])
col_annotation <- pData(HSMM_SHH_Stem_cells)[,c("Pseudotime","regroup","CellType")]
col_annotation <- as.data.frame(col_annotation)
a<-plot_genes_branched_heatmap(HSMM_SHH_Stem_cells[marker,], branch_point = 1,
                               num_clusters = 2, cores = 10, use_gene_short_name = T,
                               show_rownames = TRUE,return_heatmap = TRUE)
Cairo(file="Plot_heatmap_BEAM_Marker_gene_DNA_replication_and_Spindel_organization.png",type="png",units="in",bg="white",width=7,height=5,pointsize=10,dpi=300)
a$ph_res
dev.off()
#Plot marker genes for each enriched GO terms.
Markers_Enriched <- rownames(fData(HSMM_SHH_Stem_cells)
                             [which(fData(HSMM_SHH_Stem_cells)$gene_short_name %in% 
                                      c("Cdc20","Cenpe","Ccnb1","Cenpf","Cenpa","Cenpe","Map1b","Ptn","Lgals1","Mcm2","Msh6","Msh2",
                                        "Ezh2","Nsd2","Hist1h1c","Pcna","Rif1","Brca2")),])
col_vector <- c("#00D4F5","#18F500","#EC5700")
Cairo(file="Plot_marker_gene_branches_in_GO_terms.png",type="png",units="in",bg="white",width=7,height=7,pointsize=10,dpi=300)
plot_genes_branched_pseudotime(HSMM_SHH_Stem_cells[Markers_Enriched,],color_by = "State",
                               branch_point = 1,branch_labels = c("Branch 1","Branch 2"),
                               ncol = 4) + 
  scale_color_manual(values = col_vector) +
  #scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint = Average) +
  #facet_wrap(~Group) +
  #stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#-------- End. Charpter 8 ---------

#-------- Chapter 9. Comparing Piezo2 wt and Piezo KO in Peudotime ------
#Comparing the pseudotime and disctance on trajectory from quiescence root between Piezo2 Wt and Piezo2 KO.
#Regenerate the pseudotime trajecotry in both genotypes (Figure 3d)
HSMM_SHH_Stem_cells <- orderCells(HSMM_SHH_Stem_cells, root_state = GM_state(HSMM_SHH_Stem_cells))
Average <- (max(pData(HSMM_SHH_Stem_cells)$Pseudotime)+min(pData(HSMM_SHH_Stem_cells)$Pseudotime))/2
Cairo(file="Plot_Trajectory_Pseudotime_facet_1.png",type="png",units="in",bg="white",width=6,height=4.0,pointsize=10,dpi=300)
plot_cell_trajectory(HSMM_SHH_Stem_cells, 1, 2, color_by = "Pseudotime", cell_size = 1,nrow = 2)+ 
  scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint = Average) +
  facet_wrap(~Group) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#Ridge plot of pseudotime based on each branch.
library(ggplot2)
library(ggridges)
theme_set(theme_ridges())
pData(HSMM_SHH_Stem_cells)$Group_Num <- as.numeric(dplyr::recode(pData(HSMM_SHH_Stem_cells)$Group,
                                                                 "MSP2"="2",
                                                                 "MS"="1"))
pre_branch <- subset(pData(HSMM_SHH_Stem_cells),State=="3")
branch_1 <- subset(pData(HSMM_SHH_Stem_cells),State=="1")
branch_2 <- subset(pData(HSMM_SHH_Stem_cells),State=="2")
#Ridge plot of pre_branch,branch_1,branch_2. 
#Change pre_branch dataframe to branch_1 and branch_2 to regenerate Extended Data Fig 7.
#For pre-branch, Figure 3e uppper will be regenerated.
Cairo(file="Branch_Distribution_pre_branch.png",type="png",units="in",bg="white",width=6,height=2,pointsize=10,dpi=300)
ggplot(pre_branch,aes(x=Pseudotime,y=Group)) +
  geom_density_ridges_gradient(aes(fill=..x..),scale=2.0,size = 0.2, rel_min_height = 0.0001, gradient_lwd = 1.0) +
  scale_fill_gradientn(
    colours = c("orange", "magenta2", "red"),
    name = "Pseudotime"
  ) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_discrete(expand = c(0.5, 1)) +
  scale_fill_viridis_c(name="Pseudotime",option="C") +
  expand_limits(x=c(0,12)) +
  labs(title = 'Density of Pseudotime in Pre-branch',
       subtitle = 'Pseudotime of Cells in Each Genotype', 
       x = "Pseudotime") +
  theme_ridges(font_size = 10, grid = TRUE, line_size = 0.1, center_axis_labels=TRUE) +
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'))
dev.off()
#Violin plot of Pseudotime comparison in pre-branch,branch-1 and branch-2
#Change pre_branch dataframe to branch_1 and branch_2 to regenerate Extended Data Fig 7.
#For pre-branch, whole branch pseudotime will be compared.
Cairo(file="Pseudotime_Violin_pre_branch.png",type="png",units="in",bg="white",width=2.5,height=2.5,pointsize=10,dpi=300)
ggplot(pre_branch,aes(x=Group,y=Pseudotime,fill=Group)) +
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Pre-branch",x="Genotype",y="Pseudotime") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
#statistics
a<-aov(Pseudotime~Group,data=pre_branch)
b<-wilcox.test(Pseudotime~Group,data=pre_branch)
summary(a)
b
#Statistics of Pseudotime in over 0.25 course of pre-branch.(Figure 3d lower)
quantile_pseudotime <- quantile(pre_branch$Pseudotime,0.25)
pre_branch_pseudotime_1 <- subset(pre_branch,Pseudotime > quantile_pseudotime)
Cairo(file="Pseudotime_Violin_pre-branch-pseudotime-025.png",type="png",units="in",bg="white",width=2.5,height=2.5,pointsize=10,dpi=300)
ggplot(pre_branch_pseudotime_1,aes(x=Group,y=Pseudotime,fill=Group)) +
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Pre-branch",x="Genotype",y="Pseudotime") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
a<-aov(Pseudotime~Group,data=pre_branch_pseudotime_1)
b<-wilcox.test(Pseudotime~Group,data=pre_branch_pseudotime_1)
summary(a)
b
#Plot complicate tree structure to generate the trajectory model. (Figure 3f)
Cairo(file="Plot_Tree_by_Pseudotime.png",type="png",units="in",bg="white",width=3,height=5,pointsize=10,dpi=300)
plot_complex_cell_trajectory(HSMM_SHH_Stem_cells, 1, 2, color_by = "Pseudotime", cell_size = 2,cell_link_size = 0.5)+ 
  #scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint =) +
  #facet_wrap(~Group) +
  #stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "white", fill=NA, size=1))
dev.off()
#------ End. Charpter 9.--------

#------ Charpter 10. ---------
#Cell cycle analysis for Sox2+ MB cells. Identify cell cycle phase and 
#Recreat cell cycle gene signature AUCell enrichment score matrix for clustering by gitools.
cc_matrix <- pData(HSMM_SHH_Stem_cells)[,14:20]
a1<-(cc_matrix$M..65g.-0.175)/sd(cc_matrix$M..65g.)
b1<-(cc_matrix$G1..31g.-0.150)/sd(cc_matrix$G1..31g.)
c1<-(cc_matrix$S..72g.-0.098)/sd(cc_matrix$S..72g.)
d1<-(cc_matrix$G2..52g.-0.15)/sd(cc_matrix$G2..52g.)
norm_cc_matrix <- cbind(a1,b1,c1,d1)
rownames(norm_cc_matrix) <- rownames(cc_matrix)
colnames(norm_cc_matrix) <- c("M","G1","S","G2")
norm_cc_matrix <- as.matrix(norm_cc_matrix)
norm_cc_matrix <- t(norm_cc_matrix)
Cairo(file="Plot_heatmap_cell_cycle.png",type="png",units="in",bg="white",width=7,height=4,pointsize=10,dpi=300)
pheatmap(norm_cc_matrix)
dev.off()
write.table(norm_cc_matrix,"cell_cylce_normalized.txt",quote=FALSE,sep='\t',row.names = TRUE,col.names = TRUE)
#The pheatmap of normalized cc-matrix would generate a heatmap based on AUCell enrichment score of Cell cycle gene signature. 
#A group of cells with low enrichment score could be observed. To identify G0 cells, we generated the heatmap of cell cycle normalized matrix heatmap in gitools. 
#Cells in seurate identified G1 phase were clustered the into two clusters in gitools by K-means.
#For cells in Seurat identified G1, cells were reclustered into two clusters and those with low G1,S,G2,M AUCell enrichment score were identified as G0 cells.
#Save the pData(HSMM_SHH_stem_cells) and merge the clustering information generated from gitools in excel.
#Merge the cell_cycling_clustering and pData(HSMM_SHH_Stem_cells) in Excel to generate pData_Sox2_positive_cells.txt file.
#We provide pData_Sox2_positive_cells.txt here to reproduce the results.
pData(HSMM_SHH_Stem_cells) <- read.delim("pData_Sox2_positive_cells.txt",header=TRUE,row.names = 1)
#Plot non-cycling cells in trajectory
Cairo(file="Plot_Trajectory_by_Cell_Non_Cycling_gitools_facet.png",type="png",units="in",bg="white",width=4.95,height=3.75,pointsize=5,dpi=300)
plot_cell_trajectory(HSMM_SHH_Stem_cells, 1, 2, color_by = "Cell_cycle_final", cell_size = 0.5,state_number_size = 0.1,show_branch_points = TRUE) + 
  facet_wrap(~Group) +
  #scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint = 0.75) +
  #stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  scale_color_manual(values = c("white","Green")) +
  theme(legend.text=element_text(size=3)) + 
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#Plot cycling cells in trajectory
Cairo(file="Plot_Trajectory_by_Cell_Cycling_gitools_facet.png",type="png",units="in",bg="white",width=4.95,height=3.75,pointsize=5,dpi=300)
plot_cell_trajectory(HSMM_SHH_Stem_cells, 1, 2, color_by = "Cell_cycle_final", cell_size = 0.5,state_number_size = 0.1,show_branch_points = TRUE) + 
  facet_wrap(~Group) +
  #scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint = 0.75) +
  #stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  scale_color_manual(values = c("red","white")) +
  theme(legend.text=element_text(size=3)) + 
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#Plot cycling cells in trajectory. Reproduce figure 3b left.
Cairo(file="Plot_Trajectory_by_Cell_phase_gitools_facet.png",type="png",units="in",bg="white",width=4.95,height=3.75,pointsize=5,dpi=300)
plot_cell_trajectory(HSMM_SHH_Stem_cells, 1, 2, color_by = "Cell_cycle_phase", cell_size = 0.5,state_number_size = 0.1,show_branch_points = TRUE) + 
  facet_wrap(~Group) +
  #scale_colour_gradient2(low="yellow",mid="Magenta2",high="red",midpoint = 0.75) +
  #stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  scale_color_manual(values = c("green","blue","red","orange")) +
  theme(legend.text=element_text(size=3)) + 
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
#------- End. Charpter 10.---------

#------- Charpter 11. Combined analysis with Pseudotime, cell cycle phase, cell cycle signature in pre-branch.
#G1.31g in pre-branch G0 Loess cross validation
#Loess cross validation
Sox2MB <-read.delim("pData_Sox2_positive_cells.txt",header=TRUE,row.names = 1)
spans <- c(0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.4,0.35,0.3,0.25,0.2)
k <- length(spans)
viewports <- lapply(1:k,function(i)
  grid::viewport(width=1/k,height=1,x=(i-1/2)/k,y=1/2))
names(viewports) <- spans
G0 <- subset(Sox2MB,Cell_cycle_phase == "G0")
g <- ggplot(data=subset(G0,State == 3),aes(x=Pseudotime,y=G1.31g,group=Group,col=Group)) +
  coord_cartesian(ylim=c(0.03,0.15))
#Plot Loess Smoothing Fitting (all spans)
Cairo(file="Plot_G1.31g_loess.png",type="png",units="in",bg="white",width=60,height=3,pointsize=10,dpi=300)
for (i in 1:k){
  print(g + 
          stat_smooth(method = 'loess',span = spans[i],level = 0.90) +
          scale_color_manual(values = c("red","yellow")) +
          labs(title=paste("G0 Sox2 cells at pre-branch"," Span=",spans[i]),x="Pseudotime",y="G1.31g") +
          theme_classic(),
        vp=viewports[[i]])
}
dev.off()
#Plot Loess Smoothing Fitting (Figure 4c)
Cairo(file="Plot_G1.31g_G0_cells_prebranch.png",type="png",units="in",bg="white",width=5,height=3.75,pointsize=10,dpi=300)
ggplot(data=subset(G0,State == 3),aes(x=Pseudotime,y=G1.31g,group=Group,col=Group)) +
  #geom_point() +
  stat_smooth(method = 'loess') +
  scale_color_manual(values = c("red","yellow")) +
  #facet_wrap(~State) +
  labs(title="G0 Sox2 cells at pre-branch",x="Pseudotime",y="G1.31g") +
  theme_classic()
dev.off()

#G2M.score in pre-branch G2M Loess cross validation
G2M <- subset(Sox2MB,Cell_cycle_phase == "G2M")
#Plot Loess Smoothing Fitting (Figure 4c)
Cairo(file="Plot_G2M.Score_G2M_cells_prebranch.png",type="png",units="in",bg="white",width=5,height=3.75,pointsize=10,dpi=300)
ggplot(data=subset(G2M,State == 3),aes(x=Pseudotime,y=G2M.Score,group=Group,col=Group)) +
  #geom_point() +
  stat_smooth(method = 'loess') +
  scale_color_manual(values = c("red","yellow")) +
  #facet_wrap(~State) +
  labs(title="G2M Sox2 cells at pre-branch",x="Pseudotime",y="G1.31g") +
  theme_classic()
dev.off()

#------- End. Charpter 11.---------
#------- Charpter 12. Combined analysis with Pseutotime, cell cycle phase and Pi3k-Akt-mTOR hallmark gene signature in pre-branch.
#Akt_UP_MTOR_DN.V1.DN.184.g in pre-branch G0 Loess cross validation
#Loess cross validation
spans <- c(0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.4,0.35,0.3,0.25,0.2)
k <- length(spans)
viewports <- lapply(1:k,function(i)
  grid::viewport(width=1/k,height=1,x=(i-1/2)/k,y=1/2))
names(viewports) <- spans
#Plot Loess Smoothing Fitting (all spans)
Sox2MB<-read.delim("pData_Sox2_positive_cells.txt",header=TRUE,row.names = 1)
G0 <- subset(Sox2MB,Cell_cycle_phase == "G0")
g <- ggplot(data=subset(G0,State == 3),aes(x=Pseudotime,y=AKT_UP_MTOR_DN.V1_DN.184g,group=Group,col=Group)) +
  coord_cartesian(ylim=c(0,0.03))
Cairo(file="Plot_AKT_UP_MTOR_DN.V1_DN.184g_loess.png",type="png",units="in",bg="white",width=60,height=3,pointsize=10,dpi=300)
for (i in 1:k){
  print(g + 
          stat_smooth(method = 'loess',span = spans[i],level = 0.90) +
          scale_color_manual(values = c("red","yellow")) +
          labs(title=paste("G1 cells in Pre-branch"," Span=",spans[i]),x="Pseudotime",y="AKT_UP_MTOR_DN.V1_DN.184g") +
          theme_classic(),
        vp=viewports[[i]])
}
dev.off()
#Plot Loess Smoothing Fitting (span = 0.75, Figure 4c)
Cairo(file="Plot_AKT_UP_MTOR_DN.V1_DN.184g_G0_cells_prebranch.png",type="png",units="in",bg="white",width=5,height=3.75,pointsize=10,dpi=300)
ggplot(data=subset(G0,State == 3),aes(x=Pseudotime,y=AKT_UP_MTOR_DN.V1_DN.184g,group=Group,col=Group)) +
  #geom_point() +
  stat_smooth(method = 'loess') +
  scale_color_manual(values = c("red","yellow")) +
  #facet_wrap(~State) +
  labs(title="G0 cells in Pre-branch",x="Pseudotime",y="AKT_UP_MTOR_DN.V1_DN.184g") +
  theme_classic()
dev.off()






