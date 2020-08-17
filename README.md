#Browsing Sox2+ BTICs in in Math1-Cre;SmoM2 or Math1-Cre;SmoM2:Piezo2-fl/fl SHH MB mice.
  This code repository is supportive for single cell RNAseq analysis for Piezo2 project from Xi Huang's Lab, Sickkids, CA. The raw-data of this project could be downloaded from GSEXXXX (link). This code could be used to reproduce the results that are published on (link). 
 The pipelines we used for single cell RNAseq analysis are majorly based on Monocle2, Monocle3-alpha, Seurat and AUCells from SCENIC. Except Monocle2, all the pinelines are excuted on Linux system. Monocle2 is majorly used for pseudotime and trajectory analysis. Monocle3-alpha is used for UMAP dimension reduction and cell type identification. Seurat is used for Cell Cycle identification. SCENIC is used for single cell gene-set enrichment analysis.
     The links of these major pinelines are listed as below:
     Monocle2, http://cole-trapnell-lab.github.io/monocle-release/docs/
     Monocle3-alpha, http://cole-trapnell-lab.github.io/monocle-release/monocle3/
     Seurat, https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
     AUCell from SCENIC, https://scenic.aertslab.org/tutorials/AUCell.html
Essential files to regenerate the results are also provided as below:
     barcodes.txt, re-annotated barcodes according to Monocle2 and Monocle3-alpha.
     features.txt, re-annotated features according to Monocle2 and Monocle3-alpha.
     MB_signature.gmt.txt, gene signature of three MB transcriptomic prgrommes.
     MB_signature_mTOR_Cell_Cycle.gmt.txt, gene signature for enrichment analysis and cell cylce phase identification.
     pData_Sox2_positive_cells.txt, pre-prepared results and information to regenrate the results of Sox2+ MB cells. Most information of this file will be created during the pipeline. For easy to regenerate the results, we provided the file directly here.
    markers_of_each_sub_cluster_17_clusters.txt, pre-prepared feature genes of each cluster. For principalGraphTest() takes more than 48 hours to generate the results, we uploaded the result file for easy.
    DE_genes_based_on_Branchpoint_1_2_branches_target, refined DE genes generated from BEAM.
    clusters_of_DE_genes_target, cluster information of branch dependent gene for enrichment analysis.

    To regenrate the results, we need to install some certain dependent libraries and packages.The packages that need to be installed are listed below:
   "Matrix","Monocle2","Monocle3-alpha","ggthemes","RColorBrewer","VGAM","Cairo","dplyr",
   "AUCell","GSEABase","biomaRt","NMF","pheatmap","Seurat","clusterProfiler","org.Mm.eg.db",
   "ggplot2","ggridges"

   The code used for this project is attached as Reproductive_Code_Final.R.