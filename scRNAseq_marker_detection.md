# Find marker genes and DE genes (across treatments)
This workflow has code for detecting marker genes (for each seurat cluster) and DE genes across treatments. 

> Prerequisite: Processed Seurat object stored in RDS format.

```
library(Seurat)
library(Seurat)
library(dplyr)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggthemes)
library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(glmGamPoi)
library(clustree)
library(ggraph)
library(SingleR)
library(presto) 
library(msigdbr) 
library(fgsea) 

library(future)
options(future.globals.maxSize = 8000 * 1024^2)

```
## Read data (Objects saved after QC and singleR cell-type identificiations)
```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )

Idents(sample) <- "orig.ident"
```
## Find marker genes for each cluster and store as CSV
This code will determine top marker genes for each seurat cluster as compared to all other clusters.
```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )
DefaultAssay(sample) <- "RNA"

markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

name = "Tumor"
outname = paste0(name,"_markers.rds")
saveRDS(markers, file = outname)

outname = paste0(name,"_markers.csv")
write.csv(markers, file = outname, row.names = F)

markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
```

## Find DE genes across treatments (irrespetive of clusters)
We assume there are two conditions `tumor-1` and `tumor-2` and determine DE genes across treatments
```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )

Idents(sample) <- "orig.ident"
DefaultAssay(sample) <- "RNA"


x = FindMarkers(sample, ident.1 = "tumor-1", ident.2 = "tumor-2", verbose = FALSE, only.pos = T)
x = rownames_to_column(x, "Gene_ID")
name = "tumor-1_vs_tumor-2"
outname = paste0(name,"_markers.txt")
write.table(x, file = outname, sep = "\t", quote = F, row.names = F)

```

## Find DE genes across cell-types and each treatment
Assuming you have singleR cell-types stored in `sample@meta.data$celltype` and treatment stored in `sample$orig.ident` after SCTransform.  
Also, we restrict this analysis when the number of cells are **greater than 20**.  

```
name = "Tumor"
outname = paste0(name,"_singleR_ImmGenData.rds")
sample = readRDS(file = outname )

DefaultAssay(sample) <- "integrated"


sample$celltype.treat <- paste(sample$celltype, sample$orig.ident,  sep = "_")
my_order = sort(unique(sample$celltype.treat))
Idents(sample) <- factor(sample$celltype.treat, levels = my_order)

setwd("/depot/pccr/data/Low_lab_projects/hr03128_single-cell-RNAseq-study_2021/2_Seurat/Integrate/Tumor/Integrate_res0.6/ImmGenData_DE_Analysis")

treatment = sort(unique(sample@meta.data$orig.ident))
N1 = c(1,1,1,2,2,3)
N2 = c(2,3,4,3,4,4)
cell_type =  sort(unique(sample@meta.data$celltype))

cell_counts = table(sample@meta.data$celltype.treat)

for (cell in sort(unique(sample@meta.data$celltype))) {
    

      for(i in 1:6) {
      ID1 = paste(cell, treatment[N1[i]], sep = "_")
      ID2 = paste(cell, treatment[N2[i]], sep = "_")
      comp_name = paste0(ID1, "_VS_", ID2)
      
      print(comp_name)
      
      if(ID1 %in% names(cell_counts) == FALSE)  {next}
      if(ID2 %in% names(cell_counts) == FALSE)  {next}
          
      if(cell_counts[ID1] < 20) {next}
      if(cell_counts[ID2] < 20 ) {next}
      
      x = FindMarkers(sample, ident.1 = ID1,  ident.2 = ID2, verbose = FALSE, only.pos = T)
      x = rownames_to_column(x, "Gene_ID")
      outname = paste0(comp_name,".txt")
      write.table(x, file = outname, sep = "\t", quote = F, row.names = F)
      
      }

  
}

```

As hypothetical example, if you have two assigned cell-types (Macrophages and Fibroblasts) and four treatments (1-tumor, 2-tumor, 3-tumor and 4-tumor), you will get following comparisons:  

1. Fibroblasts_1-tumor_VS_Fibroblasts_2-tumor.txt
2. Fibroblasts_1-tumor_VS_Fibroblasts_3-tumor.txt
3. Fibroblasts_1-tumor_VS_Fibroblasts_4-tumor.txt
4. Fibroblasts_2-tumor_VS_Fibroblasts_3-tumor.txt
5. Fibroblasts_2-tumor_VS_Fibroblasts_4-tumor.txt
6. Fibroblasts_3-tumor_VS_Fibroblasts_4-tumor.txt
7. Macrophages_1-tumor_VS_Macrophages_2-tumor.txt
8. Macrophages_1-tumor_VS_Macrophages_3-tumor.txt
9. Macrophages_1-tumor_VS_Macrophages_4-tumor.txt
10. Macrophages_2-tumor_VS_Macrophages_3-tumor.txt
11. Macrophages_2-tumor_VS_Macrophages_4-tumor.txt
12. Macrophages_3-tumor_VS_Macrophages_4-tumor.txt


## Find DE genes across clusters and within each treatment

```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )
Idents(sample) <- "orig.ident"


sample$cluster.treat <- paste(sample$seurat_clusters, sample$orig.ident,  sep = "_")
my_order = sort(unique(sample$cluster.treat))
Idents(sample) <- factor(sample$cluster.treat, levels = my_order)

treatment = sort(unique(sample@meta.data$orig.ident))

N1 = c(1,1,1,2,2,3)
N2 = c(2,3,4,3,4,4)

cluster =  sort(unique(sample@meta.data$seurat_clusters))

cell_counts = table(sample@meta.data$cluster.treat)

for (cell in sort(unique(sample@meta.data$seurat_clusters))) {
    

      for(i in 1:6) {
      ID1 = paste(cell, treatment[N1[i]], sep = "_")
      ID2 = paste(cell, treatment[N2[i]], sep = "_")
      comp_name = paste0(ID1, "_VS_", ID2)
      
      print(comp_name)
      
      if(ID1 %in% names(cell_counts) == FALSE)  {next}
      if(ID2 %in% names(cell_counts) == FALSE)  {next}
      #
      if(cell_counts[ID1] < 20) {next}
      if(cell_counts[ID2] < 20 ) {next}
      #
      x = FindMarkers(sample, ident.1 = ID1,  ident.2 = ID2, verbose = FALSE, only.pos = T)
      x = rownames_to_column(x, "Gene_ID")
      outname = paste0(comp_name,".txt")
      write.table(x, file = outname, sep = "\t", quote = F, row.names = F)

      }

  
}

```
As hypothetical example, if you have two seurat assigned clusters (0 and 1) and and four treatments (1-tumor, 2-tumor, 3-tumor and 4-tumor), you will get following comparisons:

1. 0_1-tumor_VS_0_2-tumor.txt
2. 0_1-tumor_VS_0_3-tumor.txt
3. 0_1-tumor_VS_0_4-tumor.txt
4. 0_2-tumor_VS_0_3-tumor.txt
5. 0_2-tumor_VS_0_4-tumor.txt
6. 0_3-tumor_VS_0_4-tumor.txt
7. 1_1-tumor_VS_1_2-tumor.txt
8. 1_1-tumor_VS_1_3-tumor.txt
9. 1_1-tumor_VS_1_4-tumor.txt
10. 1_2-tumor_VS_1_3-tumor.txt
11. 1_2-tumor_VS_1_4-tumor.txt
12. 1_3-tumor_VS_1_4-tumor.txt


