# Integrate scRNAseq samples with scTransform
```
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

library(future)
options(future.globals.maxSize = 8000 * 1024^2)

```
## Read data (These are objects saved after QC or DoubletFInder)
```
# Read data after QC

tumor1   = readRDS(file = "tumor1_QC.rds")
tumor2   = readRDS(file = "tumor2_QC.rds")
tumor3   = readRDS(file = "tumor3_QC.rds")
tumor4   = readRDS(file = "tumor4_QC.rds")

sample_list = list(tumor1, tumor2, tumor3, tumor4)
names(sample_list) = c("tumor1", "tumor2", "tumor3", "tumor4")
sample_list
```

## Process through Seurat defults
1. SelectIntegrationFeatures (top 3000 or could be custom number)
2. PrepSCTIntegration
3. FindIntegrationAnchors 
4. IntegrateData
5. RunPCA, RunUMAP, FindNeighbors

> For large number of samples implementation of `library(future)` is essential to incarporate parallel processing.

```
sample_list <- lapply(X = sample_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)
sample_list <- PrepSCTIntegration(object.list = sample_list, anchor.features = features)


anchors <- FindIntegrationAnchors(object.list = sample_list , normalization.method = "SCT",
    anchor.features = features)

combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:50)
combined.sct <- FindNeighbors(combined.sct, dims = 1:50)

```

## Optimal resolution detection through Clustree
Clustering is a core tool for analysing single-cell RNA-sequencing (scRNA-seq) datasets. The clustering is primarily controlled by two parameters, number of principle components and then resolution.
A clustering tree visualises the relationships between at a range of resolutions.
```
fortree = readRDS(file = "tumor_SCT.rds")

res.vec <- seq(0,1, by = 0.1)

for (r in res.vec){
  fortree <- FindClusters(
    fortree , 
    resolution = r)
}

outname = paste0(name,"_clustree.png")

png(file=outname  ,units = "in", width = 12, height = 10, bg = "white", pointsize = 1, res=300)
clustree(fortree, assay="SCT",prefix="SCT_snn_res.")
dev.off()

```
![**Figure A**](/images/ISP_0.png)  

## Cluster visualization at various resolutions (PDF output)
```
res = seq(0.1, 1, by = 0.1)

name = "Tumor"
outname = paste0(name,"_SCT.pdf")

pdf(outname,  width = 10, height = 8)

for (r in as.numeric(res)) {
    
    my_cluster <- FindClusters(combined.sct, resolution = r)

    title = paste0(" resolution = ", r)
    print(DimPlot(my_cluster, reduction = "umap", label = TRUE, label.size = 6) + ggtitle(title))

}

dev.off()

```
## Elbowplot:

Identifying the true dimensionality of a dataset and the most significant PC can be challenging/uncertain. Elbowplot method generates a ranking of principle components based on the percentage of variance explained by each one. 
In this example, we can observe an **elbow** (i.e. beginning of the straight line) somewhere between PC 40-50, suggesting that the majority of true signal is captured in the first 50 PCs.

```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )

p = 20
r = 0.3

title = paste0("nPCs = ", p,", resolution = ", r) 

ElbowPlot(sample, ndims = 50)
```
![**Figure B**](/images/ISP_1.png)  

## Show dimension plot (dimPlot)
```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )

p = 20
r = 0.3

title = paste0("nPCs = ", p,", resolution = ", r) 

DimPlot(sample , reduction = "umap", label = TRUE, label.size = 6) + ggtitle(title)
```

![**Figure C**](/images/scRNAseq_integrate_1.png.png)  


## SingleR cell type assignment
[SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html) is an automatic annotation method for single-cell RNA sequencing (scRNAseq) data. 
Given a reference dataset of samples (single-cell or bulk) with known labels, it labels new cells from a test 
dataset based on similarity to the reference.  

Cell type assignment is performed using cell type reference datasets provided in [celldex package](https://bioconductor.org/packages/3.13/data/experiment/html/celldex.html).
Most of these references are derived from bulk RNA-seq or microarray data of cell populations that (hopefully) consist of a pure cell type after sorting and/or culturing. 

## Current data is from mouse genome

### Assignment by Mouse RNA-seq
This reference consists of a collection of mouse bulk RNA-seq data sets downloaded from the gene expression omnibus. A variety of cell types are available, again mostly from blood but also covering several other tissues.  

```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )

ref <- SingleR::MouseRNAseqData()

sample_filtered <- sample@assays$RNA@data

sample_singler <- SingleR(test = sample_filtered, ref = ref, assay.type.test=1,
                     labels = ref$label.main)

sample$celltype <- sample_singler$pruned.labels


cell_counts = as.data.frame(table(sample_singler$pruned.labels))
names(cell_counts) = c("Cell_Type", "Cell_Frequency")
cell_counts = dplyr::arrange(cell_counts, desc(Cell_Frequency))


sample$celltype <- sample_singler$pruned.labels

#name = "Tumor"
#outname = paste0(name,"_singleR_MouseRNAseqData.rds")
#saveRDS(sample, file = outname )


P1 = DimPlot(sample, reduction = "umap", split.by = "orig.ident", label=T,  label.size = 6)

P2 = DimPlot(sample, reduction = "umap", split.by = "orig.ident", group.by = "celltype", label=F)

subset = subset(x = sample, subset = celltype != "NA")

P3 = DimPlot(subset, reduction = "umap", split.by = "orig.ident", group.by = "celltype", label=T,  label.size = 4)

P4 = DimPlot(subset, reduction = "umap", split.by = "orig.ident", group.by = "celltype", label=F)

P1
P2

#P3
LabelClusters(plot = P4, id = "celltype", repel=TRUE)

```

![**Figure D**](/images/scRNAseq_integrate_2.png)  


```
cell_counts = dplyr::arrange(cell_counts, desc(Cell_Frequency))

cell_counts %>%
  kable(align = "c") %>%
  kable_styling(c("striped", "bordered"), full_width = T, fixed_thead = T) %>%
  row_spec(0, bold = T, color = "white", background = "#0571b0") %>%
  scroll_box(height = "500px", width = "800px")

```


### ## Assignment by Immunological Genome Project (ImmGen)
The ImmGen reference consists of microarray profiles of pure mouse immune cells. This is currently the most highly resolved immune reference - possibly overwhelmingly so, given the granularity of the fine labels.

```
name = "Tumor"
outname = paste0(name,"_SCT.rds")
sample = readRDS(file = outname )

ref <- SingleR::ImmGenData()

sample_filtered <- sample@assays$RNA@data

sample_singler <- SingleR(test = sample_filtered, ref = ref, assay.type.test=1,
                     labels = ref$label.main)

cell_counts = as.data.frame(table(sample_singler$pruned.labels))
names(cell_counts) = c("Cell_Type", "Cell_Frequency")
cell_counts = dplyr::arrange(cell_counts, desc(Cell_Frequency))


sample$celltype <- sample_singler$pruned.labels

#name = "Tumor"
#outname = paste0(name,"_singleR_ImmGenData.rds")
#saveRDS(sample, file = outname )

P1 = DimPlot(sample, reduction = "umap", split.by = "orig.ident", label=T,  label.size = 6)

P2 = DimPlot(sample, reduction = "umap", split.by = "orig.ident", group.by = "celltype", label=F)

subset = subset(x = sample, subset = celltype != "NA")

P3 = DimPlot(subset, reduction = "umap", split.by = "orig.ident", group.by = "celltype", label=T,  label.size = 6)

P4 = DimPlot(subset, reduction = "umap", split.by = "orig.ident", group.by = "celltype", label=F)

P1
P2

#P3
LabelClusters(plot = P4, id = "celltype", repel=TRUE)

cell_counts %>%
  kable(align = "c") %>%
  kable_styling(c("striped", "bordered"), full_width = T, fixed_thead = T) %>%
  row_spec(0, bold = T, color = "white", background = "#0571b0") %>%
  scroll_box(height = "500px", width = "800px")



```
![**Figure D**](/images/scRNAseq_integrate_3.png)  





