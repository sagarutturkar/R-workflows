# DoubletFinder
[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) is an R package that predicts doublets in single-cell RNA sequencing data.

DoubletFinder can be broken up into 4 steps:  

1. Generate artificial doublets from existing scRNA-seq data  
2. Pre-process merged real-artificial data  
3. Perform PCA and use the PC distance matrix to find each cell's proportion of artificial k nearest neighbors (pANN)  
4. Rank order and threshold pANN values according to the expected number of doublets  

```
# Read seuarat object after QC
infile = paste0(name,"_SCT.rds")
sample = readRDS(file = infile)
```

# Important parameters in DoubletFinder
1. seu ~ This is a fully-processed Seurat object (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been run).  
2. PCs ~ The number of statistically-significant principal components, specified as a range (e.g., PCs = 1:10)  
3. pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 25%, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).  
4. pK ~ This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values should be estimated using the strategy described below.  
5. nExp ~ This defines the pANN threshold used to make final doublet/singlet predictions. This value can best be estimated from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated proportion of homotypic doublets.

We have selected these parameters as below:  

1. Seu fully-processed Seurat object  

2. PCs = 1:10. [DoubletFinder GitHub](https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/91) suggested adding extra PCs (1:50) doesn't add too much information. This was tested for one sample and then PCs were selected as 1:10.  

3. pN value was set to default 25%.  

4. pK - Optimal pK value was determined for each dataset using the graph below.  

5. nExp (Predicted doublet rate) was set based on doublet rate estimation table available through Cell-ranger. But these could also be caluclated [manually](https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76) for each sample.  

![**Figure A**](/images/DoubletFinder_1.PNG)  

```
sweep.res <- paramSweep_v3(sample, PCs = 1:10, sct = TRUE)
```

To maximize the accuracy of DoubletFinder predictions, we sought a ground-truth-agnostic metric that coincides with pK values that maximize AUC in Cell Hashing and Demuxlet data. 
Mean-variance normalized bimodality coefficient (BCmvn) achieves this goal, featuring a single, easily-discernible maximum at pK values that optimize AUC.


```

sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sample@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(sample@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
```
## Select the parameters as per the data **(custom for each data)**
```
pk_cutoff = 0.005
pn_cutoff = 0.25

seu <- doubletFinder_v3(sample, PCs = 1:10, pN = pn_cutoff, pK = pk_cutoff, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)

outname = paste0(name,".RData")
save.image(file = outname)

```

## Plot data
```
name = "sample"

outname = paste0(name,".RData")
load(file = outname)

bcmvn <- find.pK(sweep.stats)

x = as.data.frame(table(seu@meta.data$DF.classifications_0.25_0.005_843))

x %>%
  kable(align = "c") %>%
  kable_styling(c("striped", "bordered"), full_width = T, fixed_thead = T) %>%
  row_spec(0, bold = T, color = "white", background = "#0571b0") %>%
  scroll_box(height = "175px", width = "800px")

DimPlot(seu, split.by = "DF.classifications_0.25_0.005_843")
```

