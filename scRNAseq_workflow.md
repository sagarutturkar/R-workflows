# scRNAseq workflow

## Process raw data through Cell-Ranger
```
# Cellranger 

cellranger count --id=G1-spleen \
--transcriptome=/1_Cell_ranger/reference/refdata-gex-mm10-2020-A \
--fastqs=/1_Cell_ranger/input/G1-spleen \
--sample=G1 \
--localcores 32  \
--expect-cells=5000

```

## QC through Seurat

Seurat allows to easily explore QC metrics and filter cells based on various quality criteria. A few commonly used QC metrics are given below:

1. Transcript (nCount) and gene (nFeature) abundance
	- Low-quality cells or empty droplets will often have very few genes/transcripts
	- Cell doublets or multiplets may exhibit an aberrantly high gene count
2. The percentage of reads that map to the mitochondrial genome
	- Low-quality / dying cells often exhibit extensive mitochondrial contamination
	- Mitochondrial QC metrics is calculated as percentage of counts originating from a set of mitochondrial genes (i.e. all genes starting with MT)
3. Ribosomal contents
	- Less than 50% ribosomal content is often preferred

```
# import cell-ranger matrix
sample <- Read10X(data.dir = "/1_Cell_ranger/output/G1-spleen/outs/filtered_feature_bc_matrix")
sample <- CreateSeuratObject(counts = sample, project = "Tumor", min.cells = 3, min.features = 200)

# Determine mitochondrial and ribosomal percentage
sample[["percent.mt"]] = PercentageFeatureSet(sample, pattern = "^mt-")
sample[["percent.ribo"]] = PercentageFeatureSet(sample, pattern='RP[SL][[:digit:]]')

```

## visualization of QC matrices
```
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

VlnPlot(sample, features = c('percent.mt','percent.ribo'), pt.size = 0)

FeatureScatter(object=sample, feature1='percent.mt', feature2='nCount_RNA', pt.size = 0.1, plot.cor = F)

```
### Transcript (nCount), gene (nFeature) and Low-quality / dying cells (mt) abundance.
 
![**Figure 1**](/images/Seurat_1.png)  

### To identify dead cells/debris and doublets, mitochondrial contents and transcript abundance are plotted together

![**Figure 2**](/images/Seurat_2.png)  





