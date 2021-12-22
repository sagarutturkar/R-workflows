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


## Determine filtering parameters (custom for each data):
```

ori_cell_count = dim(sample@assays$RNA)[2]

mylist = list()
mt_cutoff = c(10, 15, 20, 25)
temp_table = data.frame()

for (i in mt_cutoff) {
  
  x <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < i)
  cell_count = dim(x@assays$RNA)[2]
  myvec = c(ori_cell_count, cell_count)
  print(myvec)
  temp_table = rbind(temp_table, myvec)
  
}

colnames(temp_table) = c("Original Cell Count", "Filtered Cell Count")
rownames(temp_table) = c("MT cutoff 10%", "MT cutoff 15%", "MT cutoff 20%", "MT cutoff 25%")


temp_table %>%
  kable(align = "c") %>%
  kable_styling(c("striped", "bordered"), full_width = T, fixed_thead = T) %>%
  row_spec(0, bold = T, color = "white", background = "#0571b0") %>%
  scroll_box(height = "250px", width = "450px")

```

![**Figure 3**](/images/Seurat_3.png)  

```
# After testing various cutoffs (table above) - we decided to use following filtering parameters
upper_cutoff = 6000
lower_cutoff = 200
MT_cutoff    = 15  

infile = paste0(name,"_base.rds")

sample = readRDS(file = infile)

sample <- subset(sample, subset = nFeature_RNA > lower_cutoff & nFeature_RNA < upper_cutoff & percent.mt < MT_cutoff)

outfile = paste0(name,"_QC.rds")
saveRDS(sample, file = outfile)


```



