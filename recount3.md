# Download custom datasets (GTEX/TCGA) using recount3 package:
This workflow will define how to download custom datasets from [recount3: uniformly processed RNA-seq](http://rna.recount.bio/) project. Recount3 has pre-processed RNAseq data from common databases like GTEX/TCGA as well as several public SRA projects. The currentflow can be customized to select any data of interest and same can be downloaded easily.

```
# Library
library("recount3")
library(recount)
library(dplyr)
library(tidyverse)
library(stringr)
library(xlsx)
```

## Explore the recount3 and find the data of interest:
Here, we are selected the GTEX and TCGA-BLCA data to download. This includes 2
```
human_projects <- available_projects()


BLCA <- subset(
    human_projects,
    project == "BLCA" & project_type == "data_sources"
)


rse_obj <- create_rse(BLCA)

assay(rse_obj, "counts") <- transform_counts(rse_obj)

assayNames(rse_obj)

assays(rse_obj)$TPM = recount::getTPM(rse_obj)

BLCA_names = rse_obj$tcga.tcga_barcode
```

## Download data of interest
### GTEX (tissue type bladder) raw-counts
```
# GTEX raw-counts

gtex_bladder_counts = assays(rse_obj)$counts
gtex_bladder_counts = as.data.frame(gtex_bladder_counts)
gtex_bladder_counts = rownames_to_column(gtex_bladder_counts, "Gene_ID")
gtex_bladder_counts[1:5,1:5]
dim(gtex_bladder_counts)

write.table(gtex_bladder_counts, file = "gtex_bladder_Rawcounts.txt", sep = "\t", quote = F, row.names = F)
```

### GTEX (tissue type bladder) Transcripts-per-million (TPM)
```
gtex_bladder_TPM = assays(rse_obj)$TPM
gtex_bladder_TPM = as.data.frame(gtex_bladder_TPM)
gtex_bladder_TPM = rownames_to_column(gtex_bladder_TPM, "Gene_ID")
gtex_bladder_TPM[1:5,1:5]
dim(gtex_bladder_TPM)


gtex_bladder_TPM$rowsum = rowSums(gtex_bladder_TPM[2:length(gtex_bladder_TPM)])
gtex_bladder_TPM = separate(gtex_bladder_TPM, "Gene_ID", c("Gene_ID"))
gtex_bladder_TPM[1:5,1:5]

# sort  data by rowsums
# remove duplicate gene IDs (low counts one will be removed)
gtex_bladder_TPM = gtex_bladder_TPM %>%
  dplyr::arrange(desc(rowsum)) %>%
  dplyr::distinct(Gene_ID, .keep_all = T) %>%
  dplyr::select(-c("rowsum"))

saveRDS(gtex_bladder_TPM, file = "gtex_bladder_TPM_noDups.rds")
```

### TCGA-BLCA raw-counts
```
BLCA_counts = assays(rse_obj)$counts
BLCA_counts = as.data.frame(BLCA_counts)
BLCA_counts = rownames_to_column(BLCA_counts, "Gene_ID")
BLCA_counts[1:5,1:5]
dim(BLCA_counts)

colnames(BLCA_counts) = c("Gene_ID", BLCA_names)

human_metadata <- read.xlsx("Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/human_canine_metadata_ddhawan.xlsx", sheetName = "Human_metadata", header = T, stringsAsFactors = F)

BLCA_counts = dplyr::select(BLCA_counts, all_of(c("Gene_ID",human_metadata$Barcode)))


write.table(BLCA_counts, file = "BLCA_Rawcounts_234PT.txt", sep = "\t", quote = F, row.names = F)


system("RScript Z:/Scripts/R_scripts/filter_countMatrix_edgeR.R   Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/GTEX_BLADDER/BLCA_Rawcounts_234PT.TXT  BLCA_Rawcounts_234PT.filtered.TXT")

```
### TCGA-BLCA raw-counts (filter and remove duplicate genes)
```
gtex_filtered =  read.table(file = "gtex_bladder_Rawcounts.filtered.TXT", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "",  check.names = F, as.is = T)

gtex_filtered = separate(gtex_filtered, Gene_ID, into = c("Gene_ID"), sep = "\\.", extra = "drop")



BLCA_filtered =  read.table(file = "BLCA_Rawcounts_234PT.filtered.TXT", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "",  check.names = F, as.is = T)

BLCA_filtered = separate(BLCA_filtered, Gene_ID, into = c("Gene_ID"), sep = "\\.", extra = "drop")

gtex_filtered[1:5,1:5]
BLCA_filtered[1:5,1:5]

dim(gtex_filtered)
dim(BLCA_filtered)

length(unique(gtex_filtered$Gene_ID))
length(unique(BLCA_filtered$Gene_ID))


# Change BLCA barcodes to first 12 letters only
BLCA_barcodes = colnames(BLCA_filtered[2:length(BLCA_filtered)])
BLCA_barcodes = str_sub(BLCA_barcodes, 1,12)
colnames(BLCA_filtered) = c("Gene_ID", BLCA_barcodes)

#get intersection of GTEX BLCA
gtex_BLCA_intersect = intersect(BLCA_filtered$Gene_ID, gtex_filtered$Gene_ID)
length(gtex_BLCA_intersect)

BLCA_filtered = dplyr::filter(BLCA_filtered, Gene_ID %in% gtex_BLCA_intersect)
gtex_filtered = dplyr::filter(gtex_filtered, Gene_ID %in% gtex_BLCA_intersect)
dim(gtex_filtered)
dim(BLCA_filtered)

gtex_BLCA_counts = dplyr::left_join(BLCA_filtered, gtex_filtered)
gtex_BLCA_counts[1:5,1:5]
dim(gtex_BLCA_counts)

write.table(gtex_BLCA_counts, file = "GTEX_BLCA_Rawcounts_255PT_18451Genes.txt", sep = "\t", quote = F, row.names = F)

```


