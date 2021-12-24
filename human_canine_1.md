# Goal
Goal of this analysis is to compare the human and canine transcriptomic data. Here, we would like to compare the  gene-expression between two species instead of looking at the individual pathways. 

```
require(data.table)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(pals)
library(xlsx)
library(gplots)
library(knitr)
library(kableExtra)
library(ggpubr)
library(janitor)
library(tables)
library(rJava)
```


```
# colors
RD   =  brewer.pal(n = 5, name = "RdBu")
BLUE = brewer.pal(n = 5, name = "Blues")
RYD  = brewer.pal(n = 5, name = "RdYlBu")

# https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

my_color  =  colorRampPalette(c("#4d9221", "black", "#e31a1c"))(10)

```

# Approach
We first started with combining the gene-expression data (TPM, considering it is normalized value by gene) across human and canine (using orthologs). However, several other factors affect the gene-expression values (see this paper)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7373998/]. We did run CLUST tool on combined (human and canine) TPM matrix, which generated few patterns, but most were skewed and did not stand when looked through heatmap. The reason was the vast differece in the TPM values across two species.

Therefore, another workaround was tested, where we compared each species data against the respective normal using CLUST. We primarily focused on comparison between three groups:  
1. Canine_Early  
2. human_luminal  
3. human_basal  

We ran CLUST for individual species with above groups and respective normal. For canine data, Nor-early group was added as CLUST works best with three groups data. Below are the CLUST profiles for each individual species.


```
#load_TPM
canine_TPM = read.table(file = "Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/FPKM/TPM_copy.TXT", header = T, sep = "\t", quote = "", check.names = F)

canine_TPM[1:5,1:5]
dim(canine_TPM)

canine_metadata <- read.xlsx("Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/human_canine_metadata_ddhawan.xlsx", sheetName = "Canine_metadata", header = T, stringsAsFactors = F)

canine_32PT = canine_metadata$sample

canine_TPM = canine_TPM %>%
  dplyr::select(all_of(c("Gene_ID", canine_32PT)))


BLCA_TPM = readRDS(file = "Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/Human_clust/GTEX_BLCA_TPM_18451genes.rds") 

BLCA_TPM[1:5,1:5]
dim(BLCA_TPM)

```


```
# load_CLUST_data
#human
Clust_human = read.table(file = "Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/Human_clust/By_subtype/Clusters_Objects.tsv", sep = "\t", header = F, check.names = F, as.is = T, skip = 2, na.strings = c(c("","NA")))
names_human = paste0("C",seq(from=0, to = length(Clust_human)-1))
names(Clust_human) = names_human
head(Clust_human)
human_genes = na.omit(Clust_human$C0)

#dog
Clust_dog = read.table(file = "Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/Canine_clust/By_stage/Clusters_Objects.tsv", sep = "\t", header = F, check.names = F, as.is = T, skip = 2, na.strings = c(c("","NA")))
names_dog = paste0("C",seq(from=0, to = length(Clust_dog)-1))
names(Clust_dog) = names_dog
head(Clust_dog)

dog_genes = na.omit(c(Clust_dog$C0, Clust_dog$C1, Clust_dog$C2, Clust_dog$C3))



ortho_table = read.table(file = "Z:/PCCR/Lilly/Universal_data/For_GSEA/dog_to_human_orthologs_12142020.TXT", header = T, sep = "\t", quote = "", na.strings = c(c("","NA")))

ortho_table = ortho_table %>%
  drop_na(Human_ENSEMBL)

dog_table = as.data.frame(dog_genes)
dog_table = dplyr::left_join(dog_table, ortho_table, by = c("dog_genes" = "Dog_ENSEMBL"))
dog_table = dog_table %>%
  drop_na(Human_ENSEMBL) %>%
  dplyr::arrange(dog_genes)

human_canine_intersect_clust2 = intersect(human_genes ,dog_table$Human_ENSEMBL)

ortho_table_358 = dog_table %>%
  dplyr::filter(Human_ENSEMBL %in% human_canine_intersect_clust2)

write.table(ortho_table_358, file = "CLUST2_intersect.txt", sep = "\t", quote = F, row.names = F)
```

### CLUST results for human data
![**Figure A**](/images/human_clust2.png)  

### CLUST results for canine data
![**Figure B**](/images/canine_clust2.png)  


### Selecting patterns:
We are interested in comparison of (Canine_Early, human_luminal, human_basal) groups. We determined the CLUST patterns that are up-regulated in group of interest as compared to normal data. From human data, **Cluster C0** comprised of 2544 genes that are up-regulated in (human_luminal, human_basal) groups. Similarly, for canine data, **Cluster C0, C1, C2, C3** contain 1956 genes that are up-regulated in canine_early group.

### Comparing data across species
So far, the data has been compared within the species and we have determined genes that show up-regulation pattern against respective normal. To compare the data between species, we need to determine the genes that are shared across species (i.e. orthologs), and then check how many of those orthologs are present in above two gene-lists (i.e. 2544 human genes and 1956  canine genes).

We have 358 genes that have shared orthologs between two species. These genes will be further probed in Heatmap.
```
library('VennDiagram')

my_list = list(human = human_genes, dog = dog_table$Human_ENSEMBL)

    venn.diagram(my_list, 
                 fill = c("#FEAD72", "#FED976"),
                 category.names = names(my_list),
                 output = TRUE ,
                 imagetype="png" ,
                 height = 2000 , 
                 width = 2000 , 
                 resolution = 300,
                 filename = "Venn.png",
                 main = "Venn",
                 main.cex = 1.5, main.just = c(1,-2),
                 alpha = c(0.5),
                 scaled = FALSE, euler.d = FALSE, # This is needed if sets are inclusive
                 cat.pos  = c(0, 180),
                 #cat.just     = list(c(0, -0.7), c(0, 0)),
                 cex = 1.5, cat.cex = 1.5,
                 print.mode = c("raw", "percent")
    )
```
![**Figure C**](/images/HC1_Venn1.png)  

## Heatmap of 358 orthologous genes in human data
```
human_metadata <- read.xlsx("Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/human_canine_metadata_ddhawan.xlsx", sheetName = "Human_metadata", header = T, stringsAsFactors = F)

human_anno = human_metadata %>%
  dplyr::select(all_of(c("sample", "TCGA408")))

human_anno = column_to_rownames(human_anno, "sample")

human_Basal_squamous = human_metadata %>%
  dplyr::filter(TCGA408 == "Basal_squamous") %>%
  dplyr::select(sample)


human_Luminal = human_metadata %>%
  dplyr::filter(TCGA408 == "Luminal") %>%
  dplyr::select(sample)


human_normal = human_metadata %>%
  dplyr::filter(TCGA408 == "Normal") %>%
  dplyr::select(sample)



breaksList = seq(-1, 1, by = 0.2)

# Human
heat_data = BLCA_TPM %>%
  dplyr::filter(Gene_ID %in% ortho_table_358$Human_ENSEMBL) %>%
  dplyr::select(all_of(c("Gene_ID", human_normal$sample, human_Basal_squamous$sample, human_Luminal$sample)))

write.table(heat_data, file = "human_heat_data_clust2.txt", sep = "\t", quote = F, row.names = F)


heat_data = column_to_rownames(heat_data, "Gene_ID")


my_gaps = c(21, 105)

png(filename = "human_clust2_heatmap.png", units = "in", width = 12, height = 6, bg = "white", res = 300)


pheatmap(heat_data, show_rownames = F, main = "Heatmap - BLCA", show_colnames = F, 
         cluster_rows = T, cluster_cols = F,
         color = my_color, breaks = breaksList,
         cellheight=1, cellwidth=4, na_col = "grey",
         scale = "row", gaps_col = my_gaps, 
         treeheight_row = 0, annotation_col = human_anno
)

dev.off()

```

![**Figure D**](/images/human_clust2_heatmap.png)  

## Heatmap of 358 orthologous genes in canine data
```
canine_metadata <- read.xlsx("Z:/PCCR/knappd/Screening_study_2021/Human_canine_clust_December2021/human_canine_metadata_ddhawan.xlsx", sheetName = "Canine_metadata", header = T, stringsAsFactors = F)

# dog
canine_early = canine_metadata %>%
  dplyr::filter(Stage == "Early") %>%
  dplyr::select(sample)


canine_normal = canine_metadata %>%
  dplyr::filter(Stage == "Normal") %>%
  dplyr::select(sample)


heat_data = canine_TPM %>%
  dplyr::filter(Gene_ID %in% unique(ortho_table_358$dog_genes)) %>%
  dplyr::select(all_of(c("Gene_ID", canine_normal$sample, canine_early$sample)))

write.table(heat_data, file = "canine_heat_data_clust2.txt", sep = "\t", quote = F, row.names = F)


heat_data = column_to_rownames(heat_data, "Gene_ID")

my_gaps = c(4)

canine_anno = canine_metadata %>%
  dplyr::select(all_of(c("sample", "Stage")))

canine_anno = column_to_rownames(canine_anno, "sample")

png(filename = "canine_clust2_heatmap.png", units = "in", width = 12, height = 6, bg = "white", res = 300)

pheatmap(heat_data, show_rownames = F, main = "Heatmap - Canine", show_colnames = F, 
         cluster_rows = T, cluster_cols = F,
         color = my_color, breaks = breaksList,
         cellheight=1, cellwidth=4, na_col = "grey",
         scale = "row", gaps_col = my_gaps,
         annotation_col = canine_anno, treeheight_row = 0
)

dev.off()
```

![**Figure E**](/images/canine_clust2_heatmap.png)  

## Human and Canine figure side-by-side (excluding normals)  
> Note: **The gene-order is not same yet**, but we can keep the gene order same. I only want to show a conceptual example here.
![**Figure F**](/images/HC1_Clust2_heatmap_combined.png) 
Most of these genes are up-regulated in both human and canine as compared to their normal data. We get some outliers, but once we set the sample order, and further select the genes, the Heatmap would look much cleaner.


# Future directions:
Remember, I started this analysis with the goal of finding the human-canine homologous genes that show similar expression (up-regulation) across three groups (Canine_Early, human_luminal, human_basal), and we get 358 genes with such pattern. In this way, we could probe different desired pattern (e.g. we could look for genes that show up-regulation in Canine_Early and human_basal, but down-regulated in human_luminal)  provided such pattern is deermined by CLUST.

> Note: Please note this is a supervised analysis where we are probing the data as per the desired condition. Our starting point was 2544 human genes and 1956  canine genes, but we could only determined 358 (~ 8%) genes with desired pattern. If our starting gene-set through CLUST is small, we may only get a handful of such genes. It really depends on how many orthologs are shared.

We could perform pathway analysis on these 358 genes, select the specific pathways of interest and then make  heatmap with  genes within specific pathways (hoping that pattern stays intact).  

