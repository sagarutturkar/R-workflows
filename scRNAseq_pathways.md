# Enrichment Analysis Results for Cluster Biomarkers

```
require(data.table)
library(clusterProfiler)
library(org.Cf.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(cowplot)
library(ReactomePA)
library(enrichplot)
library(org.Mm.eg.db)
library(tidyverse)
library(dplyr)
```
## Load data and marker genes
```
name = "Tumor"
outname = paste0(name,"_markers.rds")
markers = readRDS(outname)

cid = unique(markers$cluster)

my_list = list()

for (c in cid) {
  
  x = markers %>%
  dplyr::filter(cluster == c) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
  .$gene %>% unique()
  
  ENT_ID = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  
  my_list[[c]] = ENT_ID$ENTREZID
}

lapply(my_list, length)

```

## KEGG Enrichment
```
cid_KEGG <- compareCluster(geneCluster = my_list, organism = "mmu",
                     fun = "enrichKEGG")

dotplot(cid_KEGG)

write.csv(cid_KEGG, file="KEGG_by_cluster.csv", row.names = FALSE, quote = F)

```
![**KEGG**](/images/scRNAseq_pathways_1.png) 


## REACTOME Enrichment
```
cid_RCT <- compareCluster(geneCluster = my_list, organism = "mouse",
                     fun = "enrichPathway")

dotplot(cid_RCT) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))

write.csv(cid_RCT, file="REACTOME_by_cluster.csv", row.names = FALSE, quote = F)
```
![**REACTOME**](/images/scRNAseq_pathways_2.png) 

## GO.BP Enrichment
```
cid_Go.BP <- compareCluster(geneCluster = my_list,  OrgDb = org.Mm.eg.db, ont = "BP",
                     fun = "enrichGO")

dotplot(cid_Go.BP) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))

write.csv(cid_Go.BP, file="GO.BP_by_cluster.csv", row.names = FALSE, quote = F)
```
![**GO.BP**](/images/scRNAseq_pathways_3.png) 
