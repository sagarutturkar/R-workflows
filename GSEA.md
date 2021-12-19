# Gene Set Enrichment Analysis (GSEA):
One of the most popular methods to perform gene set enrichment analysis is [GSEA](http://www.gsea-msigdb.org/gsea/index.jsp) from the Broad Institute. GSEA is a computational 
method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two 
biological states (e.g. phenotypes). The method derives its power by focusing on gene sets, that is, groups of genes that share common 
biological function, chromosomal location, or regulation. The goal of GSEA is to determine whether members of a gene set S tend to occur 
towards the top or bottom of the ranked list L, which suggest gene set is correlated with a phenotypic (observed expression) signature.

## GSEA vs Pathway Analysis
GSEA algorithm is conceptually different than the pathway analysis performed through the tools like DAVID or IPA. Typically, the pathway 
analysis tools determine the overlap between user-supplied gene lists and the curated databases, looking for the overlaps that are bigger 
than that expected by random chance. The user-supplied gene list to DAVID or IPA is usually generated after selecting the genes that are 
differentially expressed at certain significance threshold. In contrast, the GSEA considers every gene that is expressed in the dataset.

## How to perform GSEA analysis:
GSEA can be run on [GenePattern web server](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/GSEA/14) 
or can be downloaded as [standalone application](http://software.broadinstitute.org/gsea/login.jsp) on your computer. Of course, the 
best place to learn about GSEA is the [user guide](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html)
described on the web site.  

GSEA analysis has been well described in various blog posts. Two most useful reference blogs posts are
1. [Mark Ziemann: Data analysis step 8: 
Pathway analysis with GSEA](http://genomespot.blogspot.com/2014/09/data-analysis-step-8-pathway-analysis.html).  
2. [GSEA explained by Tommy Tang](http://crazyhottommy.blogspot.com/2016/08/gene-set-enrichment-analysis-gsea.html).  

GSEA has two run options:  
1. Use the raw counts and allow GSEA to perform the ranking  OR
2. Provide a pre-ranked list of genes

## Why GSEA pre-ranked is preferred:
![**Figure A**](/images/GSEA_3.PNG)  
A complete description is [available here.](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/GSEAPreranked/1)  

In general, providing the pre-ranked list is preferred option for RNA-Seq analysis. It is recommended to perform ranking of all the genes 
detected in DESeq2 **(without any significance filtering)**. Although multiple methods to perform ranking are available, the most commonly 
used method for ranking is calculation of expression matrix with formulae below and then sorted in descending order.

>signed fold-change * -log10pvalue  

Creating the pre-ranked GSEA list can be performed withing DESeq2 code as shown below:
```
library(DESeq2)

# Prepare DESeq2 Object
dds<- DESeqDataSetFromMatrix(countData = countData,
                             colData = samples,
                             design = ~ treatment)
                             
                             
dds <- DESeq(dds)
resultsNames(dds)

#Generate results for given comparison
res1 <- results(dds,contrast=list("treatmentTreatment", "treatmentControl"))
res1$Gene_ID = rownames(res1)

#Generate Results Table
res1_table<-as.data.frame(res1)
res1_table$fc <- 2^res1_table$log2FoldChange
res1_table  <-  res1_table[c(1,2,3,4,8,5,6,7)]
head(res1_table)

#Generate input file for GSEA
# Extracted columns (GeneID, log2foldchange, pvalue)
# Please note the order could be different in your data
GSEA = res1_table[c(11,7,3)] 

#calculate expression matric with recommended formula
GSEA$Rank = sign(GSEA$log2FoldChange) * -log10(GSEA$pvalue)
GSEA_export = GSEA[c(1,5)]
GSEA_export = GSEA_export[order(GSEA_export$Rank, decreasing = TRUE),]
GSEA_export = GSEA_export[complete.cases(GSEA_export),] # Remove NA values (when pvalue is NA)

write.table(GSEA_export, file = "GSEA.txt", sep = "\t", row.names = FALSE, quote = FALSE)

```

The above code will genarte a ranked genelist for GSEA. This should be saved without headers in the file `GSEA.rnk`.  


Example **GSEA.rnk** file
```
IFIT1	230.739612420274
OASL	220.750809016291
IFIT3	163.61074021669
IFIT2	156.489754913417
CCL20	136.358435563505
BIRC3	100.289570931789
HIST1H2AC	84.3196781855268
OAS1	83.9008611266189
TUBA3D	71.3197748570235
ISG15	69.7715037404975
CDKN1A	68.2509517932692
HIST2H2BE	61.0363170205825
HMOX1	57.8609010201463
DDX58	55.1975444793238
```
A pre-ranked gene-list supplied to GSEA is compared against the MSigDB gene sets that are divided into 8 major collections.  

![**Figure 1**](/images/MSigDB_collections.png)  

## GSEA Results:
The best way to navigate through GSEA results is to start with the report provided in the HTML format. The report file is named as index.html for each gene set. The basic elements of the report are summarized below.
![**Figure 2**](/images/GSEA_1.png)  

> **Note:** The report described in Figure 2 is an example report and descriptions are intended only for the most important elements of the report. Please follow the [Guide to interpret results](http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results) for complete report description.

## Important components of GSEA results:
1. **The Enrichment Score (ES)**, reflects the degree to which a gene set is overrepresented at the top or bottom of a ranked list of genes. A positive ES indicates gene set enrichment at the top of the ranked list; a negative ES indicates gene set enrichment at the bottom of the ranked list.
2. **The Normalized Enrichment Score (NES)** is the primary statistic for examining gene set enrichment results. By normalizing the enrichment score, GSEA accounts for differences in gene set size and in correlations between gene sets and the expression dataset; therefore, the normalized enrichment scores (NES) can be used to compare analysis results across gene sets. 
3. **The nominal p-value** estimates the statistical significance of the enrichment score for a single gene set. However, when evaluating multiple gene sets, a correction needs to perform for gene set size and multiple hypothesis testing. Because p-value is not adjusted for either, it is of limited value when comparing gene sets.
4. **The False Discovery Rate (FDR)** is the estimated probability that a gene set with a given NES represents a false positive finding. For example, an FDR of 25% indicates that the result is likely to be valid 3 out of 4 times. The GSEA analysis report highlights enrichment gene sets with an FDR of less than 25% as those most likely to generate interesting hypotheses and drive further research but provides analysis results for all analyzed gene sets.  

## Enrichment plot example:
An enrichment plot is created for each gene set within the collection and accessible through snapshot of enrichment results. GSEA calculates the ES by walking down the ranked list of genes, increasing a running-sum statistic when a gene is in the gene set and decreasing it when it is not. The magnitude of the increment depends on the correlation of the gene with the phenotype. The ES is the maximum deviation from zero encountered in walking the list. A positive ES indicates gene set enrichment at the top of the ranked list; a negative ES indicates gene set enrichment at the bottom of the ranked list.  

In the analysis results, the enrichment plot provides a graphical view of the enrichment score for a gene set. An example enrichment plot is shown in Figure 3.  

![**Figure 3**](/images/GSEA_2.png) 

Enrichment plot description:  
1.	The top portion of the plot shows the running ES for the gene set as the analysis walks down the ranked list. The score at the peak of the plot (the score furthest from 0.0) is the ES for the gene set. Gene sets with a distinct peak at the beginning (such as the one shown here) or end of the ranked list are generally the most interesting.
2.	The leading edge subset of a gene set is the subset of members that contribute most to the ES. 
3.	The bottom portion of the plot shows the value of the ranking metric as you move down the list of ranked genes. The ranking metric measures a gene’s correlation with a phenotype. The value of the ranking metric goes from positive to negative as you move down the ranked list.

## Custom GSEA plot:
![**Figure 4**](/images/GSEA_4.png)  

## References:
1.	Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A 102, 15545-15550, doi:10.1073/pnas.0506580102 (2005).  

2.	Huang da, W., Sherman, B. T. & Lempicki, R. A. Systematic and integrative analysis of large gene lists using DAVID bioinformatics resources. Nat Protoc 4, 44-57, doi:10.1038/nprot.2008.211 (2009).  

3.	Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550, doi:10.1186/s13059-014-0550-8 (2014).  


