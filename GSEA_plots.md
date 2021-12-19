# GSEA
GSEA is commonly employed pathway analysis technique. I have a [separate write-up on principles of GSEA](https://github.com/sagarutturkar/R-workflows/blob/main/GSEA.md). 

# GSEA custom plots
I generated a custom R-script that reads the output from GSEA standalone tool and generates the **barplot** and **dotplot** for top enriched gene-sets **(pvalue < 0.05).**

# Quickstart
```
# Syntax:
Rscript <GSEA_plot.R>  <GSEA STANDARD OUTPUT DIRECTORY PATH>  <TAG>

# Example
Rscript    GSEA_plot.R    C:/Users/sutturka/gsea_home/output/oct16/GSEA_test.GseaPreranked.1634401613095/  C1

```

# Code:
```
require(data.table)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(xlsx)
library(knitr)
library(kableExtra)
library(ggplot2)

## Syntax Example
## system("Rscript Z:/Scripts/R_scripts/GSEA_plot_2021.R    C:/Users/sutturka/gsea_home/output/oct16/GSEA_test.GseaPreranked.1634401613095/  C1")


## Accept the arguments and store in the variables
args <- commandArgs(trailingOnly = TRUE)

dir_path = args[1]
tag      = args[2] 

outname_barplot = paste0(tag,"_barplot.PNG")
outname_dotplot = paste0(tag,"_dotplot.PNG")


# function for GSEA data processing
# Takes GSEA output path as input 
# Reads the TSV files (^gsea_report_for_na.*tsv$)
# Return the GSEA data frame with select columns

process_GSEA <- function(dir_path) {
  
  file_paths = list.files(path = dir_path, pattern = "^gsea_report_for_na.*tsv$", full.names = T)
  GSEA_df = data.frame()
  
  cat("\n\nRead following files: \n")
  
  for (myfile in file_paths) {
    
    print(myfile)
    my_df = read.table(file = myfile, sep = "\t", header = T, stringsAsFactors = F)
    
    my_cols = c("NAME", "SIZE", "ES", "NES", "NOM.p.val", "FDR.q.val")
    my_df = dplyr::select(my_df, all_of(my_cols))
    
    names(my_df) = c("pathway", "SIZE", "ES", "NES", "pval", "qval")
    dim(my_df)
    
    GSEA_df = rbind(GSEA_df, my_df)
    
  }
  
  cat("\n\nDimensions of original GSEA table is: \n")
  print(dim(GSEA_df))
  return(GSEA_df)
  
}


# function for GSEA barplot
# Takes GSEA table (data frame) as input and generate the barplot at pval < 0.05
plot_gsea_barplot <- function(GSEA.results) {
  
  GSEA.results$direction = if_else(GSEA.results$NES >= 0, "Positive", "Negative")
  GSEA.results$direction = factor(GSEA.results$direction, levels = c("Positive", "Negative"))
  
  fgsea.sig.pos = GSEA.results %>%
    dplyr::filter(pval <= 0.05, direction == "Positive")
  
  fgsea.sig.neg = GSEA.results %>%
    dplyr::filter(pval <= 0.05, direction == "Negative")
  
  fgsea.sig = rbind(fgsea.sig.pos, fgsea.sig.neg)
  
  cat("\n\nDimensions of filtered (pval < 0.05) GSEA table is: \n")
  print(dim(fgsea.sig))
  
  write.table(fgsea.sig, file = "test.txt", sep = "\t", row.names = F)
  
  p = ggplot(fgsea.sig, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(aes(fill = direction)) +
    coord_flip( ylim = c(-2,2) ) +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Gene Set Enrichment (Pvalue < 0.05)", fill = "NES") + 
    theme_minimal() + 
    #scale_fill_manual(values=c("#999999", "#E69F00")) +
    scale_fill_manual(values=c("#00C0B8", "#F8766D")) +
    theme(legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20))
  
  return(p)
  
}



plot_gsea_dotplot <- function(GSEA.results) {
  
  GSEA.results$direction = if_else(GSEA.results$NES >= 0, "Positive", "Negative")
  GSEA.results$direction = factor(GSEA.results$direction, levels = c("Positive", "Negative"))
  
  fgsea.sig.pos = GSEA.results %>%
    dplyr::filter(pval <= 0.05, direction == "Positive")
  
  fgsea.sig.neg = GSEA.results %>%
    dplyr::filter(pval <= 0.05, direction == "Negative")
  
  fgsea.sig = rbind(fgsea.sig.pos, fgsea.sig.neg)
  
  cat("\n\nDimensions of filtered (pval < 0.05) GSEA table is: \n")
  print(dim(fgsea.sig))
  
  fgsea.sig$significance = 1 / fgsea.sig$pval
  
  write.table(fgsea.sig, sep = "\t", row.names = F, file = "test.txt")
  
  p = ggplot(data=fgsea.sig) + 
    aes(x=pathway,  y=NES, size=pval, color=direction) + 
    geom_point(alpha = 0.8)+
    ylab("Normalized Enrichment Score") +
    ggtitle("Gene Set Enrichment (Pvalue < 0.05)") +
    xlab("Pathway") +
    coord_flip() +
    theme_minimal() +
    scale_color_manual(values = c("Positive" = "#1d91c0", "Negative" = "#fe9929")) +
    scale_y_continuous(limits = c(-2, 2)) +
    scale_size_continuous(trans = "reverse") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed', colour = "grey")
    ) +
    theme(axis.line = element_line(color="black", size = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
   return(p)
  
}



cat("Reading data from Directory: \n")
print(dir_path)
GSEA = process_GSEA(dir_path)
barplot = plot_gsea_barplot(GSEA)
dotplot = plot_gsea_dotplot(GSEA)


cat("\n\nWrote GSEA barplot to file:", outname_barplot, "\n")

png(filename = outname_barplot, width = 3000, height = 3000, pointsize = 1,  bg = "white", res = 300)
print(barplot)
dev.off()

cat("\n\nWrote GSEA dotplot to file:", outname_dotplot, "\n")

png(filename = outname_dotplot, width = 3000, height = 3000, pointsize = 1,  bg = "white", res = 300)
print(dotplot)
dev.off()

```

# Output:
![**Dotplot**](/images/C1_dotplot.PNG) 

![**Barplot**](/images/C1_barplot.PNG) 



