# Venn diagram
This script will generate Venn diagram with custom lists.

# Pre-requisites:
1. Create empty directory.
2. Create file with extension `*.venn` containing a venn list.
3. Generate multiple `*.venn` files, one for each list to compare.


# Quickstart:
```

# Syntax:
Rscript <VennDiagram_with_lists.R>  <DIRECTORY PATH>

# Example
Rscript    VennDiagram_with_lists.R    C:/Venn_lists/

```

# Code:
```
library('VennDiagram')
library(tidyverse)
library("optparse")
library(gplots)

## Accept the arguments and store in the variables
# args <- commandArgs(trailingOnly = TRUE)
# 
# path  = args[1]
# count = args[2]
# title = args[3]


option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path for the .venn files", metavar="character"),
  
  make_option(c("-c", "--count"), type="character", default="2", 
              help="count of .venn files", metavar="character"),
  
  make_option(c("-t", "--title"), type="character", default="Overlap", 
              help="Main title", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


setwd(opt$path)

myfiles = list.files(path=opt$path, pattern="*.venn", full.names=TRUE)

file_list = list()
mynames = vector()

for (myfile in myfiles) {
  
  mydata = read.table(file = myfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
  
  file_list = append(file_list, as.vector(mydata))
  
  mynames = append(mynames, colnames(mydata))
  
}

# Export the venn components
ItemsList  <- venn(file_list, show.plot=FALSE)
x = attr(ItemsList,"intersections")
ItemNames = names(x)

for (i in ItemNames) {
  
  data = x[[i]]
  data_name = str_replace(i, ":", "-")
  outname = paste(length(data),data_name, "items.txt", sep = "_")
  write.table(data, file = outname, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

if (opt$count == 2)
{

    venn.diagram(file_list, 
                 fill = c("#FEAD72", "#FED976"),
                 category.names = mynames,
                 output = TRUE ,
                 imagetype="png" ,
                 height = 2000 , 
                 width = 2000 , 
                 resolution = 300,
                 filename = "Venn.png",
                 main = opt$title,
                 main.cex = 1.5, main.just = c(1,-2),
                 alpha = c(0.5),
                 scaled = FALSE, euler.d = FALSE, # This is needed if sets are inclusive
                 cat.pos  = c(0, 180),
                 #cat.just     = list(c(0, -0.7), c(0, 0)),
                 cex = 1.5, cat.cex = 1.5,
                 print.mode = c("raw", "percent")
    )
  
}

if (opt$count == 3)
{
  
  venn.diagram(file_list, 
               fill = c("gray69", "#FEAD72", "#FED976"),
               category.names = mynames,
               imagetype="png" ,
               height = 2500 , 
               width = 2500 , 
               resolution = 300,
               filename = "Venn.png",
               main = opt$title,
               main.cex = 1.5, main.just = c(1,-2),
               alpha = c(0.5),
               scaled = FALSE, euler.d = FALSE, # This is needed if sets are inclusive
               cat.pos  = c(360, 0, 180),
               cat.just     = list(c(1, 3), c(0, 3), c(0.5, 0)),
               cex = 1.5, cat.cex = 1.5,
               print.mode = c("raw", "percent")
  )
  
}


  

```

