#installing packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
BiocManager::install("org.Hs.eg.db")
BiocManager::install("annotate")
#Load packages
library('org.Hs.eg.db')

library(annotate)

#Load Magma results

genes <- read.csv("MAGMA_Results.txt", sep="")
genes
genes$ZSTAT
#load the keys?
## select() interface: 
## Objects in this package can be accessed using the select() interface 
## from the AnnotationDbi package. See ?select for details. 
## Bimap interface: 
x <- org.Hs.egSYMBOL 
# Get the gene symbol that are mapped to an entrez gene identifiers 
mapped_genes <- mappedkeys(x) 
# Convert to a list 
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) { 
  # Get the SYMBOL for the first five genes
  xx[1:5]
  # Get the first one 
  xx[[1]] } 
#choose the column to be converted
a<- genes$GENE
a = as.character(a)
a

d= getSYMBOL(na.omit(a), data='org.Hs.eg')
d

#convert p value to z score
pval = genes$P
Zscore=abs(scale(qnorm(pval/2)))

#Check if converted correctly
Zscore

NodeFile = data.frame(d,Zscore)
write.table(NodeFile, file='NodeFile.txt', append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, quote = FALSE)
