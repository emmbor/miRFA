# miRFA

MiRFA is a pipeline built for functional analysis of a mature miRNA in a pancreatic cancer context. The pipeline automatically performs the following steps:

1) MiRNA target prediction


2) Correlation analyses between miRNA-target genes on mRNA and protein expression levels

3) Functional enrichment of significantly correlated targets on mRNA and protein expression levels

## Getting started

``` R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi",
                       "MODA", 
                       "STRINGdb", 
                       "limma", 
                       "org.Hs.eg.db"), 
                     version = "3.8")
```R

install.packages(c("foreach", 
                   "doParallel", 
                   "Rcpp", 
                   "dynamicTreeCut", 
                   "flashClust", 
                   "reticulate", 
                   "plyr", 
                   "parallel", 
                   "igraph", 
                   "WGCNA"))

### Installing

Clone or download the repository

## Running the pipeline

### Input
A demo file shows the structure of an accepted input. Notice that the input miRNA has to be a mature miRNA name, e.g. 'hsa-miR-885-5p'. 

### Change parameters

The file demorun.R shows how to call the mirna_search.R script. The following parameters should be adjusted:

```{r}
path = "path/to/data"

mirnas<- read.table("your_mirna_list.txt",header=F)

## Optional - add a name for your run
runname<-name_run(mirnas,name_of_run="Run1")
```
Now the pipeline can be run.

The script is applied to every row on the miRNA input list. Another option is to only run one miRNA, e.g.:

```{r}
mirna_search(miRNA="hsa-miR-885-5p")
```

### Output

The following lists are generated from the pipeline: 

1. Tarbase miRNA target genes
2. TargetScan miRNA target genes
3. DIANA-microT-CDS miRNA target genes
4. Union of miRNA target genes (1-3)
5. A Venn Diagram showing the overlap of miRNA target genes generated from the target prediction databases
6. Positive, significant correlations on mRNA expression levels between the miRNA and its predicted target. 
7. Negative, significant correlations on mRNA expression levels between the miRNA and its predicted target. 
8. Positive, significant correlations on protein expression levels between the miRNA and its predicted target. 
9. Negative, significant correlations on protein expression levels between the miRNA and its predicted target. 
10. Significant correlations on both mRNA and protein expression levels between the miRNA and its predicted target. 
11. Significantly enriched GO terms.
12. Significantly enriched KEGG pathways. 














