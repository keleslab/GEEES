# GEEES
## Method Overview
**GEEES** is a novelly proposed approach for inferring cell-specific **G**ene-**E**nhanc**E**r Int**E**ractions from Multi-modal **S**ingle Cell Data with transcriptome and chromatin accessibility profiles. GEEES estimates gene-enhancer associations at the single
cell level by considering a cell neighbourhood defined by both the expression of the gene and the accessibility of the enhancer in the gene-enhancer pair.

![alt text](https://github.com/Shuyang12138/GEEES/blob/main/Figures/GEEESFramework.jpg?raw=true)

The codes implementing GEEES are provided in `R/Functions/GEEES.R` which could be used as below.
```r
# Input:
# object: a seurat object with normalized RNA and ATAC assay
# gene.use: a list of gene name for detecting regulatory enhancers
# nclust: the number of clusters used for parallel computing
GEEES_result <- SNframe(object,"Bench_log.txt","benchSN_K562_30.Rdata",genes.use = gene.use,nclust = nclust)
```

## Benchmark Study
GEEES is benchmarked against the state-of-the-art methods and a number of multivariate regression approaches we devised on a diverse set of multi-modal single cell datasets using a wide variaty of gold standard gene-enhancer interaction datasets. The datasets involved and summarized result table is shown as below.

![alt text](https://github.com/Shuyang12138/GEEES/blob/main/Figures/Benchmark.jpg?raw=true)

The codes implementing the multivariate regression approaches are provided in `R/Functions/Regression.R`
