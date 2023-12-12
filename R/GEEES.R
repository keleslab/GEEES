#' Compute cell-specific GEEES statistics for a specific gene-enhancer pair
#'
#' @param peak_data a peak accessibility vector.
#' @param gene_data a gene expression vector.
#' @param num size of the cell neighbor.
#'
#' @return a vector containing GEEES statistics for gene-enhancer pair on each single cell
#' @importFrom stats binom.test
#'
#' @export
edge <- function(peak_data,gene_data,num){
  h = floor(num/2)
  peak_data <- as.numeric(peak_data>0)
  sig <- sapply(seq_along(gene_data),FUN = function(i){
    nx <- (rank(abs(gene_data-gene_data[i]),ties.method = "random")<=num)
    #if(is.null(phat)){
      phat <- mean(peak_data==1)
    #}
    nxy <- sum(nx&peak_data==1)
    #print(c(nxy,sum(nx),num))
    if(gene_data[i]>=mean(gene_data)){
      p <- binom.test(x = nxy,n = sum(nx),p = phat,alternative = "greater")$p.value
    }
    else{
      p <- binom.test(x = nxy,n = sum(nx),p = phat,alternative = "less")$p.value
    }

    return(p)
  })


  return(sig)
}

# find_neighbor <- function(data,h,i){
#   s = sum(data==data[i])
#   if(s>=h){
#     return(data==data[i])
#   }
#   else{
#     n2 <- length(data)
#     n3 <- sum(data==0)
#     r <- rank(data,ties.method = "random")
#     m <- min(r[data==data[i]])
#     u <- data[which(r==min(s+m+h,n2))]
#     l <- data[which(r==max(n3*(n3>h)+1,m-h))]
#     return(data<=u&data>=l)
#   }
#
# }
#' Compute cell-specific GEEES statistics for all specified gene with their cis-peaks
#'
#' @param object a Seurat object
#' @param peak.assay the name of the assay containing enhancer accessibility
#' @param expression.assay the name of the assay containing gene expression
#' @param expression.slot the name of the slot to use in expression assay
#' @param peak.slot the name of the slot to use in enhancer accessibility assay
#' @param distance the window size of peak considered as cis-peak for genes
#' @param min.distance the minimum distance between gene and peak for a peak to be considered as a candidate enhancer
#' @param min.cells gene and peak filtering where only genes and peaks that are expressed/accessible in at least min.cells cells are considered
#' @param genes.use a vector of gene names for enhancer detection (all genes in the object will be used if not given)
#' @param verbose if the total numbers of genes and peaks are reported
#' @param cl a list of clusters for parallel computing
#' @param num number of cells in a cell neighbor
#'
#' @return a list of matrices where each element is a cell by peak matrix containing the cell-specific GEEES statistics for a specific gene and its cis-peaks
#'
#' @import Seurat
#' @importFrom Signac Annotation FindTopFeatures
#' @importFrom SeuratObject GetAssayData
#' @importFrom GenomicRanges granges
#' @importFrom parallel clusterEvalQ clusterExport parLapply
#' @import Matrix
#' @importFrom pbapply pblapply
#'
#' @export

NetworkInCell <- function(
    object,
    peak.assay,
    expression.assay,
    expression.slot = "data",
    peak.slot = "counts",
    distance = 5e+05,
    min.distance = NULL,
    min.cells = 1,
    genes.use = NULL,
    verbose = T,
    cl,
    num
) {
  if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) {
      stop("min.distance should be a numeric value")
    }
    if (min.distance < 0) {
      warning("Requested a negative min.distance value, setting min.distance to zero")
      min.distance <- NULL
    } else if (min.distance == 0) {
      min.distance <- NULL
    }
  }


    gene.coords <- CollapseToLongestTranscript(
      ranges = Annotation(object = object[[peak.assay]])
    )

  meta.features <- GetAssayData(
    object = object, assay = peak.assay, slot = "meta.features"
  )

  peak.data <- GetAssayData(
    object = object, assay = peak.assay, slot = peak.slot
  )
  if (!("count" %in% colnames(x = meta.features))) {
    # compute total count
    hvf.info <- FindTopFeatures(object = peak.data)
    hvf.info <- hvf.info[rownames(x = meta.features), ]
    meta.features <- cbind(meta.features, hvf.info)
  }
  expression.data <- GetAssayData(
    object = object, assay = expression.assay, slot = expression.slot
  )
  peakcounts <- meta.features[rownames(x = peak.data), "count"]
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  if (is.null(x = genes.use)) {
    expression.data <- expression.data[genes.keep, ]
  } else {
    genes.keep <- intersect(
      x = names(x = genes.keep[genes.keep]), y = genes.use
    )
    expression.data <- expression.data[genes.keep, ]
  }
  if (verbose) {
    message(
      "Testing ",
      nrow(x = expression.data),
      " genes and ",
      sum(peaks.keep),
      " peaks"
    )
  }
  genes <- rownames(x = expression.data)
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
  peaks <- granges(x = object[[peak.assay]])
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- DistanceToTSS(
    peaks = peaks,
    genes = gene.coords.use,
    distance = distance
  )
  if (!is.null(x = min.distance)) {
    peak_distance_matrix_min <- DistanceToTSS(
      peaks = peaks,
      genes = gene.coords.use,
      distance = min.distance
    )
    peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
  }
  if (sum(peak_distance_matrix) == 0) {
    stop("No peaks fall within distance threshold\n",
         "Have you set the proper genome and seqlevelsStyle for ",
         peak.assay,
         " assay?")
  }
  genes.use <- colnames(x = peak_distance_matrix)[colSums(peak_distance_matrix)>0]
  all.peaks <- rownames(x = peak.data)
  print(sum(rowSums(peak_distance_matrix)>0))

  peak.data <- t(x = peak.data)
  peaks <- rownames(peak_distance_matrix)
  # if (nbrOfWorkers() > 1) {
  #   mylapply <- future_lapply
  # } else {
  #   mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  # }

  # run in parallel across genes
  if(!is.null(cl)){
  l <- sapply(1:length(genes.use),list)
  clusterExport(cl,varlist = c("peak_distance_matrix","genes.use","expression.data","peak.data","edge","peaks","num","find_neighbor"),envir=environment())
  # #
  clusterEvalQ(cl,library(Matrix))
  #clusterEvalQ(cl,library(glmnet))


  tmp <- parLapply(
    cl,
    l,
    #tmp <- pblapply(l,
    fun = function(i) {
      if(i %% 30 == 0) print(paste("echo 'now processing:",i,"'"))
      #print(paste("echo 'now processing:",i,"'"))
      #peaks <- rownames(peak_distance_matrix)
      # print(str(peaks))
      # print(str(peak_distance_matrix))
      ind <- as.logical(peak_distance_matrix[, genes.use[[i]]])
      # print(str(ind))
      peak.use <- peaks[ind]
      # print(genes.use[[i]])
      # print(str(peak.use))
      # print(ind)

      gene.expression <- expression.data[genes.use[[i]], , drop = FALSE]
      #print(i)
      peaks_result <- sapply(seq_along(peak.use),FUN = function(j){
        peak_data <- peak.data[,peak.use[j]]

        gene_data <- as.vector(gene.expression)
        #print(j)

        return(edge(peak_data,gene_data,num = num))
      })
      if(length(peak.use)==1){
        peaks_result <- matrix(peaks_result,ncol = 1)
      }
      else if(length(peak.use)==0){
        return(NULL)
      }
      colnames(peaks_result) <- peak.use
      #print(str(peaks_result))
      return(peaks_result)

    }
  )
  }
  else{
    l <- sapply(1:length(genes.use),list)
    tmp <- pblapply(l,
      FUN = function(i) {
        ind <- as.logical(peak_distance_matrix[, genes.use[[i]]])
        peak.use <- peaks[ind]

        gene.expression <- expression.data[genes.use[[i]], , drop = FALSE]
        peaks_result <- sapply(seq_along(peak.use),FUN = function(j){
          peak_data <- peak.data[,peak.use[j]]

          gene_data <- as.vector(gene.expression)

          return(edge(peak_data,gene_data,num = num))
        })
        if(length(peak.use)==1){
          peaks_result <- matrix(peaks_result,ncol = 1)
        }
        else if(length(peak.use)==0){
          return(NULL)
        }
        colnames(peaks_result) <- peak.use
        return(peaks_result)

      }
    )
  }

  names(tmp) <- genes.use
  return(tmp)

}
#' Run GEEES in a parallel setting with log file
#'
#' @param data a seurat object
#' @param logname the path of the log file
#' @param netname the path of the file saving results
#' @param genes.use a vector of genes to be used
#' @param nclust number of clusters to be used in the parallel computation
#' @param print if the result should be returned or saved
#' @param num number of cells in the neighborhood
#' @param peak.slot the name of the slot to be used in accessibility assay
#' @param expression.slot the name of the slot to be used in expression assay
#' @param peak.assay the name of the assay containing enhancer accessibility
#' @param expression.assay the name of the assay containing gene expression
#'
#' @return a list of matrices where each element is a cell by peak matrix containing the cell-specific GEEES statistics for a specific gene and its cis-peaks
#'
#' @importFrom parallel makeCluster stopCluster
#' @export
GEEESframe <- function(data,logname=NULL,netname=NULL,genes.use,nclust=10,print=F,num=30,peak.slot = "counts",expression.slot = "counts",peak.assay = "ATAC",expression.assay = "SCT"){
  cl <- makeCluster(nclust,outfile = logname)
  SpecificNetwork <- NetworkInCell(data,peak.assay = peak.assay,expression.assay=expression.assay,peak.slot = peak.slot,expression.slot = expression.slot,genes.use = genes.use,cl = cl,num=num)
  #SpecificNetwork <- NetworkInCell(data,"ATAC","SCT",genes.use = genes.use,whole = whole,beta = beta,type = type,peak_for_beta = peak_for_beta,celltypes = celltypes)
  stopCluster(cl)
  if(print){
    return(SpecificNetwork)
  }
  save(SpecificNetwork,file = netname)
}

#' Aggregating cell-specific GEEES statistic over the population for a summarized statistic
#' @param data The data frame of cell-specific gene-enhancer GEEES statistics returned by NetworkInCell
#'
#' @return a data frame with aggregated GEEES statistics for all gene-enhancer pairs
#'
#' @importFrom pbapply pblapply
#' @importFrom MatrixGenerics colMedians
#' @importFrom BiocGenerics do.call
#'
#' @return a data frame with 3 columns: gene (gene name), peak (the paired cis-peak to the gene), max_rate (aggregated GEEES statistic for this pair)
#' @export
generate_pair <- function(data){
  # generate rank in different criteria
  cd4_gene_cis <- pblapply(1:length(data),FUN = function(j){
    i <- data[[j]]
    if(is.null(i)) return(NULL)
    # if(is.null(dim(i))){
    #   stat <- mean(i)
    #   peak <- names(i)
    #   print(j)
    #   gene <- names(cd4SpecificNetwork)[j]
    #   return(data.frame(gene,peak,stat))
    # }
    i[is.na(i)] <- 1

    #adj <- sapply(1:nrow(i),FUN = function(k) return(p.adjust(i[k,],method = "BH")))
    #adj <- t(adj)
    #stat.BH <- colMeans(adj<0.05)
    # fishers <- sapply(1:ncol(i),FUN = function(k) return(fisher(i[,k])$p))
    # stat.fisher <- p.adjust(fishers,method = "BH")
    stat.median <- colMedians(-log10(i),useNames = T)
    # stat.mean <- colMeans(-i)
    # stat <- colMeans(i<0.05)
    peak <- colnames(i)
    gene <- rep(names(data)[j],length(peak))
    return(data.frame(gene,peak,max_rate=stat.median))

  })
  cd4_gene_cis <- do.call("rbind",cd4_gene_cis)
  return(cd4_gene_cis)
}
