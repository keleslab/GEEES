edge <- function(peak_data,gene_data,num,whole=F,phat=NULL,tfidf = F){
  h = floor(num/2)
  if(tfidf == F){peak_data <- as.numeric(peak_data>0)
  sig <- sapply(seq_along(gene_data),FUN = function(i){
    nx <- (rank(abs(gene_data-gene_data[i]),ties.method = "random")<=num)
    if(is.null(phat)){
      phat <- mean(peak_data==1)
    }
    nxy <- sum(nx&peak_data==1)
    #print(c(nxy,sum(nx),num))
    if(gene_data[i]>=mean(gene_data)){
      p <- binom.test(x = nxy,n = sum(nx),p = phat,alternative = "greater")$p.value
    }
    else{
      p <- binom.test(x = nxy,n = sum(nx),p = phat,alternative = "less")$p.value
    }
    if(whole){
      p <- ifelse(peak_data[i]==1,p,1)
    }
    
    return(p)
  })}
  else{
    n2 <- length(peak_data)
    sig <- sapply(seq_along(gene_data),FUN = function(i){
      nx <- find_neighbor(gene_data,h,i)
      ny <- find_neighbor(peak_data,h,i)
      numx <- sum(nx)
      numy <- sum(ny)
      numxy <- sum(nx&ny)
      d <- (numxy*n2-numx*numy)/(sqrt(numx*numy)*sqrt((n2-numx)*(n2-numy))/sqrt(n2-1)+1e-6)
      p <- pnorm(d,lower.tail = F)
      return(p)
    })
  }
  
  return(sig)
}
find_neighbor <- function(data,h,i){
  s = sum(data==data[i])
  if(s>=h){
    return(data==data[i])
  }
  else{
    n2 <- length(data)
    n3 <- sum(data==0)
    r <- rank(data,ties.method = "random")
    m <- min(r[data==data[i]])
    u <- data[which(r==min(s+m+h,n2))]
    l <- data[which(r==max(n3*(n3>h)+1,m-h))]
    return(data<=u&data>=l)
  }
  
}
NetworkInCell <- function(
    object,
    peak.assay,
    expression.assay,
    expression.slot = "data",
    peak.slot = "counts",
    gene.coords = NULL,
    distance = 5e+05,
    min.distance = NULL,
    min.cells = 1,
    genes.use = NULL,
    n_sample = 100,
    pvalue_cutoff = 1.1,
    score_cutoff = -0.1,
    verbose = TRUE,
    promoter.data = NULL,
    cl,
    whole = F,
    beta=F,
    AB=NULL,
    num,
    tfidf = F
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
  
  if (is.null(x = gene.coords)) {
    gene.coords <- CollapseToLongestTranscript(
      ranges = Annotation(object = object[[peak.assay]])
    )
  }
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
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  
  # run in parallel across genes
  l <- sapply(1:length(genes.use),list)
  clusterExport(cl,varlist = c("peak_distance_matrix","genes.use","expression.data","peak.data","edge","AB","beta","whole","peaks","num","tfidf","find_neighbor"),envir=environment())
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
        if(beta){
          xi <- sum(peak_data>0)
          phat <- (xi+AB[1,peak.use[j]])/(AB[1,peak.use[j]]+AB[2,peak.use[j]]+length(peak_data))
          if(AB[1,peak.use[j]]==0){
            phat=mean(peak_data>0)
          }
          return(edge(peak_data,gene_data,whole = whole,phat = phat,num = num,tfidf = tfidf))
        }
        return(edge(peak_data,gene_data,whole = whole,num = num,tfidf = tfidf))
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
  
  names(tmp) <- genes.use
  return(tmp)
  
}
SNframe <- function(data,logname=NULL,netname=NULL,genes.use,whole = F,nclust=10,print=F,beta=F,AB=NULL,num=30,peak.slot = "counts",expression.slot = "counts",tfidf = F,peak.assay = "ATAC",expression.assay = "SCT"){
  cl <- makeCluster(nclust,outfile = logname)
  SpecificNetwork <- NetworkInCell(data,peak.assay = peak.assay,expression.assay=expression.assay,peak.slot = peak.slot,expression.slot = expression.slot,genes.use = genes.use,cl = cl,whole = whole,beta = beta,AB=AB,num=num,tfidf = tfidf)
  #SpecificNetwork <- NetworkInCell(data,"ATAC","SCT",genes.use = genes.use,whole = whole,beta = beta,type = type,peak_for_beta = peak_for_beta,celltypes = celltypes)
  stopCluster(cl)
  if(print){
    return(SpecificNetwork)
  }
  save(SpecificNetwork,file = netname)
  # cl <- makeCluster(2,outfile = "cd4toplog.txt")
  # system.time(
  # cd4SpecificNetwork <- NetworkInCell(cd4,"ATAC","SCT",genes.use = top_genes[1:5])
  # )
  # stopCluster(cl)
  # promoter.pbmc.marker <- count.promoter(data,1000,500,names(SpecificNetwork))
  # SpecificNetPromoter <- promoter.p(data@assays$SCT@data,promoter.pbmc.marker,names(SpecificNetwork))
  # 
  # save(SpecificNetPromoter,file = promname)
}
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
    fishers <- sapply(1:ncol(i),FUN = function(k) return(fisher(i[,k])$p))
    stat.fisher <- p.adjust(fishers,method = "BH")
    stat.median <- colMedians(-log10(i))
    stat.mean <- colMeans(-i)
    stat <- colMeans(i<0.05)
    peak <- colnames(i)
    gene <- rep(names(data)[j],length(peak))
    return(data.frame(gene,peak,stat,stat.fisher,stat.median,stat.mean))
    
  })
  cd4_gene_cis <- do.call("rbind",cd4_gene_cis)
  return(cd4_gene_cis)
}
