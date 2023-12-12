#' Generate adaptive regression results for gene-enhancer association
#'
#' @param object a data frame has 3 columns: gene (gene name), peak (matched cis-peak), coaccess (coaccessibility calculated by cicero between the peak and the promoter of the gene)
#' @param rep3_data a seurat object
#' @param PFER upper bound per family error rate
#' @param q number of (unique) selected variables (or groups of variables depending on the model) that are selected on each subsample.
#' @param peak.assay the name of the assay containing enhancer accessibility
#' @param expression.assay the name of the assay containing gene expression
#' @param expression.slot the name of the slot to use in expression assay
#' @param peak.slot the name of the slot to use in enhancer accessibility assay
#'
#' @return a data frame with adaptive regression results: gene (gene name), peak (paired cis-peak to this gene), max_rate (adaptive regression result for the gene-peak pair)
#'
#' @importFrom SeuratObject GetAssayData
#' @importFrom pbapply pblapply
#' @importFrom stabs stabsel glmnet.lasso
#'
#' @export

adaptive.regress.stabs <- function(
    object,
    rep3_data,
    PFER,
    q,
    peak.assay = "ATAC",
    expression.assay = "SCT",
    expression.slot = "data",
    peak.slot = "counts"
){
  gene <- unique(object$gene)
  expression.data <- GetAssayData(
    object = rep3_data, assay = expression.assay, slot = expression.slot
  )
  gene <- gene[rowSums(expression.data[gene,])>0]
  peak.data <- GetAssayData(
    object = rep3_data, assay = peak.assay, slot = peak.slot
  )
  object$penalty <- 1/(object$coaccess-min(object$coaccess)+0.5)-min(1/(object$coaccess-min(object$coaccess)+0.5))+0.01
  tmp <- pblapply(gene,FUN = function(i){
    coord <- which(object$gene==i)
    if(length(coord)<5){
      return(vector(mode = "list",length = 6))
    }
    #print(i)
    subset <- object[coord,]
    peak <- subset$peak
    clust_count <- t(peak.data[peak, ,drop = F])
    penalty <- subset$penalty
    poisson.fit.myc <- stabsel(x = as.matrix(clust_count),y = as.matrix(expression.data[i,]),fitfun = glmnet.lasso,PFER = PFER,q = q, args.fitfun = list(penalty.factor = penalty),sampling.type = "MB" )
    max_rate <- poisson.fit.myc$max
    max_rate_peak <- max_rate[subset$peak]
    selected <- rep(F,length(subset$peak))
    names(selected) <- subset$peak
    selected[names(poisson.fit.myc$selected)] <- T
    PFER <- rep(poisson.fit.myc$PFER,length(subset$peak))
    cutoff <- rep(poisson.fit.myc$cutoff,length(subset$peak))
    #result <- list(gene  = subset$gene,peak = subset$peak,max_rate = max_rate_peak,selected = selected,PFER = PFER,cutoff = cutoff)
    result <- list(gene  = subset$gene,peak = subset$peak,max_rate = max_rate_peak)
    return(result)
  })
  gene.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 1))
  peak.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 2))
  max.rate.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 3))
  # selected.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 4))
  # PFER <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 5))
  # cutoff <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 6))
  #links <- data.frame(gene = gene.vec,peak = peak.vec,max_rate = max.rate.vec,selected = selected.vec,PFER = PFER,cutoff = cutoff)
  links <- data.frame(gene = gene.vec,peak = peak.vec,max_rate = max.rate.vec)
  return(links)
}


#' Generate multi-response regression results for gene-enhancer association
#'
#' @param object a data frame has 2 columns: gene (gene name), peak (matched cis-peak)
#' @param rep3 a seurat object
#' @param PFER upper bound per family error rate
#' @param q number of (unique) selected variables (or groups of variables depending on the model) that are selected on each subsample.
#' @param peak.assay the name of the assay containing enhancer accessibility
#' @param expression.assay the name of the assay containing gene expression
#' @param expression.slot the name of the slot to use in expression assay
#' @param peak.slot the name of the slot to use in enhancer accessibility assay
#'
#' @return a data frame with adaptive regression results: gene (gene name), peak (paired cis-peak to this gene), max_rate (multi-response regression result for the gene-peak pair)
#'
#' @importFrom SeuratObject GetAssayData
#' @importFrom pbapply pblapply
#' @importFrom stabs stabsel
#'
#' @export
remMap.regression.stabs <- function(
    object,
    rep3,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    expression.slot = "data",
    peak.slot = "counts",
    PFER,
    q
){
  gene <- unique(object$gene)
  expression.data <- GetAssayData(
    object = rep3, assay = expression.assay, slot = expression.slot
  )
  peak.data <- GetAssayData(
    object = rep3, assay = peak.assay, slot = peak.slot
  )
  gene <- gene[rowSums(expression.data[gene,])>0]
  tmp <- pblapply(gene,FUN = function(i){
    # print(i)
    coord <- which(object$gene==i)
    if(length(coord)<5){
      return(list(gene  = NULL,peak = NULL,mean_rate_gene = NULL,max_rate_gene = NULL,mean_rate_prom = NULL,max_rate_prom = NULL))
    }
    subset <- object[coord,]
    peak <- subset$peak
    clust_count <- t(peak.data[peak, ,drop = F])
    promo.peak <- gsub("_","-",subset[1,subset$ind[1]])
    promot <- peak.data[promo.peak,]
    poisson.fit.myc <- stabsel(x = as.matrix(clust_count),y = as.matrix(cbind(expression.data[i,],promot)),fitfun = rem.function,q = q,PFER = PFER, sampling.type = "MB" )
    max_rate <- poisson.fit.myc$max
    max_rate_peak <- max_rate[subset$peak]
    PFER <- rep(poisson.fit.myc$PFER,length(subset$peak))
    selected <- rep(F,length(subset$peak))
    names(selected) <- subset$peak
    selected[names(poisson.fit.myc$selected)] <- T
    cutoff <- rep(poisson.fit.myc$cutoff,length(subset$peak))
    result <- list(gene  = subset$gene,peak = subset$peak,max_rate = max_rate_peak,PFER = PFER,selected = selected,cutoff = cutoff)
    return(result)
  })
  gene.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 1))
  peak.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 2))
  max.rate.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 3))
  PFER <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 4))
  selected <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 5))
  cutoff <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 6))
  #links <- data.frame(gene = gene.vec,peak = peak.vec,max_rate = max.rate.vec,PFER = PFER,selected = selected,cutoff = cutoff)
  links <- data.frame(gene = gene.vec,peak = peak.vec,max_rate = max.rate.vec)
  return(links)
}

#' multi-response regression for one gene
#' @param x peak accessibility matrix
#' @param y gene expression vector
#' @param q the number of selected peaks
#' @param ... other arguments can be passed
#'
#' @return a list of a logical vector indicating which peaks are selected at this round of stability selection
#' @import remMap
#'

rem.function <- function(x,y,q,...){
  lambda <- max(abs(t(as.matrix(x))%*%y))
  if(lambda==0){
    ret <- logical(ncol(x))
    sequence <- matrix(logical(ncol(x)*21),ncol = 21)
    names(ret) <- colnames(x)
    return(list(selected = ret))
  }
  lamL1.v = 10^seq((log10(lambda)-2),log10(lambda),0.1) # from small to large until q selected
  index <- c()
  for(i in 1:21){
    index <- c(index,lapply(1:21,FUN = function(j) return(list(x1 = i, x2 = j))))
  }

  coef <- lapply(index,FUN = function(j){
    lam1 <- lamL1.v[j$x1]
    lam2 <- lamL1.v[j$x2]
    fit <- remMap(X.m = x,Y.m = y,lamL1 = lam1,lamL2 = lam2)
    a <- fit$phi!=0
    return(a)
  })
  gene.coef <- lapply(coef, FUN = function(i){
    return(i[,1])
  })
  selected.num <- sapply(gene.coef,FUN = function(i) return(sum(i)))
  good.num <- which(selected.num<=q)
  good.num.max <- good.num[selected.num[good.num]==max(selected.num[good.num])]
  ranks <- order(sapply(index[good.num.max],FUN = function(i) return(sum(lamL1.v[c(i$x1,i$x2)]))))
  selected <- gene.coef[good.num.max[ranks[1]]][[1]]
  ret <- logical(ncol(x))
  ret[selected] <- T
  names(ret) <- colnames(x)
  # sequence <- gene.coef[good.num]
  # sequence <- matrix(unlist(sequence),nrow = ncol(x))
  # rownames(sequence) <- colnames(x)
  # colnames(sequence) <- paste0("s",seq(ncol(sequence)))
  return(list(selected = ret))
}


#' Generate sequential regression results for gene-enhancer association
#'
#' @param object a data frame has 5 columns: Peak1 (the first peak as input for cicero), Peak2 (the second peak as input for cicero), gene (gene name), peak (matched cis-peak), ind (indicator to be 1 or 2 indicating which peak before is the promoter for the gene)
#' @param rep3_data a seurat object
#' @param PFER upper bound per family error rate
#' @param q number of (unique) selected variables (or groups of variables depending on the model) that are selected on each subsample.
#' @param peak.assay the name of the assay containing enhancer accessibility
#' @param expression.assay the name of the assay containing gene expression
#' @param expression.slot the name of the slot to use in expression assay
#' @param peak.slot the name of the slot to use in enhancer accessibility assay
#'
#' @return a data frame with adaptive regression results: gene (gene name), peak (paired cis-peak to this gene), max_rate (sequential regression result for the gene-peak pair)
#'
#' @importFrom SeuratObject GetAssayData
#' @importFrom pbapply pblapply
#' @importFrom stabs stabsel
#'
#' @export
sequential.regression.stabs <- function(
    object,
    rep3_data,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    expression.slot = "data",
    peak.slot = "counts",
    PFER,
    q
){
  gene <- unique(object$gene)
  expression.data <- GetAssayData(
    object = rep3_data, assay = expression.assay, slot = expression.slot
  )
  peak.data <- GetAssayData(
    object = rep3_data, assay = peak.assay, slot = peak.slot
  )
  gene <- gene[rowSums(expression.data[gene,])>0]
  tmp <- pblapply(gene,FUN = function(i){
    # print(i)
    coord <- which(object$gene==i)
    if(length(coord)<5){
      return(list(gene  = NULL,peak = NULL,mean_rate = NULL,max_rate = NULL,count = NULL,curoff=NULL))
    }
    subset <- object[coord,]
    peak <- subset$peak
    clust_count <- t(peak.data[peak, ,drop = F])
    promo.peak <- gsub("_","-",subset[1,subset$ind[1]])
    promot <- peak.data[promo.peak,]

    poisson.fit.myc <- stabsel(x = as.matrix(clust_count),y = as.matrix(cbind(promot,expression.data[i,])),fitfun = sequential.function,q = q,PFER = PFER, sampling.type = "MB" )
    max_rate <- poisson.fit.myc$max
    max_rate_peak <- max_rate[subset$peak]
    PFER <- rep(poisson.fit.myc$PFER,length(subset$peak))
    selected <- rep(F,length(subset$peak))
    names(selected) <- subset$peak
    selected[names(poisson.fit.myc$selected)] <- T
    cutoff <- rep(poisson.fit.myc$cutoff,length(subset$peak))
    result <- list(gene  = subset$gene,peak = subset$peak,max_rate = max_rate_peak,PFER = PFER,selected = selected,cutoff = cutoff)
    return(result)
  })
  gene.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 1))
  peak.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 2))
  max.rate.vec <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 3))
  PFER <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 4))
  selected <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 5))
  cutoff <- do.call(what = "c", args = lapply(X = tmp, FUN = `[[`, 6))
  links <- data.frame(gene = gene.vec,peak = peak.vec,max_rate = max.rate.vec)
  return(links)
  #return(tmp)
}

#' sequential regression for one gene
#' @param x peak accessibility matrix
#' @param y gene expression vector
#' @param q the number of selected peaks
#'
#' @return a list of a logical vector indicating which peaks are selected at this round of stability selection
#'
#' @importFrom glmnet glmnet predict.glmnet coef.glmnet
#' @importFrom stats predict coef
sequential.function <- function(x,y,q){
  poisson.fit.myc <- glmnet(y = y[,1], x = x,family = "poisson",alpha = 1)
  if(poisson.fit.myc$lambda[1]==Inf){
    ret <- logical(ncol(x))
    names(ret) <- colnames(x)
    return(list(selected = ret))
  }
  positive.peak <- rownames(coef(poisson.fit.myc)[-1, ,drop = F])[rowSums(coef(poisson.fit.myc)[-1, ,drop = F]>0)>rowSums(coef(poisson.fit.myc)[-1, ,drop = F]<=0)]
  if(length(positive.peak)<2){
    ret <- logical(ncol(x))
    names(ret) <- colnames(x)
    ret[positive.peak] <- T
    return(list(selected = ret))
  }
  poisson.fit.pro <- glmnet(y = y[,-1], x = x[,positive.peak],alpha = 1,pmax = q)
  if(poisson.fit.pro$lambda[1]==Inf){
    ret <- logical(ncol(x))
    names(ret) <- colnames(x)
    return(list(selected = ret))
  }
  selected <- predict(poisson.fit.pro, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  names(ret) <- colnames(x)
  ret[positive.peak[selected]] <- T
  return(list(selected = ret))
}
