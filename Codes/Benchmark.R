Benchmark <- function(bench,data,bench_name,not_pair,badsep = c(":","-"),goodsep = c(":","-")){
  # not_pair: all peaks involved in evaluation (list)
  for(i in 1:length(bench)){
    good <- bench[[i]]
    bad <- not_pair[[i]]
    data <- cbind(data,TF = NA)
    bad_peak <- find_peak(StringToGRanges(data$peak),StringToGRanges(bad,sep = badsep))
    bad_peak <- rownames(bad_peak)[rowSums(bad_peak)>0]
    data$TF[is.element(data$peak,bad_peak)] <- F
    data.list <- pblapply(unique(data$gene),FUN = function(k) return(data[data$gene==k,]))
    names(data.list) <- unique(data$gene)
    good_gene <- intersect(data$gene,good$gene)
    
    data.list[good_gene] <- pblapply(good_gene, FUN = function(j){
      subset <- data.list[[j]]
      peaks <- subset$peak
      # if(negative){
      #   bad <- good[good$Significant=="False",]
      #   good <- good[good$Significant=="True",]
      #   bad_peak <- bad$peak[bad$gene==j]
      #   bad_peak <- find(StringToGRanges(peaks),StringToGRanges(bad_peak,sep=badsep))
      #   bad_peak <- rownames(bad_peak)[rowSums(bad_peak)>0]
      #   subset$TF[is.element(subset$peak,bad_peak)] <- F
      # }
      good_peak <- good$peak[good$gene==j]
      good_peak <- find_peak(StringToGRanges(peaks),StringToGRanges(good_peak,sep = goodsep))
      good_peak <- rownames(good_peak)[rowSums(good_peak)>0]
      subset$TF[is.element(subset$peak,good_peak)] <- T
      
      return(subset)
    })
    data.list <- do.call("rbind",data.list)
    colnames(data.list) <- c(colnames(data)[1:(ncol(data)-1)],bench_name[i])
    data <- data.list
    
  }
  
  return(data)
}