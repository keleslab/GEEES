#' Find overlap between two lists of chromatin regions
#'
#' @param peaks first list of chromatin regions
#' @param genes second list of chromatin regions
#' @param sep the vector of separator to use for genome string
#'
#' @importFrom Signac GRangesToString
#' @importFrom IRanges findOverlaps
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methods as
#' @importFrom stats end start
#' @return  a binary matrix indicating if two pairs of peaks are overlapped with first list on the row and second list on the column
#' @keywords internal
find_peak <- function(peaks, genes, sep = c("-", "-")) {
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- GRangesToString(grange = genes, sep = sep)
  return(hit_matrix)
}

#' @title CollapseToLongestTranscript
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom data.table as.data.table
#'
#' @export
#'
CollapseToLongestTranscript <- function(ranges) {
  seqnames=start=end=strand=gene_biotype=gene_name=NULL
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , list(unique(seqnames),
           min(start),
           max(end),
           strand[[1]],
           gene_biotype[[1]],
           gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

#' @importFrom Signac findOverlaps Extend GRangesToString
#' @importFrom IRanges resize
#' @importFrom Matrix sparseMatrix
DistanceToTSS <- function(peaks, genes, distance = 200000, sep = c("-", "-")) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}

#' @importFrom IRanges resize
#' @importFrom Signac StringToGRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom methods as
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb seqnames
#' @keywords internal
LinksToGRanges <- function(linkmat, gene.coords, sep = c("-", "-")) {
  # get TSS for each gene
  tss <- resize(gene.coords, width = 1, fix = 'start')
  gene.idx <- sapply(
    X = rownames(x = linkmat),
    FUN = function(x) {
      which(x = x == tss$gene_name)[[1]]
    }
  )
  tss <- tss[gene.idx]

  # get midpoint of each peak
  peak.ranges <- StringToGRanges(
    regions = colnames(x = linkmat),
    sep = sep
  )
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2)

  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "dgTMatrix")

  # create dataframe
  df <- data.frame(
    chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
    tss = start(x = tss)[dgtm@i + 1],
    pk = midpoints[dgtm@j + 1],
    score = dgtm@x,
    gene = rownames(x = linkmat)[dgtm@i + 1],
    peak = colnames(x = linkmat)[dgtm@j + 1]
  )

  # work out start and end coords
  df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
  df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
  df$tss <- NULL
  df$pk <- NULL

  # convert to granges
  gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  return(sort(x = gr.use))
}

#' Distance adjusted statistics
#' 
#' @param result result statistics without distance adjustment
#' @param gene.coords data frame for gene coordinates
#' 
#' @return data frame with one additional column named distance_adjusted for distance adjusted statistics
#' 
#' @importFrom Signac StringToGRanges
#' @export
Distance_adjust <- function(result,gene.coords){
  result$TSS <- gene.coords[result$gene,"start"]
  result$peak_start <- data.frame(StringToGRanges(result$peak))$start
  result$dist <- -abs(result$TSS-result$peak_start)
  result$distance_adjusted <- (result$max_rate+0.05)*exp(result$dist/200000)
  return(result)
}