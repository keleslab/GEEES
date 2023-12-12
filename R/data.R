#' A list containing normalized gene expression and peak accessibility for K562 multi-modal single-cell data with 200 cells and 5 genes
#' @name K562_mini
#' @docType data
#' @description
#' A list containing normalized gene expression (in SCT) and peak accessibility (in ATAC) for K562 multi-modal single-cell data with 200 cells and 5 genes.
#' 
#'
NULL


#' Coaccessibility from cicero for running a demo
#'
#' @format a dataset with 510 rows and 3 columns:
#' \describe{
#'   \item{Peak1}{First peak input for cicero}
#'   \item{Peak2}{Second peak input for cicero}
#'   \item{gene}{gene name}
#'   \item{peak}{cis-peak region}
#'   \item{coaccess}{coaccessibility computed by cicero between the gene promoter and peak}
#'   \item{ind}{indicator to be 1 or 2 indicating which peak before is the promoter of the gene}
#' }
"coaccess.adapt.data.K562"

