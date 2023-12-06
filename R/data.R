#' Eisenberg's housekeeping genes
#'
#' @format ## `eisenbert`
#' a vector with house keeping genes
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/23810203/>
"eisenberg"

#' HSAIO's housekeeping genes
#'
#' @format ## `hsaio`
#' a vector with house keeping genes
#'
#' @source <https://www.gsea-msigdb.org/gsea/msigdb/cards/HSIAO_HOUSEKEEPING_GENES>
"hsaio"

#' LIN500's housekeeping genes
#'
#' @format ## `lin500`
#' a vector with house keeping genes
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/23810203/>
"lin500"

#' Example expression data
#'
#' @format ## `example_mat`
#' a numeric matrix containing transformed values with samples in columns and gene names in the rows
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112876>
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907>
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37448>
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60336>
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60337>
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75202>
"example_mat"

#' Example meta data
#'
#' @format ## `example_meta`
#' a data frame containing sample level information line name and study
#' \describe{
#'   \item{Name}{sample name}
#'   \itme{Study}{GSE study that the sampel came from}
#'   \item{CellType}{cell type that was assigned to the sample}
#' }
#'
"example_meta"

