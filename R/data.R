#' Eisenberg's housekeeping genes
#'
#' @format ## `eisenberg`
#' a vector with house keeping genes
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/23810203/>
"eisenberg"

#' HSAIO's housekeeping genes
#'
#' @format ## `hsiao`
#' a vector with house keeping genes
#'
#' @source <https://www.gsea-msigdb.org/gsea/msigdb/cards/HSIAO_HOUSEKEEPING_GENES>
"hsiao"

#' LIN500's housekeeping genes
#'
#' @format ## `lin500`
#' a vector with house keeping genes
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/23810203/>
"lin500"

#' Eisenberg's mouse housekeeping genes
#'
#' @format ## `eisenberg_mouse`
#' a vector with house keeping genes
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/23810203/>
"eisenberg_mouse"

#' HSAIO's mouse housekeeping genes
#'
#' @format ## `hsiao_mouse`
#' a vector with house keeping genes
#'
#' @source <https://www.gsea-msigdb.org/gsea/msigdb/cards/HSIAO_HOUSEKEEPING_GENES>
"hsiao_mouse"

#' LIN500's mouse housekeeping genes
#'
#' @format ## `lin500_mouse`
#' a vector with house keeping genes
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/23810203/>
"lin500_mouse"

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
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major_Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"example_meta"

#' GSE15907_matrix
#'
#' @format ## `GSE15907_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907>
#'
"GSE15907_matrix"

#' GSE15907_meta
#'
#' @format ## `GSE15907_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE15907_meta"

#' GSE37448_matrix
#'
#' @format ## `GSE37448_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37448>
#'
"GSE37448_matrix"

#' GSE37448_meta
#'
#' @format ## `GSE37448_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE37448_meta"

#' GSE60336_matrix
#'
#' @format ## `GSE60336_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60336>
#'
"GSE60336_matrix"

#' GSE60336_meta
#'
#' @format ## `GSE60336_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE60336_meta"

#' GSE60337_matrix
#'
#' @format ## `GSE60337_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60337>
#'
"GSE60337_matrix"

#' GSE60337_meta
#'
#' @format ## `GSE60337_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE75195_meta"

#' GSE75195_matrix
#'
#' @format ## `GSE75195_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75195>
#'
"GSE75195_matrix"

#' GSE75195_meta
#'
#' @format ## `GSE75195_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE75195_meta"

#' GSE75202_matrix
#'
#' @format ## `GSE75202_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75202>
#'
"GSE75202_matrix"

#' GSE75202_meta
#'
#' @format ## `GSE75202_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE75202_meta"

#' GSE75203_matrix
#'
#' @format ## `GSE75203_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75203>
#'
"GSE75203_matrix"

#' GSE75203_meta
#'
#' @format ## `GSE75203_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE75203_meta"

#' GSE75203_matrix
#'
#' @format ## `GSE75203_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75203>
#'
"GSE75203_matrix"

#' GSE75203_meta
#'
#' @format ## `GSE75203_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE75203_meta"

#' GSE112876_matrix
#'
#' @format ## `GSE112876_matrix`
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112876>
#'
"GSE112876_matrix"

#' GSE112876_meta
#'
#' @format ## `GSE112876_meta`
#' a data frame containing sample level information
#' \describe{
#'   \item{SampleID}{sample ID}
#'   \item{Study}{GSE study that the sample came from}
#'   \item{MajorCellType}{Major cell type that was assigned to the sample}
#'   \item{Major Lineage}{Major lineage that was assigned to the sample}
#'   \item{Cell Family}{Cell Family that was assigned to the sample}
#'   \item{GEOname}{Unique GEO name for each sample}
#'   \item{Organism}{Species sample was taken from}
#'   \item{GeneticBackground}{Organism strain}
#'   \item{Sex}{Sex of organism that the sample was taken from}
#'   \item{Tissue Source}{Tissue source for the sample}
#'   \item{Exclusion Marker}{Any exclusion markers used to differentiate cells}
#'   \item{Treatment}{Treatment administered to the organism}
#'   \item{CellType}{Specific cell type that was assigned to the sample}
#'   \item{Phenotype_Marker}{Phenotypic marker used to identify the cell family}
#'   \item{Separation Method}{Method of separation used for each sample}
#' }
#'
"GSE112876_meta"
