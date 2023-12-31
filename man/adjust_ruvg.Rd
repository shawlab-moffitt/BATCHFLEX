% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust_ruvg.R
\name{adjust_ruvg}
\alias{adjust_ruvg}
\title{Adjust RUVg}
\usage{
adjust_ruvg(mat, housekeeping, k, drop, center, round, tolerannce)
}
\arguments{
\item{mat}{Numeric matrix with features as rownames and sample names as the column names}

\item{housekeeping}{Name of housekeeping gene set or character vector of housekeeping genes}

\item{k}{Used in the RUVg method, the number of factors of unwanted variation to be estimated from the data}

\item{drop}{Used in the RUVg method, the number of singular values to drop in the estimation of the factors of unwanted variation. This number is usually zero, but might be set to one if the first singular value captures the effect of interest. It must be less than k}

\item{center}{Used in the RUVg method, if TRUE, the counts are centered, for each gene, to have mean zero across samples. This is important to ensure that the first singular value does not capture the average gene expression}

\item{round}{Used in the RUVg method, if TRUE, the normalized measures are rounded to form pseudo-counts}

\item{tolerance}{Used in the RUVg method, tolerance in the selection of the number of positive singular values, i.e., a singular value must be larger than tolerance to be considered positive}
}
\value{
an adjusted expression matrix
}
\description{
Wrapper for \code{RUVSeq::RUVg}
}
\examples{
set.seed(333)
}
