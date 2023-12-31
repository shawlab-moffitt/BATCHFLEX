% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust_limma.R
\name{adjust_limma}
\alias{adjust_limma}
\title{Adjust Limma}
\usage{
adjust_limma(
  mat,
  meta,
  treatment,
  batch.1,
  batch.2,
  log2_transformed = TRUE,
  ...
)
}
\arguments{
\item{mat}{Numeric matrix with features as rownames and sample names as the column names}

\item{meta}{Data frame of sample data with the first column being sample names that match the column names of the matrix}

\item{treatment}{Column name from the meta file of the column that will be used for treatment information}

\item{batch.1}{Column name from the meta file of the column that will be used for batch one information}

\item{batch.2}{Column name from the meta file of the column that will be used for batch two information}

\item{log2_transformed}{whether data is log2 transformed already or not}

\item{...}{other parameters to pass to \code{limma::removeBatchEffects}}
}
\value{
an adjusted expression matrix
}
\description{
Wrapper for \code{limma::removeBatchEffects}
}
\examples{
set.seed(333)
}
