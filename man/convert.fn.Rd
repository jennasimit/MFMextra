% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sims.R
\name{convert.fn}
\alias{convert.fn}
\title{Convert genotype calls, as output from hapgen2, to a genotype score matrix}
\usage{
convert.fn(gcalls)
}
\arguments{
\item{gcalls}{data.frame of genotype calls where rows are SNPs and columns are individuals}
}
\value{
genotype score matrix where rows are SNPs and columns are individuals
}
\description{
Convert genotype calls, as output from hapgen2, to a genotype score matrix
}
\author{
Jenn Asimit
}
