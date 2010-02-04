\name{trio.gen1}
\docType{data}
\alias{trio.gen1}
\alias{trio.gen2}
\alias{trio.gen.err}
\alias{trio.ped1}
\alias{trio.ped2}
\alias{trio.ped.err}
\title{Genotype Data for Trios}
\description{
Different example data set containing genotype data for 10 SNPs in 100 trios in different formats. 

\code{trio.gen1}, \code{trio.gen2}, and \code{trio.gen.err} consist of 12 columns and 300 rows, where 
the first two columns contain family identifier and individual identifier. In the columns afterwards, 
each SNPs is encoded in one variable denoting the number of minor alleles.

\code{trio.ped1}, \code{trio.ped2}, and \code{trio.ped.err} consist of 26 columns and 300 rows, where
the first six columns identify the family structure of the data, and the phenotype. Besides the variables
providing information on the family structure and the phenotypes
(columns 1 to 6), each SNPs is encoded in two variables denoting the alleles.

Contrary to the other data sets, \code{trio.gen.err} and \code{trio.ped.err} contain Mendelian errors. 
}

\author{Qing Li, mail2qing@yahoo.com}

\seealso{\code{\link{trio.check}}, \code{\link{trio}}}

\keyword{datasets}

