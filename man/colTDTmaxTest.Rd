\name{colTDTmaxTest}
\alias{colTDTmaxTest}
\alias{colTDTmaxStat}
\alias{print.maxTestTrio}
\alias{print.maxStatTrio}

\title{
Maximum Genotypic TDT
}
\description{
Computes the maximum over the gTDT statistics for an additive, dominant, and recessive model.
\code{colTDTmaxTest} additionally computes permutation-based p-values.
}
\usage{
colTDTmaxTest(geno, perm = 1e+05, size = 50, minimum = 1e-06, verbose = FALSE)
colTDTmaxStat(geno, size = 50)

\method{print}{maxTestTrio}(x, top = 5, digits = 4, ...)
\method{print}{maxStatTrio}(x, top = 5, digits = 4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{a numeric matrix in which each column represents a SNP. Each column must be
    a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix might be generated from a ped-file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{perm}{number of permutations of the response for which the permuted values of the test statistic should be computed.
}
  \item{size}{number of SNPs that should be considered simultaneously when estimating the parameter.
}
  \item{minimum}{minimum value that a test statistic must show that for the corresponding SNP the p-value is computed.
}
  \item{verbose}{logical indicating whether some information on what is currently computed should be printed.
}
  \item{x}{an object of class \code{maxTestTrio} or \code{maxTestStat}, i.e. the output of \code{colTDTmaxTest} or of \code{colTDTmaxStat}.}
  \item{digits}{number of digits that should be printed.}
  \item{top}{number of interactions that should be printed. If the number of interactions is smaller than or equal to 
    \code{top}, then the statistics for all interactions are printed in the order of their computation. Otherwise, they 
    the \code{top} Top interactions are printed.
  }
  \item{...}{ignored.}
}
\value{
For \code{colTDTmaxStat}, an object of class \code{maxStatTrio} consisting of a vector \code{stat} containing the values of the Max
statistic for the SNPs in \code{geno}, a matrix \code{max.stat} containing the values of the gTDT statistic for testing an additive, a dominant,
and a recessive effect, and additional information required by \code{colTDTmaxTest}. 

For \code{colTDTmaxTest}, an object of class \code{maxTestTrio} consisting of \code{stat}, \code{max.stat}, and the unadjusted
p-values \code{pval} corresponding to \code{stat}. 
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{tdt}}
}

\keyword{htest}
\keyword{models}
\keyword{array}