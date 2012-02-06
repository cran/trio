\name{tdt}
\alias{tdt}
\alias{colTDT}
\alias{print.tdt}
\alias{print.colTDT}

\title{
Genotypic TDT
}
\description{
Computes the genotypic TDT for a SNP or for each column of a matrix representing a SNP.
}

\usage{
tdt(snp, model = c("additive", "dominant", "recessive"))

colTDT(mat.snp, model = c("additive", "dominant", "recessive"), 
   size = 50)

\method{print}{tdt}(x, digits = 4, ...)

\method{print}{colTDT}(x, top = 5, digits = 4, ...)
}

\arguments{
  \item{snp}{a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks (i.e. \code{snp[1:3]}, \code{snp[4:6]}, ...) must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    The vector must thus have the same structure as the output of \code{\link{trio.check}}, or the genotype example data
    sets such as \code{trio.gen1} (see \code{data(trio.gen1)}), and can be generated from a ped-file by, e.g., employing 
    \code{\link{ped2geno}}.
}
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each of the SNPs must have the same structure
    as \code{snp}, and can, e.g., be generated from a ped-file by employing \code{\link{ped2geno}}.
}
  \item{model}{type of model that should be fitted. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will
    fit a dominant model, and \code{model = "r"} an recessive model. 
}
  \item{size}{the number of SNPs considered simultaneously when computing the parameter estimates. Ignored if \code{fast = FALSE}.
}
  \item{x}{an object of class \code{tdt} or \code{colTDT}, i.e. the output of the function \code{tdt} (or \code{tdtGxG})
    or the function \code{colTDT}.
}
  \item{digits}{number of digits that should be printed.}
  \item{top}{number of interactions that should be printed. If \code{top} is less than or equal to zero, set to \code{NA}, or
    larger than the number of SNPs, then the statistics for all SNPs are printed in the order as they were in the genotype matrix
    used as input into \code{colTDT}. Otherwise, the \code{top} interactions with the smallest p-values are printed.
  }
  \item{...}{ignored.}
}
\value{
An object of class \code{tdt} or \code{colTDT} consisting of the following numeric values or vectors, respectively:

\item{coef}{the estimated parameter,}
\item{se}{the estimated standard deviation of the parameter estimate,}
\item{stat}{Wald statistic,}
\item{OR}{the odds ratio, i.e.\ \code{exp(coef)}},
\item{lowerOR}{the lower bound of the 95\% confidence interval for \code{OR},}
\item{upperOR}{the upper bound of the 95\% confidence interval for \code{OR},}
\item{usedTrios}{the number of trios affecting the parameter estimation (only for \code{colTDT}),}
\item{...}{further internal parameters}
}

\references{
Schwender, H., Taub, M.A., Beaty, T.H., Marazita, M.L., and Ruczinski, I. (2011). 
Rapid Testing of SNPs and Gene-Environment Interactions in Case-Parent Trio Data Based on 
Exact Analytic Parameter Estimation. \emph{Biometrics}. DOI: 10.1111/j.1541-0420.2011.01713.x.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{tdt2way}}, \code{\link{ped2geno}}
}
\keyword{htest}
\keyword{models}
\keyword{array}
