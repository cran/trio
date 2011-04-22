\name{tdt}
\alias{tdt}
\alias{colTDT}
\alias{print.tdt}
\alias{print.colTDT}

\title{
Genotypic TDT
}
\description{
Computes the genotic TDT for a SNP or for each column of a matrix representing a SNP.
}

\usage{
tdt(snp, model = c("additive", "dominant", "recessive"))

colTDT(mat.snp, model = c("additive", "dominant", "recessive"), 
   fast = TRUE, size = 50, warnError = TRUE)

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
  \item{fast}{should closed-form parameter estimates be used? If \code{FALSE}, the function \code{clogit} is used to
    fit the conditional logistic regression model. The latter ``traditional" way of performing a gTDT will be removed
    in one of the next versions of \code{trio}.
}
  \item{size}{the number of SNPs considered simultaneously when computing the parameter estimates. Ignored if \code{fast = FALSE}.
}
  \item{warnError}{logical indicating whether the statistics for the gTDT should be returned as \code{NA} if the fitting of the
    conditional logistic regression model fails. This might in particular happen when the considered SNP is almost monomorphic,
    i.e. if the minor allele frequency of this SNP is very low.
}

  \item{x}{an object of class \code{tdt} or \code{colTDT}, i.e. the output of the function \code{tdt} (or \code{tdt2way})
    or the function \code{colTDT}.
}
  \item{digits}{number of digits that should be printed.}
  \item{top}{number of interactions that should be printed. If the number of interactions is smaller than or equal to 
    \code{top}, then the statistics for all interactions are printed in the order of their computation. Otherwise, 
    the \code{top} interactions with the smallest p-values are printed.
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
\item{...}{further internal parameters}
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
