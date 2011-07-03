\name{colGxE}
\alias{colGxE}
\alias{print.colGxE}
\title{
Genotypic TDT for Gene-Environment Interactions
}
\description{
Performs a genotypic TDT for gene-environment interactions for each SNP represented by a column of a matrix in genotype format
and a binary environmental factor.
}
\usage{
colGxE(mat.snp, env, model = c("additive", "dominant", "recessive"), 
   size = 50)

\method{print}{colGxE}(x, top = 5, digits = 4, onlyGxE = FALSE, ...)
}

\arguments{
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each column must be
    a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix might be generated from a ped-file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{env}{a vector of length \eqn{t} (see \code{mat.snp}) containing for each offspring the value of a binary environmental
    variable, which must take the values 0 and 1.
}
  \item{model}{type of model that should be fitted. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will
    fit a dominant model, and \code{model = "r"} an recessive model. 
}
  \item{size}{the number of SNPs considered simultaneously when computing the parameter estimates.
}
  \item{x}{an object of class \code{colGxE}, i.e. the output of the function \code{colGxE}.
}
  \item{top}{number of interactions that should be printed. If the number of interactions is smaller than or equal to 
    \code{top}, then the statistics for all interactions are printed in the order of their computation. Otherwise, 
    the \code{top} interactions with the smallest p-values are printed.
}
  \item{onlyGxE}{logical indicating whether only the statistics for the parameter of the GxE interaction should be printed.
    If \code{FALSE}, the statistics for both parameters in the model are shown.
}
  \item{digits}{number of digits that should be printed.}
  \item{...}{ignored.}
}
\details{
A conditional logistic regression model including two parameters, one for \eqn{G}, and the other for \eqn{GxE}, is fitted, where
\eqn{G} is specified according to \code{model}.
}
\value{
An object of class \code{colGxE} consisting of the following numeric matrices with two columns (one for each parameter):

\item{coef}{the estimated parameter,}
\item{se}{the estimated standard deviation of the parameter estimate,}
\item{stat}{Wald statistic,}
\item{OR}{the odds ratio, i.e.\ \code{exp(coef)}},
\item{lowerOR}{the lower bound of the 95\% confidence interval for \code{OR},}
\item{upperOR}{the upper bound of the 95\% confidence interval for \code{OR},}
\item{usedTrios}{the number of trios affecting the parameter estimation,}
\item{type}{\code{model}.}
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{tdt}}, \code{\link{ped2geno}}
}
\keyword{htest}
\keyword{models}
\keyword{array}
