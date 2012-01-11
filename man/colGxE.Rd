\name{colGxE}
\alias{colGxE}
\alias{colGxEPerms}
\title{
Genotypic TDT for Gene-Environment Interactions
}
\description{
Performs a genotypic TDT for gene-environment interactions for each SNP represented by a column of a matrix in genotype format
and a binary environmental factor. While \code{colGxE} computes the p-values based on asymptotic ChiSquare-distributions, 
\code{colGxEPerms} can be used to determine permutation-based p-values.
}
\usage{
colGxE(mat.snp, env, model = c("additive", "dominant", "recessive"), 
   size = 50, addGandE = TRUE, whichLRT = c("both", "2df", "1df", "none"),
   add2df = TRUE, addCov = FALSE, famid=NULL)
   
colGxEPerms(mat.snp, env, model = c("additive", "dominant", "recessive"),
   B = 10000, size = 20, addPerms = TRUE, famid = NULL, rand = NA)
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
  \item{addGandE}{should the ORs and their confidence intervals for the exposed cases be added to the output?}
  \item{whichLRT}{character string specifying which likelihood ratio test should be added to the output. If \code{"2df"}, 2 degree of
    freedom likelihood ratio tests comparing the fitted models (containing one parameter for the SNP and one for the gene-environment
    interaction) with models containing no factor will be performed. If \code{"1df"}, one degree of freedom likelihood ratio
    tests comparing the fitted model (containing two parameters, one for the SNP and the other for the interaction) with models
    only containing the respective SNP will be added to the output. If \code{"both"} (default), both tests will be performed,
    whereas none test will be done, if \code{whichLRT = "none"}.}
  \item{add2df}{should the results of a 2 df Wald test for testing both the SNP and the interaction effect simultaneously be
    added to the model?}
  \item{addCov}{should the covariance between the parameter estimations for the SNP and the gene-environment interaction be added
    to the output? Default is \code{addCov = FALSE}, as this covariance is given by the negative variance of the parameter estimate
    for the SNP.}
  \item{famid}{a vector of the same length as \code{env} specifying the family IDs for the corresponding values of the environmental
    variable in \code{env}. Can be used to reorder the vector \code{env} when the order of the trios differs between \code{env} and
    \code{mat.snp}.
}
  \item{B}{number of permutations.}
  \item{addPerms}{should the matrices containing the permuted values of the test statistics for the SNP and the gene-environment
    interaction be added to the output?}
  \item{rand}{integer for setting the random number generator into a reproducible state.}
}
\details{
A conditional logistic regression model including two parameters, one for \eqn{G}, and the other for \eqn{GxE}, is fitted, where
\eqn{G} is specified according to \code{model}.
}
\value{
For \code{colGxE}, an object of class \code{colGxE} consisting of the following numeric matrices with two columns (one for each parameter):

\item{coef}{the estimated parameter,}
\item{se}{the estimated standard deviation of the parameter estimate,}
\item{stat}{Wald statistic,}
\item{OR}{the odds ratio, i.e.\ \code{exp(coef)}},
\item{lowerOR}{the lower bound of the 95\% confidence interval for \code{OR},}
\item{upperOR}{the upper bound of the 95\% confidence interval for \code{OR},}
\item{usedTrios}{the number of trios affecting the parameter estimation,}
\item{env}{vector containing the values of the environmental factor,}
\item{type}{\code{model},}
\item{addGandE}{the value of \code{addGandE},}
\item{addOther}{a logical vector specifying which of the likelihood ratio tests and if the 2 df Wald test was performed,}

and depending on the specifications in \code{colGxE}

\item{cov}{numeric vector containing the covariances,}
\item{lrt2df}{a numeric matrix with two columns, in which the first column contains the values of the 1 df likelihood ratio test statistic
  and the second the corresponding p-values,}
\item{wald2df}{a numeric matrix with two columns, in which the first column contains the values of the 2 df Wald test statistics
  and the second the corresponding p-values,}
\item{lrt1df}{a numeric matrix with two columns, in which the first column contains the values of the 2 df likelihood ratio test statistic
  and the seocnd the corresponding p-values.}

For \code{colGxEPerms}, 
\item{stat}{a matrix with two columns containing the values of gTDT statistics for the main effects of the SNPs and the
  gene-environment interactions when considering the original, unpermuted case-pseudo-control status,}
\item{pval}{a matrix with two columns comprising the permutation-based p-values corresponding to the test statistics in \code{stat},}

and if \code{addPerms = TRUE}

\item{matPermG}{a matrix with \code{B} columns containing the values of the gTDT statistic for the SNPs when considering
  the \code{B} permutations of the case-pseudo-control status,}
\item{matPermGxE}{a matrix with \code{B} columns containing the values of the gTDT statistic for the gene-environment interactions
  when considering the \code{B} permutations of the case-pseudo-control status.}
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{colTDT}}, \code{\link{ped2geno}}
}
\keyword{htest}
\keyword{models}
\keyword{array}
