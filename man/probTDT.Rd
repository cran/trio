\name{probTDT}
\alias{probTDT}

\title{
TDT on genotype probabilities matrix
}
\description{
Computes the genotypic TDT for a a matrix representing SNP genotype probabilities.
}

\usage{
probTDT(mat.geno, model = c("additive", "dominant", "recessive"), 
   size = 50)
}

\arguments{
  \item{mat.geno}{a numeric matrix with one row for each SNP and \eqn{9 * t} columns representing genotype probabilities for \eqn{t} trios. Each of the \eqn{t}
    blocks (i.e. \code{snp[1:9]}, \code{snp[10:18]}, ...) must consist of sets of the three genotype probabilities for AA, AB and BB calls, of father, mother, and offspring
    (in this order), as would be output by BEAGLE, for example. The genotype probabilities must sum to 1 (up to slight imprecision) in each individual. Missing values are allowed and need to be coded by \code{NA}.
    Note that the order of the columns is not checked to be in terms of minor allele -- any dominant or recessive tests are for allele B, 
	as ordered in the \code{mat.geno}, not necessarily for the minor allele.
}
  \item{model}{type of model that should be fitted. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will
    fit a dominant model, and \code{model = "r"} an recessive model. See description of \code{mat.geno} for a caveat about allele ordering.
}
  \item{size}{the number of SNPs considered simultaneously when computing the parameter estimates. Ignored if \code{fast = FALSE}.
}
  \item{...}{ignored.}
}
\value{
An object of class \code{colTDT} consisting of the following numeric values or vectors, respectively:

\item{coef}{the estimated parameter,}
\item{se}{the estimated standard deviation of the parameter estimate,}
\item{stat}{Wald statistic,}
\item{OR}{the odds ratio, i.e.\ \code{exp(coef)}},
\item{lowerOR}{the lower bound of the 95\% confidence interval for \code{OR},}
\item{upperOR}{the upper bound of the 95\% confidence interval for \code{OR},}
\item{usedTrios}{the number of trios affecting the parameter estimation,}
\item{pMendelErr}{the sum across families of probabilities of Mendelian errors,}
\item{...}{further internal parameters}
}

\references{
Schwender, H., Taub, M.A., Beaty, T.H., Marazita, M.L., and Ruczinski, I. (2011). 
Rapid Testing of SNPs and Gene-Environment Interactions in Case-Parent Trio Data Based on 
Exact Analytic Parameter Estimation. \emph{Biometrics}. DOI: 10.1111/j.1541-0420.2011.01713.x.

Taub M.A., Schwender H., Beatty T.H., Louis T.A., Ruczinski I. (2012).  
Incorporating genotype uncertainties into the genotypic TDT for main effects and gene-environment 
interactions.  \emph{Genetic Epidemiology}, 36(3):225-34. DOI: 10.1002/gepi.21615
}

\author{
Margaret Taub, \email{mtaub@jhsph.edu}
}

\seealso{
\code{\link{tdt}}
}

\examples{
# Load the simulated data.
data(trio.data)

# All SNPs in prob.mat.test can be tested by
prob.tdt.out <- probTDT(prob.mat.test)

# By default, an additive mode of inheritance is considered.
# If another mode, e.g., the dominant mode, should be 
# considered, then this can be done by
prob.tdt.out2 <- probTDT(prob.mat.test, model = "dominant")

\dontrun{
# By default, statistics for the top 5 SNPs are displayed.
# If another number of SNPs, say 10, should be displayed,
# then this can be done by
print(prob.tdt.out2, top = 10)

# The statistics for all SNPs (not ordered by their 
# significance) can be obtained by
print(prob.tdt.out2, top = 0)
}
}

\keyword{htest}
\keyword{models}
\keyword{array}
