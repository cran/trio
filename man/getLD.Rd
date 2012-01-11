\name{getLD}
\alias{getLD}
\alias{mat}

\title{
 Computation of LD Measures
}
\description{
 Computes the value of D' and r^2 for each pair of SNPs in a matrix.
}
\usage{
  getLD(x, which = c("both", "rSquare", "Dprime"), parentsOnly = FALSE,
     iter = 50, snp.in.col = TRUE, asMatrix = FALSE, addVarN = FALSE)
}

\arguments{
  \item{x}{a numeric matrix consisting of 0, 1, and 2, where it is assumed that
     the values represent the numbers of minor alleles that the SNPs show. Missing values
     are allowed. By default, each column represents a SNP and each row a subject.
     This can be changed by setting \code{snp.in.col = FALSE}. The SNPs should be ordered
     by their position on the considered chromosome.
}
  \item{which}{which LD measures should be computed? Either \code{"rSquare"}, or \code{"Dprime"},
    or the values of \code{"both"} measures are computed. The latter is the default.
}
  \item{parentsOnly}{logical indicating whether only the genotypes of the parents, i.e.\ rows 1, 2,
    4, 5, ... of \code{x}, should be used in the computation of the LD measures when \code{x}
    is in genotype format and contains case-parent trio data (see \code{\link{ped2geno}} and \code{\link{read.pedfile}}).
    If \code{FALSE} (default), all rows are used in the determination of the pairwise LD measure.
} 
  \item{iter}{integer specifying how many iterations are used in the procedure of Hill (1974)
    which is used to estimate D.
}
  \item{snp.in.col}{logical indicating whether each column of \code{x} represents a SNP (and each row a
    subject). If \code{FALSE}, each row represents a SNP (and each column a subject).
}
  \item{asMatrix}{logical indicating whether the LD values are returned as a \eqn{m} x \eqn{m} matrix,
    where \eqn{m} is the number of SNPs. If \code{FALSE}, the LD values are returned as a vector
    of length \eqn{m * (m - 1) / 2}.
}
  \item{addVarN}{logical indicating whether for each pair of SNPs the number of non-missing values
    and the variance estimates of D' proposed by Zabaleta et al. (1997) should be added to the output. The
    variance estimates are required for the identification of LD-blocks with \code{\link{findLDblocks}}.
}
}

\value{
An object of class \code{getLD} consisting (depending of the specification of \code{which}) the
D' (\code{Dprime}) or r^2 (\code{rSquare}) values for each SNP pair, and (depending of the specification
of \code{addVarN}) the variance estimates for D' (\code{varDprime}) and the numbers of non-missing values
(\code{n}). Furthermore, the names of the SNPs (\code{rn}) will be added if \code{asMatrix = FALSE}.

}
\references{
  Hill, W.O. (1974). Estimation of Linkage Disequilibrium in Randomly Mating Populations. 
    \emph{Heredity}, 33, 229-239.
  
  Zapata, C., Alvarez, G., and Carollo, C. (1997). Approximate Variance of the Standardized Measure of Gametic
     Disequilibrium D'. \emph{American Journal of Human Genetics}, 61, 771-774.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{plot.getLD}}, \code{\link{findLDblocks}}
}

\keyword{array}
\keyword{univar}
