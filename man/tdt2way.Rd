\name{tdt2way}
\alias{tdt2way}
\alias{colTDT2way}
\alias{print.tdtEpi}
\alias{print.colTDTepi}
\alias{colTDTinter2way}

\title{
Genotypic TDT for Two-Way Interactions
}
\description{
Computes the genotypic TDT for the interactions of two SNPs or of each pair of columns of a matrix representing SNPs.
Alternatively, \code{colTDTinter2way} can be used to compute all possible interactions between each SNP from one matrix
\code{mat.snp1} with another SNP \code{snp} or each SNP from another matrix \code{mat.snp2}, 
where the two matrices contain SNPs, e.g., from two different genes.

This test makes use of the 16 possible genotypes that can be obtained from combining the parents' genotypes of the
two considered SNPs. Thus, for each family, genotypes for one case (i.e. the offspring) and 15 pseudo-controls are used.
}

\usage{
tdt2way(snp1, snp2, epistatic = TRUE, model = c("additive", "dominant", "recessive"),
	 add = FALSE, warnError = TRUE)
	 
colTDT2way(mat.snp, epistatic = TRUE, genes = NULL, maf = FALSE,
   model = c("additive", "dominant", "recessive"), add = FALSE, 
   warnError = TRUE)
   
colTDTinter2way(mat.snp1, mat.snp2, snp=NULL, epistatic = TRUE, 
   maf = FALSE, model = c("additive", "dominant", "recessive"), 
   add = FALSE, warnError = TRUE)
}

\arguments{
  \item{snp1, snp2}{numeric vectors of length \eqn{3 * t} representing two SNPs genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks (i.e. \code{snp1[1:3]}, \code{snp1[4:6]}, ..., and \code{snp2[1:3]}, \code{snp2[4:6]}, ...) must consist of the 
    genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    The vectors must thus have the same structure as the output of \code{\link{trio.check}}, or the genotype example data
    sets such as \code{trio.gen1} (see \code{data(trio.gen1)}), and can be generated from a ped-file by, e.g., employing 
    \code{\link{ped2geno}}.
}
  \item{mat.snp, mat.snp1}{a numeric matrix in which each column represents a SNP. Each of the SNPs must have
    the same structure as \code{snp}, and can, e.g., be generated from a ped-file by employing \code{\link{ped2geno}}.
}
  \item{mat.snp2}{a numeric matrix similar to \code{mat.snp1} (or \code{mat.snp}). The interaction between each SNP in this
    matrix and in each SNP in \code{mat.snp1} will be tested. If the interaction between each SNP in \code{mat.snp1} and a 
    single SNP should be tested, \code{snp} needs to be specified instead of \code{mat.snp2}.
}
  \item{snp}{a numeric vector representing a SNP. This vector must have the same format as each columns in \code{mat.snp1}.
    If the interaction between each SNP in \code{mat.snp1} and a single variable should be tested, \code{snp} needs to
    be specified instead of \code{mat.snp2}.
} 
  \item{epistatic}{logical value specifying whether a conditional logistic regression version of the test proposed by Cordell (2002)
    is used to test for epistatistical interactions. If \code{FALSE}, the respective two vectors containing the case and the 15 pseudo-controls for
    each trio are simply multiplied (considering \code{model}).
}
  \item{genes}{a character vector containing the names of the genes to which the SNPs belong. If specified, only the two-way interactions
    between SNPs from different genes are tested. If \code{NULL}, all two-way interactions between all possible pairs of SNPs are tested.
}
  \item{maf}{logical indicating whether the minor allele frequency (computed by considering the genotypes of only the parents) should
    be added to the output.
}
  \item{model}{type of model that should be considered. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will
    consider a dominant model for each of the respective two SNPs, and \code{model = "r"} an recessive model. Ignored if 
    \code{epistatic = TRUE}. 
}
  \item{add}{should the main effects be added to the model, if \code{epistatic = FALSE}? If \code{FALSE}, the conditional
    logistic regression model contains only the interaction as predictor. Otherwise, the two SNPs themselves are also included 
    in the model. Ignored if \code{epistatic = TRUE}.
}
  \item{warnError}{logical indicating whether the statistics for the gTDT should be returned as \code{NA} if the fitting of the
    conditional logistic regression model fails. This might in particular happen when the two considered SNPs are in (strong) LD.
}
}
\value{
Depending on \code{epistatic}, the output contains statistics and p-values either of a likelihood ratio test (\code{epistatic = TRUE})
or the Wald statistic for the conditional logistic regression model (\code{epistatic = FALSE}). If \code{maf = TRUE}, a vector \code{maf}
containing the minor allele frequencies of each SNP and a matrix \code{mat.maf} with two columns containing the SNP-wise minor allele
frequencies for each tested pair of SNPs are added to the output.
}

\references{
Cordell, H. J. (2002). Epistasis: What it Means, what it Doesn't mean, and Statistical Methods to Detect it in Humans.
Human Molecular Genetics, 11, 2463-2468.
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

