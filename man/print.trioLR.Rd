\name{print.trioLR}
\alias{print.trioLR}
\title{Printing of trioLR Objects}
\description{
Prints information on the trio logic regression model(s) fitted with \code{\link{trioLR}}. 
}

\usage{
\method{print}{trioLR}(x, asDNF=FALSE, posBeta=FALSE, digits = 3, ...)
}

\arguments{
  \item{x}{an object of class \code{trioLR}, i.e.\ the output of \code{\link{trioLR}}.
}
  \item{asDNF}{should the disjunctive normal form of the logic expression represented by the logic tree be printed?
    If \code{FALSE}, the logic expression is printed as found by the search algorithm in trio logic regression.
    An advantage of the disjunctive normal form representation is that the interactions comprised by the logic
    expression are given by the AND-combinations in the disjunctive normal form. Note that not necessarily the
    minimum disjunctive normal form is printed so that all interactions comprised by the model are shown, even if
    some of the interactions are redundant for the evaluating the logic tree.
}
  \item{posBeta}{should the disjunctive normal form be determined as if the sign of the coefficient in trio logic
    regression model is positive? If \code{FALSE}, the sign is ignored when transforming the logic tree into
    its disjunctive normal form. If \code{TRUE} and the coefficient is negative, the complement of the logic expression is
    transformed into its disjunctive normal form and the coefficient is multiplied by -1. Ignored if \code{asDNF = FALSE}
    or the fitted logic tree only contains one leaf.
}
  \item{digits}{number of digits used in the printing of the score and the parameter estimate of the fitted trio
    logic regression model(s).
} 
  \item{\dots}{ignored.
}
}

\references{
Kooperberg, C. and Ruczinski, I. (2005). Identifying Interacting SNPs Using Monte Carlo Logic Regression.
\emph{Genetic Epidemiology}, 28, 157-170.

Li, Q., Fallin, M.D., Louis, T.A., Lasseter, V.K., McGrath, J.A., Avramopoulos, D., Wolyniec, P.S., Valle, D., 
Liang, K.Y., Pulver, A.E., and Ruczinski, I. (2010). Detection of SNP-SNP Interactions in Trios of Parents 
with Schizophrenic Children. \emph{Genetic Epidemiology}, 34, 396-406.

Ruczinski, I., Kooperberg, C., and LeBlanc, M.L. (2003). Logic Regression. \emph{Journal of Computational and 
Graphical Statistics}, 12, 475-511.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}, based on the \code{plot} functions
implemented by Ingo Ruczinski and Charles Kooperberg in the \code{R} package \code{LogicReg}.
}


\seealso{\code{\link{trioLR}}
}
\keyword{print}
