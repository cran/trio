\name{trio.permTest}
\alias{trio.permTest}

\title{
Permutation Tests for Trio Logic Regression
}
\description{
Performs either a null-model or a conditional permutation test for a trio logic regression analysis. 
}
\usage{
trio.permTest(object, conditional = FALSE, n.perm = 10, nleaves = NULL, 
   control = NULL, rand = NA)
}

\arguments{
  \item{object}{an object of class \code{trioLR}, i.e. the output of the function \code{\link{trioLR}}. This
     object must be the result of a trio logic regression analysis in which a single model has been fitted (i.e.
     in \code{\link{trioLR}}, \code{search} must have been set to \code{"sa"} and \code{nleaves} must have been
     a single integer).
}
  \item{conditional}{should the conditional permutation test be performed? If \code{FALSE}, a null-model permutation test
    is done analogously to the null-model permutation test for a standard logic regression for population-based data implemented
    in the function \code{logreg} of the \code{R} package \code{LogicReg}. If \code{TRUE}, a test analogous to the conditional 
    permutation test for a standard logic regression is performed.
}
  \item{n.perm}{integer specifying the number of permutations.
}
  \item{nleaves}{integer specifying the maximum number of leaves that the logic tree in the trio logic regression model is allowed to have.
    If \code{NULL}, the maximum number of leaves saved in \code{object} is used.
}
  \item{control}{a list containing the control parameters for the search algorihms and the logic tree considered in \code{\link{trioLR}}, where
    the parameters for an MCMC run and the logic tree are ignored. If \code{NULL} (i.e. by default), the same values for the parameters are
    used that have been employed in the original analysis with \code{\link{trioLR}}. If other values should be used, it is highly recommended
    to specify \code{control} by employing \code{\link{lrControl}}.
}
  \item{rand}{an integer. If specified, the random number generator will be set into a reproducible state.
}
}

\value{
A list consisting of 
\item{origScore}{\code{NA}, if \code{conditional = FALSE}, and otherwise, the score, i.e.\ the value of the partial likelihood, 
   of the original model saved in \code{object}},
\item{permScore}{a vector of length \code{n.perm} containing the scores for the trio logic regression models built in the iterations
  of the permutation test.}
}
\references{
Li, Q., Fallin, M.D., Louis, T.A., Lasseter, V.K., McGrath, J.A., Avramopoulos, D., Wolyniec, P.S., Valle, D., Liang, K.Y., Pulver, A.E., 
and Ruczinski, I. (2010). Detection of SNP-SNP Interactions in Trios of Parents with Schizophrenic Children. Genetic Epidemiology, 34, 396-406. 
}
\author{
Qing Li, \email{mail2qing@yahoo.com}. Modified by Holger Schwender.
}

\seealso{
\code{\link{trioLR}}
}

\keyword{htest}
\keyword{models}
