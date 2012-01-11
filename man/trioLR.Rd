\name{trioLR}
\alias{trioLR}
\alias{trioLR.default}
\alias{trioLR.trioPrepare}

\title{Trio Logic Regression
}
\description{Performs a trio logic regression analysis as proposed by Li et al. (2011), where trio logic regression is an adaptation of
logic regression (Ruczinski et al., 2003) for case-parent trio data.
}
\usage{
\method{trioLR}{default}(x, y, search = c("sa", "greedy", "mcmc"), nleaves = 5, 
   penalty = 0, weights = NULL, control=lrControl(), rand = NA, ...)

\method{trioLR}{trioPrepare}(x, ...)
}

\arguments{
  \item{x}{either an object of class \code{trioPrepare}, i.e. the output of \code{\link{trio.prepare}}, or 
      a binary matrix consisting of zeros and ones. If the latter, then each column of \code{x} must correspond to a binary variable 
     (e.g., codng for a dominant or a recessive effect of a SNP), and each row to a case or a pseudo-control,
     where each trio is represented by a block of four consecutive rows of \code{x} containing the data for the case
     and the three matched pseudo-controls (in this order) so that the first four rows of \code{x} comprise
     the data for the first trio, rows 5-8 the data for the seocnd trio, and so on. Missing values are not
     allowed. A convenient way to generate this matrix is to use the function \code{\link{trio.prepare}}. Afterwards, \code{trioLR} 
     can be directly applied to the output of \code{\link{trio.prepare}}. 
}
  \item{y}{a numeric vector specifying the case-pseudo-control status for the observations in \code{x} (if \code{x} is the binary matrix).
     Since in trio logic regression, cases are coded by a \code{3} and pseudo-controls by a \code{0},
     \code{y} is given by \code{rep(c(3, 0, 0, 0), n.trios)}, where \code{n.trios} is
     the number of trios for which genotype data is stored in \code{x}. Thus, the length of \code{y}
     must be equal to the number of rows in \code{x}. No missing values are allowed in \code{y}.
     If not specified, \code{y} will be automatically generated. 
}
  \item{search}{character string naming the search algorithm that should be used in the search for the best
     trio logic regression model. By default, i.e. \code{search = "sa"}, simulated annealing, the standard
     search algorithm for a logic regression is used. In this case, depending on the length of \code{nleaves}, 
     either one trio logic regression model is fitted or several trio logic regression
     models of different sizes are fitted. For details, see \code{nleaves}.
     Alternatively, a greedy search can be used by setting \code{search = "greedy"}, or a MC logic regression
     analysis (Kooperberg and Ruczinski, 2005) for case-parent trio data can be performed by setting \code{search = "mcmc"}.
}
  \item{nleaves}{integer or vector of two integers specifying the maximum number of leaves, i.e.\ variables,
     in the logic tree of the trio logic regression model (please note in trio logic regression the model
     consists only of one logic tree).
     Must be a single integer, if \code{search = "greedy"} or \code{search = "mcmc"}. If \code{search = "sa"},
     it can also be a vector of two integers, where the second integer must be larger than the first one. In
     this case, several trio logic regression models are fitted in which the maximum numbers of leaves range from
     \code{nleaves[1]} to \code{nleaves[2]}. 
}
  \item{penalty}{a non-negative value for the \code{penalty} parameter used in logic regression. The penalty takes the
     form \code{penalty} times the number of leaves in the model. By defaulty, larger models are not penalized. 
     \code{penalty} is only relevant when one logic regression model is fitted. 
}
  \item{weights}{a numeric vector containing one weight for each trio considered in \code{x}. Thus, \code{weights}
     must contain \code{nrow(x) / 4} positive values. By default, all trios are equally weighted.  
}
  \item{control}{a list of control parameters for the search algorithms and the logic tree considered when fitting a
     (trio) logic regression model. For these parameters, see \code{\link{lrControl}}, which is the function that
     should be used to specify \code{control}.
}
  \item{rand}{integer. If specified, the random number generator will be set into a reproducible state.
}
  \item{...}{for the \code{trioPrepare} method, optional parameters to be passed to 
    the low level function \code{trioLR.default}, i.e. all arguments of \code{trioLR.default} except for
    \code{x} and \code{y}. Otherwise, ignored.
  }
}
\details{Trio logic regression is an adaptation of logic regression to case-parent trio data. Virtually all
  features for a standard logic regression analysis with the function \code{logreg} available in the \code{R}
  package \code{LogicReg} are also available for a trio logic regression analysis,
  either directly via \code{trioLR} or via the function \code{trio.permTest} for performing permutation tests.
  
  For a detailed, comprehensive description on how to perform a logic regression analysis, and thus, a trio
  logic regression analysis, see the \code{Details} section of the help page for the function \code{logreg} 
  in the \code{R} package \code{LogicReg}.  For a detailed explanation on how to specify the parameters for 
  simulated annealing, see the man page of the function \code{logreg.anneal.control} in the \code{R} package
  \code{LogicReg}. 
  
  Finally, an example for a trio logic regression analysis is given in the vignette \code{trio} available
  in the \code{R} package \code{trio}. 
}

\value{
An object of class \code{trioLR} composed of the same objects as an object of class \code{logreg}. For details,
see the \code{Value} section of the function \code{logreg} from the \code{R} package \code{LogicReg}.
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

\author{Holger Schwender, \email{holger.schwender@udo.edu}, based on \code{R} and Fortran code by Ingo Ruczinski, \email{ingo@jhu.edu},
and Charles Kooperberg for a general logic regression as well as Fortran code by Ingo Ruczinski and Qing Li for trio logic regression.
}

\seealso{
\code{\link{trio.prepare}}, \code{\link{trio.check}}, \code{\link{trio.permTest}}
}

\keyword{tree}
\keyword{multivariate}
\keyword{regression}
