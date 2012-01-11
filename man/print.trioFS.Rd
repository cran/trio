\name{print.trioFS}
\alias{print.trioFS}
\alias{plot.trioFS}
\title{Printing and plotting of a trioFS object
}
\description{
Prints or plots the most important interactions found in a trioFS analysis. 
}
\usage{
\method{print}{trioFS}(x, topX = 5, show.prop = TRUE, coded = FALSE, digits = 2, ...)

\method{plot}{trioFS}(x, topX = 15, show.prop = FALSE, coded = TRUE, cex = 0.9, 
   pch = 16, col = 1, force.topX = FALSE, include0 = TRUE, add.v0 = TRUE, 
   v0.col = "grey50", main = NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{trioFS}, i.e. the output of \code{\link{trioFS}}.
}
  \item{topX}{integer specifying how many interactions should be shown.
     If \code{topX} is larger than the number of interactions contained in
     \code{x}, all the interactions are shown. Additionally to the \code{topX} most important 
     interactions, any interaction having the same importance as the \code{topX} most
     important one are \code{print}ed or (if \code{force.topX = FALSE}) \code{plot}ted.
}
  \item{show.prop}{should the proportions of models containing the respective interactions be
    added to the output (if \code{print} is used)? If the output of \code{trioFS} should be \code{plot}ted,
    then the proportions of models can be plotted instead of the values of the importance measure
    by setting \code{show.prop = TRUE}.  
}
  \item{coded}{should the coded variable names be displayed? Might be useful
      if the actual variable names are pretty long. The coded variable name of
      the \emph{j}-th variable is \code{Xj}.
}
  \item{digits}{number of digits shown in the \code{print}ed output.
}
  \item{cex}{a numeric value specifying the relative size of the text
     and symbols.}
  \item{pch}{specifies the used symbol. See the help of \code{par} for details.}
  \item{col}{the color of the text and the symbols. See the help of \code{par} for
     how colors can be specified.}
  \item{force.topX}{if \code{TRUE} exactly \code{topX} interactions are
      \code{plot}ted. If \code{FALSE} (default) all interactions up to the \code{topX}th
      most important one and all interactions having the same importance
      as the \code{topX}th most important one are \code{plot}ted.
}
  \item{include0}{should the \eqn{x}-axis include zero regardless whether the importances
      of the shown interactions are much higher than 0?}
  \item{add.v0}{should a vertical line be drawn at \eqn{x = 0}? Ignored if 
     \code{include0 = FALSE} and all importances are larger than zero.}
  \item{v0.col}{the color of the vertical line at \eqn{x = 0}. See the help page of
     \code{par} for how colors can be specified.}
  \item{main}{character string naming the title of the plot. If \code{NULL}, a standard title
     is added to the plot.}
  \item{...}{Ignored.}
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{
   \code{\link{trioFS}}
}

\keyword{print}
\keyword{hplot}
