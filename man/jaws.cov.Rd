% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jaws_cov.R
\name{jaws.cov}
\alias{jaws.cov}
\title{The Jackstraw Weighted Shrinkage Estimation Method for High-Dimensional Covariance Matrices}
\usage{
jaws.cov(dat, r = NULL, jaws.pca.obj = NULL,
  stat.shrinkage = "F-statistics", extra.shrinkage = NULL, verbose = TRUE,
  seed = NULL)
}
\arguments{
\item{dat}{a data matrix with \code{m} rows as variables and \code{n} columns as observations.}

\item{r}{a number of significance principal components (\code{r < n}).}

\item{jaws.pca.obj}{a jaws.pca output (optional).}

\item{stat.shrinkage}{PNV shrinkage may be applied to "F-statistics" or "loadings" (default: F-statistics).}

\item{extra.shrinkage}{extra shrinkage methods may be used; see details below (optional).}

\item{verbose}{a logical specifying to print the progress (default: TRUE).}

\item{seed}{a seed for the random number generator (optional).}
}
\value{
\code{jaws.cov} returns a list consisting of
\item{PIP}{estimated covariance matrix based on posterior inclusion probabilities}
\item{PNV}{estimated covariance matrix based on proportion of null variables}

With appropriate \code{extra.shrinkage} options (for details, see the Supplementary Information of Chung and Storey (2013), the output may also include
\item{PIPhard}{estimated covariance matrix based on hard-thresholding the \code{PIP} loadings}
\item{PIPsoft}{estimated covariance matrix based on soft-thresholding the \code{PIP} loadings}
}
\description{
Estimates a covariance matrix of \code{m} variables (rows) from \code{n} samples (columns), when \code{m > n}.
Shrinkage estimators of principal component loadings are used to construct a high-dimensional covariance matrix.
Although several options are available to control characteristics of jackstraw weighted shrinkage (see \code{jaws.pca}),
the required inputs are the data matrix \code{dat} and the number of principal components \code{r} to be used.
}
\details{
By default, \code{jaws.cov} computes two canonical jackstraw weighted shrinkage estimators, namely \code{PIP} and \code{PNV}.
Additionally, extra shrinkage techniques may apply, such as soft- or hard-thresholding posterior inclusion probabilities
\code{extra.shrinkage=c("PIPsoft","PIPhard")}.
Please provide \code{r} numerical threshold values to be applied to \code{r} principal components.

This algorithm applies shrinkage to the signal component of the covariance matrix, and assumes the independently distributed noise.
Since this function relies on shrunken loadings of PCs, you may first run \code{jaws.pca} on \code{dat} with a greater control over optional arguments
and supply its output \code{jaws.pca.obj} to this function.
}
\examples{
set.seed(1234)
## simulate data from a latent variable model: Y = BX + E
B = c(rep(1,50),rep(-1,50), rep(0,900))
X = rnorm(20)
E = matrix(rnorm(1000*20), nrow=1000)
dat = B \%*\% t(X) + E
dat = t(scale(t(dat), center=TRUE, scale=FALSE))

## estimate large-scale covariance matrix
jaws.cov.out = jaws.cov(dat, r=1)
}
\author{
Neo Chung \email{nchchung@gmail.com}
}
\references{
Chung and Storey (2015) Forthcoming
}
\seealso{
\link{jaws.pca} \link{jackstraw.PCA}
}

