#' Estimate sparse loadings (i.e., coefficients) of Principal Component Analysis, Logistic Factor Analysis, and other techniques in the context of Latent Variable Models.
#' Generally, this can facilitate calculation of shrunken R^2 and related quantities that represent estimated latent variables more accurately.
#' Using systematic variation driven by latent variables, this package also estimate covariance matrices of high-dimensional data when a number of rows (variables) is exceedingly larger than a number of observations (columns).
#'
#' Two main functions are \code{\link{jaws.pca}} and \code{\link{jaws.cov}}, which estimate sparse loadings of principal components
#' and large-scale covariance matrix, respectively.
#'
#' \tabular{ll}{
#' Package: \tab jaws\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' License: \tab GPL-2\cr
#' Imports: \tab corpcor, qvalue, jackstraw, lfa\cr
#' }
#'
#'
#' @name jaws-package
#' @aliases jaws
#' @docType package
#' @title Jackstraw Weighted Shrinkage Methods
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}, John D. Storey \email{jstorey@@princeton.edu}
#' @seealso
#' \code{\link{jaws.pca}} \code{\link{jaws.cov}}

NULL
