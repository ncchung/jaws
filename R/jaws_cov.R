#' The Jackstraw Weighted Shrinkage Estimation Method for High-Dimensional Covariance Matrices
#'
#' Estimates a covariance matrix of \code{m} variables (rows) from \code{n} samples (columns), when \code{m > n}.
#' Shrinkage estimators of principal component loadings are used to construct a high-dimensional covariance matrix.
#' Although several options are available to control characteristics of jackstraw weighted shrinkage (see \code{jaws.pca}),
#' the required inputs are the data matrix \code{dat} and the number of principal components \code{r} to be used.
#'
#' By default, \code{jaws.cov} computes two canonical jackstraw weighted shrinkage estimators, namely \code{PIP} and \code{PNV}.
#' Additionally, extra shrinkage techniques may apply, such as soft- or hard-thresholding posterior inclusion probabilities
#' \code{extra.shrinkage=c("PIPsoft","PIPhard")}.
#' Please provide \code{r} numerical threshold values to be applied to \code{r} principal components.
#'
#' This algorithm applies shrinkage to the signal component of the covariance matrix, and assumes the independently distributed noise.
#' Since this function relies on shrunken loadings of PCs, you may first run \code{jaws.pca} on \code{dat} with a greater control over optional arguments
#' and supply its output \code{jaws.pca.obj} to this function.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param r a number of significance principal components (\code{r < n}).
#' @param jaws.pca.obj a jaws.pca output (optional).
#' @param stat.shrinkage PNV shrinkage may be applied to "F-statistics" or "loadings" (default: F-statistics).
#' @param extra.shrinkage extra shrinkage methods may be used; see details below (optional).
#' @param verbose a logical specifying to print the progress (default: TRUE).
#' @param seed a seed for the random number generator (optional).
#'
#' @return \code{jaws.cov} returns a list consisting of
#' \item{PIP}{estimated covariance matrix based on posterior inclusion probabilities}
#' \item{PNV}{estimated covariance matrix based on proportion of null variables}
#' @return With appropriate \code{extra.shrinkage} options (for details, see the Supplementary Information of Chung and Storey (2013), the output may also include
#' \item{PIPhard}{estimated covariance matrix based on hard-thresholding the \code{PIP} loadings}
#' \item{PIPsoft}{estimated covariance matrix based on soft-thresholding the \code{PIP} loadings}
#'
#' @export jaws.cov
#' @importFrom corpcor fast.svd
#' @import jackstraw
#' @importFrom qvalue pi0est lfdr qvalue
#' @author Neo Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Forthcoming
#'
#' @seealso \link{jaws.pca} \link{jackstraw.PCA}
#'
#' @examples
#' set.seed(1234)
#' ## simulate data from a latent variable model: Y = BX + E
#' B = c(rep(1,50),rep(-1,50), rep(0,900))
#' X = rnorm(20)
#' E = matrix(rnorm(1000*20), nrow=1000)
#' dat = B %*% t(X) + E
#' dat = t(scale(t(dat), center=TRUE, scale=FALSE))
#'
#' ## estimate large-scale covariance matrix
#' jaws.cov.out = jaws.cov(dat, r=1)
jaws.cov <- function(dat, r=NULL, jaws.pca.obj=NULL, stat.shrinkage="F-statistics", extra.shrinkage=NULL, verbose=TRUE, seed=NULL) {
  #jaws.pca.obj=NULL; extra.shrinkage=NULL; verbose=TRUE; seed=NULL
  if(!is.null(seed)) set.seed(seed)
	m=dim(dat)[1]; n=dim(dat)[2]
	if(is.null(r)) {
		warning("The number of significant PCs (r) is missing; this is strongly advised to manually determine r using various statistical and graphical criteria.")
		r = permutationPA(dat=dat, threshold=.05, verbose=TRUE)$r
		message(paste0("Permutation Parallel Analysis, with a liberal threshold of 0.5, estimated r = ", r, "."))
	}
	if(!(r > 0 && r < n)) { stop("A number of significant PCs is not in valid range between 1 and n."); }

	if(is.null(jaws.pca.obj)) {
		jaws.pca.obj = jaws.pca(dat=dat, r=r, stat.shrinkage=stat.shrinkage, extra.shrinkage=extra.shrinkage, verbose=verbose)
	}

	# if(shrink.var.res==TRUE) {
	# 	if(verbose==TRUE) message("Row-wise variances of the residual matrix are estimated by corpcor.")
	# 	if(!is.null(var.res)) stop("To shrink variances of the residual matrix using shrink.var.res=TRUE option, do not use the var.res option.")
	# 	var.res = var.shrink(t(res), verbose=verbose)
	# }

	##############################################################################################################################
 	## Posterior Inclusion Probabilities: Empirical Bayes Shrinkage
	## B.PIP.oracle = (jaws.pca.obj$svd$d[1] * jaws.pca.obj$soft.u) / sqrt(sum(M^2))
	## sqrt(n) is an estimate of the oracle, sqrt(sum(M^2))
	if(verbose==TRUE) message("Computing the Covariance Matrix from the Jaws Shrunken Coefficients")
	e.pc.PIP = jaws.pca.obj$PIP$u[,1:r] %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r) %*% t(jaws.pca.obj$svd$v[,1:r])
	res.PIP = dat - e.pc.PIP
	var.res.PIP = apply(res.PIP, 1, var)
	B.PIP = (jaws.pca.obj$PIP$u %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r))
	cov.PIP = (1/n * B.PIP %*% t(B.PIP)) + diag(var.res.PIP)

	## proportion of null variables: Pi0 Hard Thresholding
	e.pc.PNV = jaws.pca.obj$PNV$u[,1:r] %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r) %*% t(jaws.pca.obj$svd$v[,1:r])
	res.PNV = dat - e.pc.PNV
	var.res.PNV = apply(res.PNV, 1, var)
	B.PNV = (jaws.pca.obj$PNV$u %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r))
	cov.PNV = (1/n * B.PNV %*% t(B.PNV)) + diag(var.res.PNV)

	out = list(call=match.call(), PIP=cov.PIP, PNV=cov.PNV)

	##############################################################################################################################
	## Posterior Inclusion Probabilities + proportion of null variables
	if("PIPhard" %in% extra.shrinkage) {
		e.pc.PIPhard = jaws.pca.obj$PIPhard$u[,1:r] %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r) %*% t(jaws.pca.obj$svd$v[,1:r])
		res.PIPhard = dat - e.pc.PIPhard
		var.res.PIPhard = apply(res.PIPhard, 1, var)
		B.PIPhard = (jaws.pca.obj$PIPhard$u %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r))
		cov.PIPhard = (1/n * B.PIPhard %*% t(B.PIPhard)) + diag(var.res.PIPhard)
		out = c(out, list(PIPhard=cov.PIPhard))
	}

 	## Soft Thresholding PIP based on pi0
	if("PIPsoft" %in% extra.shrinkage) {
		e.pc.PIPsoft = jaws.pca.obj$PIPsoft$u[,1:r] %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r) %*% t(jaws.pca.obj$svd$v[,1:r])
		res.PIPsoft = dat - e.pc.PIPsoft
		var.res.PIPsoft = apply(res.PIPsoft, 1, var)
		B.PIPsoft = (jaws.pca.obj$PIPsoft$u %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r))
		cov.PIPsoft = (1/n * B.PIPsoft %*% t(B.PIPsoft)) + diag(var.res.PIPsoft)
		out = c(out, list(PIPsoft=cov.PIPsoft))
	}

 	## User-defined Soft Thresholding
	if(is.numeric(extra.shrinkage)) {
		e.pc.PIPsoft = jaws.pca.obj$PIPsoft$u[,1:r] %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r) %*% t(jaws.pca.obj$svd$v[,1:r])
		res.PIPsoft = dat - e.pc.PIPsoft
		var.res.PIPsoft = apply(res.PIPsoft, 1, var)
		B.PIPsoft = (jaws.pca.obj$PIPsoft$u %*% diag(jaws.pca.obj$svd$d[1:r], nrow=r))
		cov.PIPsoft = (1/n * B.PIPsoft %*% t(B.PIPsoft)) + diag(var.res.PIPsoft)
		out = c(out, list(PIPsoft=cov.PIPsoft))
	}
    return(out)
}
