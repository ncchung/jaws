#' The Jackstraw Weighted Shrinkage Estimation Method for Sparse Loadings in Principal Component Analysis
#'
#' Estimates sparse/shrunken loadings of principal component analysis.
#' Based on statistical sginificance of association between variables and principal components,
#' the sample loadings of principal components are shruken towards zeros, which improve its accuracy.
#' The only required inputs are the data matrix \code{dat} and the number of principal components \code{r} whose loadings you would like to estimate.
#'
#' By default, \code{jaws.pca} computes two canonical jackstraw weighted shrinkage estimators, namely \code{PIP} and \code{PNV}.
#' Additionally, other extra shrinkage techniques may apply, such as combining two canonical estimaotrs by setting \code{extra.shrinkage="PIPhard"}
#' and applying soft-thresholding to local fdr by setting \code{extra.shrinkage} to numerical threshold values between 0 and 1.
#' Please provide \code{r} numerical threshold values to be applied to \code{r} principal components.
#'
#' It is strongly advised that you take a careful look at your data and use appropriate graphical and statistical criteria
#' to determine a number of significant PCs, \code{r}. For example, see a contributed R package called `nFactors'.
#' In a case when you fail to specify \code{r}, \code{r} will be estimated from permutation Parallel Analysis (Buja and Eyuboglu, 1992)
#' via a function \link{permutationPA}, with a very liberal threshold.
#'
#' If \code{s} is not supplied, \code{s} is set to about 10\% of \code{m} variables.
#' If \code{B} is not supplied, \code{B} is set to \code{m*10/s}.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param p a \code{m * r} matrix of p-values for association tests between variables and \code{r} principal components, generally computed from the jackstraw method. If \code{p} is not given, \code{jackstraw.PCA} is automatically applied.
#' @param r a number (a positive integer) of significance principal components.
#' @param s a number (a positive integer) of ``synthetic'' null variables (optional).
#' @param B a number (a positive integer) of resampling iterations (optional).
#' @param stat.shrinkage PNV shrinkage may be applied to "F-statistics" or "loadings" (default: F-statistics).
#' @param extra.shrinkage extra shrinkage methods may be used; see details below (optional).
#' @param verbose a logical specifying to print the progress (default: TRUE).
#' @param save.all a logical specifying to save all objects, including a large SVD object (default: FALSE).
#' @param seed a seed for the random number generator (optional).
#'
#' @return \code{jaws.pca} returns a list consisting of
#' \item{p}{p-values for association tests between variables and each of \code{r} principal components}
#' \item{pi0}{proportion of variables not associated with \code{r} principal components, individually}
#' \item{svd}{SVD object from decomposing \code{dat}}
#' \item{PIP}{a list of outputs derived from the posterior inclusion probabilities method (including \code{pr}, \code{u}, \code{var}, \code{PVE})}
#' \item{PNV}{a list of outputs derived from the proportion of null variables method (including \code{pi0}, \code{u}, \code{var}, \code{PVE})}
#' @return With appropriate \code{extra.shrinkage} options (for details, see the Supplementary Information of Chung and Storey (2013), the output may also include
#' \item{PIPhard}{a list of outputs from hard-threshoding the \code{PIP} loadings (including \code{u}, \code{var}, \code{PVE})}
#' \item{PIPsoft}{a list of outputs from soft-threshoding the \code{PIP} loadings (including \code{pr}, \code{u}, \code{var}, \code{PVE})}
#'
#' @section Detailed explanation of the output contained in \code{PIP} and \code{PNV}:
#' \describe{
#'   \item{pr}{a matrix of posterior inclusion probabilities (equivalent to 1-lfdr) for \code{m} coefficients and \code{r} PCs.}
#'   \item{pi0}{a vector of estimated proportion of null variables for \code{r} PCs.}
#'   \item{u}{a \code{m*r} matrix of shrunken loadings.}
#'   \item{var}{a vector of shrunken variances explained by \code{r} PCs.}
#'   \item{PVE}{a vector of shrunken percent variances explained by \code{r} PCs.}
#' }
#'
#' @export jaws.pca
#' @importFrom corpcor fast.svd
#' @import jackstraw
#' @importFrom qvalue pi0est lfdr qvalue
#' @author Neo Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Forthcoming
#'
#' @seealso \link{jackstraw.PCA} \link{jaws.cov}
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
#' ## estimate sparse loadings in PCA
#' jaws.pca.out = jaws.pca(dat, r=1)
#'
jaws.pca <- function(dat, p=NULL, r=NULL, s=NULL, B=NULL, stat.shrinkage="F-statistics", extra.shrinkage=NULL, verbose=TRUE, seed=NULL, save.all=TRUE) {
	if(!is.null(seed)) set.seed(seed)
	m=dim(dat)[1]; n=dim(dat)[2]

	if(is.numeric(extra.shrinkage)) {
		if(length(extra.shrinkage) != r) stop(paste0("For r=",r," principal components, provide r threshold values."))
		if(min(extra.shrinkage) < 0 | max(extra.shrinkage) > 1) { stop(paste0("For soft-thresholding for local FDRs , threshold values must be in the valid range between 0 and 1.")) }
	}

	if(is.null(r)) {
		warning("The number of significant PCs (r) is missing; this is strongly advised to manually determine r using various statistical and graphical criteria.")
		r = permutationPA(dat=dat, threshold=.05, verbose=verbose)$r
		message(paste0("Permutation Parallel Analysis, with a liberal threshold of 0.5, estimated r = ", r, "."))
	}
	if(!(r > 0 && r < n)) { stop("A number of significant PCs is not in valid range between 1 and n."); }

	if(!is.null(p)) p = as.matrix(p)
	if(is.null(p) & verbose==TRUE) message("\nThe association significance between variables and principal components are computed by the jackstraw.")

	if(stat.shrinkage == "F-statistics") {
	  # compute p-values when stat.shrinkage=="F-statistics"
	  p = matrix(NA, nrow=m, ncol=r)
	  Fstat = matrix(NA, nrow=m, ncol=r)
	  for(i in 1:r) {
	    jackstraw.output = jackstraw.PCA(dat, r1=i, r=r, s=s, B=B, verbose=verbose)
	    p[,i] = jackstraw.output$p.value
	    Fstat[,i] = jackstraw.output$obs.stat
	    rm(jackstraw.output)
	  }
	} else if(is.null(p) & stat.shrinkage == "loadings") {
	  # compute p-values when stat.shrinkage=="loadings"
	  p = matrix(NA, nrow=m, ncol=r)
	  for(i in 1:r) {
	    p[,i] = jackstraw.PCA(dat, r1=i, r=r, s=s, B=B, verbose=verbose)$p.value
	  }
	}

	if(ncol(p) != r | nrow(p) != m) stop("The p-value matrix must match m (nrow) and r (ncol).")

	# compute pi0
	# averaging two methods to calculate pi0 seems to be stable and accurate (esp. at a very high pi0)
	pi0 = matrix(0, nrow=2, ncol=r)
	for(i in 1:r) {
		pi0[1,i] = tryCatch(pi0est(p[,i], pi0.method="smoother")$pi0, error=function(tmp) NA)
		pi0[2,i] = tryCatch(pi0est(p[,i], pi0.method="bootstrap")$pi0, error=function(tmp) NA)
	}
	pi0 = apply(pi0, 2, function(x) mean(x, na.rm=TRUE))

	svd.raw = fast.svd(dat)
	var.total = sum(svd.raw$d^2) # =sum(diag(dat %*% t(dat))) may not be computed when m is large
	r2 = r2.est(dat, svd.raw$u[,1:r,drop=FALSE] %*% diag(svd.raw$d[1:r], r) %*% t(svd.raw$v[,1:r,drop=FALSE]))

	##############################################################################################################################
 	## Posterior Inclusion Probabilities: Empirical Bayes Shrinkage
	if(verbose==TRUE) message(paste0("Computing the PIP Shrunken Loadings for the PC (r=", r, ") : "))
	PIP = list(pr=matrix(NA, nrow=m, ncol=r), u=matrix(NA, nrow=m, ncol=r), var=vector("numeric", r), PVE=vector("numeric", r))
	for(i in 1:r) {
		PIP$pr[,i] = 1-lfdr(p[,i], pi0=pi0[i])
		PIP$u[,i] = PIP$pr[,i] * svd.raw$u[,i]
		PIP$var[i] = svd.raw$d[i]^2 * sum(PIP$u[,i]^2)	#=sum(diag(e.PIP %*% t(e.PIP)))
		PIP$PVE[i] = PIP$var[i] / var.total
		if(verbose==TRUE) cat(paste(i," "))
	}
	PIP$r2 = r2.est(dat, PIP$u %*% diag(svd.raw$d[1:r], r) %*% t(svd.raw$v[,1:r,drop=FALSE]))

 	## proportion of null variables: Pi0 Hard Thresholding
	if(verbose==TRUE) message(paste0("Computing the PNV Shrunken Loadings for the PC (r=", r, ") : "))
	PNV = list(pi0=pi0, u=as.matrix(svd.raw$u[,1:r]), var=vector("numeric", r), PVE=vector("numeric", r))
	for(i in 1:r) {
	  if(stat.shrinkage == "F-statistics") {
	    PNV$u[which(rank(abs(Fstat[,i])) <= m*pi0[i]),i] = 0
	  } else if(stat.shrinkage == "loadings") {
			PNV$u[which(rank(abs(PNV$u[,i])) <= m*pi0[i]),i] = 0
	  }

		PNV$var[i] = svd.raw$d[i]^2 * sum(PNV$u[,i]^2) #=sum(diag(e.PNV %*% t(e.PNV)))
		PNV$PVE[i] = PNV$var[i] / var.total
		if(verbose==TRUE) cat(paste(i," "))
	}
	PNV$r2 = r2.est(dat, PNV$u %*% diag(svd.raw$d[1:r], r) %*% t(svd.raw$v[,1:r,drop=FALSE]))

	if(save.all == FALSE) {
		out = list(call=match.call(), p=p, pi0=pi0, PIP=PIP, PNV=PNV, r2=r2)
	} else {
		out = list(call=match.call(), p=p, pi0=pi0, svd=svd.raw, PIP=PIP, PNV=PNV, r2=r2)
	}

	##############################################################################################################################
	## Posterior Inclusion Probabilities + proportion of null variables (Pi0)
	if(verbose==TRUE) message(paste0("Computing the PIPhard Shrunken Loadings for the PC (r=", r, ") : "))
	if("PIPhard" %in% extra.shrinkage) {
		PIPhard = list(u=matrix(NA, nrow=m, ncol=r), var=vector("numeric", r), PVE=vector("numeric", r))
		for(i in 1:r) {
			PIPhard$u[,i] = PIP$pr[,i] * PNV$u[,i]
			PIPhard$var[i] = svd.raw$d[i]^2 * sum(PIPhard$u[,i]^2)
			PIPhard$PVE[i] = PIPhard$var[i] / var.total
			if(verbose==TRUE) cat(paste(i," "))
		}
		out = c(out, list(PIPhard=PIPhard))
	}

 	## Posterior Inclusion Probabilities with Additional Soft Thresholding
	if(verbose==TRUE) message(paste0("Computing the PIPsoft Shrunken Loadings for the PC (r=", r, ") : "))
	if("PIPsoft" %in% extra.shrinkage) {
 		PIPsoft = list(threshold=vector(mode="numeric", length=r), pr=matrix(NA, nrow=m, ncol=r), u=matrix(NA, nrow=m, ncol=r), var=vector("numeric", r), PVE=vector("numeric", r))
		for(i in 1:r) {
			PIPsoft$threshold[i] = max(PIP$pr[which(rank(PIP$pr[,i]) <= round(pi0[i]*m)),i])

			if(PIPsoft$threshold[i]<1) PIPsoft$pr[,i] = sapply(PIP$pr[,i], function(x) max((x - PIPsoft$threshold[i]), 0) / (1-PIPsoft$threshold[i]))
			if(PIPsoft$threshold[i]==1) PIPsoft$pr[,i] = rep(0, m)

			PIPsoft$u[,i] = PIPsoft$pr[,i] * svd.raw$u[,i]
			PIPsoft$var[i] = svd.raw$d[i]^2 * sum(PIPsoft$u[,i]^2)
			PIPsoft$PVE[i] = PIPsoft$var[i] / var.total
			if(verbose==TRUE) cat(paste(i," "))
		}
		out = c(out, list(PIPsoft=PIPsoft))
	}

	if(is.numeric(extra.shrinkage)) {
 		PIPsoft = list(threshold=extra.shrinkage, lfdr=matrix(NA, nrow=m, ncol=r), u=matrix(NA, nrow=m, ncol=r), var=vector("numeric", r), PVE=vector("numeric", r))
		for(i in 1:r) {
			if(PIPsoft$threshold[i]>0) PIPsoft$lfdr[,i] = sapply(PIP$lfdr[,1], function(x) (PIPsoft$threshold[i] - max((PIPsoft$threshold[i] - x), 0)) / PIPsoft$threshold[i])
			if(PIPsoft$threshold[i]==0) PIPsoft$lfdr[,i] = rep(1, m)
			PIPsoft$u[,i] = (1-PIPsoft$lfdr[,i]) * svd.raw$u[,i]
			PIPsoft$var[i] = svd.raw$d[i]^2 * sum(PIPsoft$u[,i]^2)
			PIPsoft$PVE[i] = PIPsoft$var[i] / var.total
		}
		out = c(out, list(PIPsoft=PIPsoft))
	}

    return(out)
}
