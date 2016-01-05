#' The Jackstraw Weighted Shrinkage Estimation Method for Sparse Coefficients in Logistic Factor Analysis
#'
#' Estimates sparse/shrunken coefficients (or loadings) of Logistic Factor Analysis (LFA).
#' Based on statistical sginificance of association between variables (i.e., SNPs) and logistic factors (LFs),
#' the conventional coefficients of LFs are shruken towards zeros, which improve its accuracy.
#' The only required inputs are the data matrix \code{dat} and the number of LFs \code{r}.
#'
#' By default, \code{jaws.lfa} computes two canonical jackstraw weighted shrinkage estimators, namely \code{PIP} and \code{PNV}.
#'
#' It is strongly advised that you take a careful look at your data and use appropriate graphical and statistical criteria
#' to determine \code{r} (Note that \code{r} does not include an intercept term).
#'
#' If \code{s} is not supplied, \code{s} is set to about 10\% of \code{m} variables.
#' If \code{B} is not supplied, \code{B} is set to \code{m*10/s}.
#'
#' NOTE that you may input a \code{m*r} matrix of p-values obtained from running \code{jackstraw.LFA}.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param p \code{m*r} matrix of p-values obtained from running \code{jackstraw.LFA}.
#' @param r a number (a positive integer) of significance logistic factors.
#' @param s a number (a positive integer) of ``synthetic'' null variables (optional).
#' @param B a number (a positive integer) of resampling iterations (optional).
#' @param verbose a logical specifying to print the progress (default: TRUE).
#' @param save.all a logical specifying to save all objects, including a large SVD object (default: FALSE).
#' @param seed a seed for the random number generator (optional).
#'
#' @return \code{jaws.lfa} returns a list consisting of
#' \item{p}{p-values for association tests between variables and each of \code{r} logistic factors}
#' \item{pi0}{proportion of variables not associated with \code{r} logistic factors, individually}
#' \item{PIP}{a list of outputs derived from the posterior inclusion probabilities method (including \code{pr}, \code{coefs}, \code{af}, \code{r2})}
#' \item{PNV}{a list of outputs derived from the proportion of null variables method (including \code{pi0}, \code{coefs}, \code{af}, \code{r2})}
#'
#' @section Detailed explanation of the output contained in \code{PIP} and \code{PNV}:
#' \describe{
#'   \item{pr}{a matrix of posterior inclusion probabilities (equivalent to 1-lfdr) for \code{m} coefficients and \code{r} LFs.}
#'   \item{pi0}{a vector of estimated proportion of null variables for \code{r} LFs.}
#'   \item{coefs}{a \code{m*r} matrix of shrunken coefficients.}
#'   \item{af}{a \code{m*n} matrix of shrunken allele frequencies.}
#'   \item{r2}{a vector of shrunken McFadden's pseudo R^2 measures for \code{r} LFs.}
#' }
#'
#' @export jaws.lfa
#' @import lfa
#' @importFrom corpcor fast.svd
#' @import jackstraw
#' @importFrom qvalue pi0est lfdr qvalue
#' @author Neo Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Forthcoming
#'
#' @seealso \link{jackstraw.LFA} \link{jaws.pca}
jaws.lfa = function(dat, p=NULL, r=NULL, s=NULL, B=NULL, verbose=TRUE, seed=NULL, save.all=FALSE) {
  	if(!is.null(seed)) set.seed(seed)
	m=dim(dat)[1]; n=dim(dat)[2]

	if(is.null(r)) {
		stop("The number of LFs (r) is missing.")
	}
	if(!(r > 0 && r < n)) { stop("A number of significant LFs is not in valid range between 1 and n-1."); }

	if(!is.null(p)) p = as.matrix(p)
	# compute p-values
	if(is.null(p)) {
		p = matrix(NA, nrow=m, ncol=r)
		if(verbose==TRUE) message("\nThe association significance between SNPs and logistic factors are computed by the jackstraw.")
		for(i in 1:r) {
			p[,i] = jackstraw.LFA(dat, r1=i, r=r, s=s, B=B, verbose=verbose)$p.value
		}
	}
	if(ncol(p) != r | nrow(p) != m) stop("The p-value matrix must match m (nrow) and r (ncol).")

	# compute pi0
	# averaging two methods to calculate pi0 seems to be stable and accurate (esp. at a very high pi0)
	pi0 = matrix(0, nrow=2, ncol=r)
	for(i in 1:r) {
		pi0[1,i] = tryCatch(qvalue(p[,i], pi0.method="smoother")$pi0, error=function(tmp) NA)
		pi0[2,i] = tryCatch(qvalue(p[,i], pi0.method="bootstrap")$pi0, error=function(tmp) NA)
	}
	pi0 = apply(pi0, 2, function(x) mean(x, na.rm=TRUE))

	# compute LFA and coefficients
	LF = lfa(dat, r)
	coefs = lfa_full(dat, LF)$coefs
	AF_null = af(dat, matrix(1, n, 1))
	coefs_marginal = lfa_full(dat, LF[,r+1,drop=FALSE])$coefs

	##############################################################################################################################
 	## Posterior Inclusion Probabilities: Empirical Bayes Shrinkage
	if(verbose==TRUE) message(paste0("Computing the PIP Shrunken Loadings for the LF (r=", r, ") : "))
	PIP = list(pr=matrix(NA, nrow=m, ncol=r), coefs=coefs)
	for(i in 1:r) {
		PIP$pr[,i] = 1-lfdr(p[,i], pi0=pi0[i])
		PIP$coefs[,i] = PIP$pr[,i] * PIP$coefs[,i]
		if(verbose==TRUE) cat(paste(i," "))
	}
	PIP$af = af_coefs(dat, LF, PIP$coefs)
	PIP$r2.m = mcfadden_Rsq_AF(dat, AF_alt=PIP$af, AF_null=AF_null)
	PIP$r2.e = efron_Rsq_AF(dat, PIP$af)

##	setting the coefficnets for the intercept term to marginal coefficients
	PIP$coefs[,r+1] = coefs_marginal
	PIP$af = af_coefs(dat, LF, PIP$coefs)
	PIP$r2.m.2 = mcfadden_Rsq_AF(dat, AF_alt=PIP$af, AF_null=AF_null)
	PIP$r2.e.2 = efron_Rsq_AF(dat, PIP$af)

 	## proportion of null variables: Pi0 Hard Thresholding
	if(verbose==TRUE) message(paste0("Computing the PNV Shrunken Loadings for the LF (r=", r, ") : "))
	PNV = list(pi0=pi0, coefs=coefs)
	for(i in 1:r) {
		PNV$coefs[which(rank(abs(PNV$coefs[,i])) <= m*pi0[i]),i] = 0	# <= choice ensures a conservative sparsity in a case of ties
		if(verbose==TRUE) cat(paste(i," "))
	}
	PNV$af = af_coefs(dat, LF, PNV$coefs)
	PNV$r2.m = mcfadden_Rsq_AF(dat, AF_alt=PNV$af, AF_null=AF_null)
	PNV$r2.e = efron_Rsq_AF(dat, PNV$af)

##	setting the coefficnets for the intercept term to marginal coefficients
	PNV$coefs[,r+1] = coefs_marginal
	PNV$af = af_coefs(dat, LF, PNV$coefs)
	PNV$r2.m.2 = mcfadden_Rsq_AF(dat, AF_alt=PNV$af, AF_null=AF_null)
	PNV$r2.e.2 = efron_Rsq_AF(dat, PNV$af)

	out = list(call=match.call(), LF=LF, coefs=coefs, p=p, pi0=pi0, PIP=PIP, PNV=PNV)

	if(save.all==FALSE) {
		PIP$af = NULL
		PNV$af = NULL
		out = list(call=match.call(), p=p, pi0=pi0, PIP=PIP, PNV=PNV)
	}

    return(out)
}
