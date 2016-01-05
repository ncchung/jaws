## dat must be row-centered
r2.est = function(dat, est) {
	m = dim(dat)[1]
	n = dim(dat)[2]

	means = rowMeans(dat)
	SSt = sapply(1:m, function(x) sum((dat[x,] - means[x])^2))

	res = dat - est
	SSr = sapply(1:m, function(x) sum((res[x,])^2))

	return(1 - SSr/SSt)
}

## include an intercept term in mod
r2.mod = function(dat, mod) {
	m = dim(dat)[1]
	n = dim(dat)[2]
	means = rowMeans(dat)
	SSt = sapply(1:m, function(x) sum((dat[x,] - means[x])^2))

	B = dat %*% mod %*% solve(t(mod) %*% mod)
	H = mod %*% solve(t(mod) %*% mod) %*% t(mod)
	fit = dat %*% H
	res = dat - fit
	SSr = sapply(1:m, function(x) sum((res[x,])^2))

	return(1 - SSr/SSt)
}

## mcfadden_Rsq_snp() in the jackstraw package.
## af_snp() in the lfa package.
mcfadden_Rsq_AF = function(X, AF_alt, AF_null=NULL){
    m=nrow(X); n=ncol(X)
    if(is.null(AF_null)){
        LF_null = matrix(1, n, 1)
        AF_null = t(apply(X, 1, lfa::af_snp, LF_null))
    }

    Rsq=sapply(1:m, function(i){jackstraw:::mcfadden_Rsq_snp(X[i,], AF_alt[i,], AF_null[i,])})
    return(Rsq)
}

## efron_Rsq_snp() in the jackstraw package.
efron_Rsq_AF = function(X, AF){
    m=nrow(X);n=ncol(X)

    Rsq=sapply(1:m, function(i){jackstraw:::efron_Rsq_snp(X[i,], AF[i,])})
    return(Rsq)
}

## Similar to af() in the lfa package, but it can accept shrunken coefficients.
af_coefs = function(dat, LF, coefs){
    m=nrow(dat)
    coefs = t(sapply(1:m, function(i){af_snp_coefs(snp=dat[i, ,drop=FALSE], LF=LF, coefs=coefs[i, ,drop=FALSE])}))

    return(coefs)
}

## Similar to af_snp() in the lfa package, but it can accept shrunken coefficients.
## lreg() in the lfa package.
af_snp_coefs = function(snp, LF, coefs=NULL){
	if(is.null(coefs)) {
		#coefficients from logistic regression assuming HWE
		coefs = lfa:::lreg(snp,LF)
	}
    est = LF %*% t(coefs)

    return(exp(est)/(1+exp(est)))
}
