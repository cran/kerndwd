kerndwd = function(x, y, kern, lambda, qval=1, wt=NULL, eps=1e-05, maxit=1e+05) {
  ####################################################################
  ## data setup
  this.call = match.call()
  if (length(levels(factor(y))) == 2)
    y = c(-1, 1)[as.factor(drop(y))]
  x = as.matrix(x)
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels.")
  nobs = as.integer(NROW(x))
  np = as.integer(NCOL(x))
  if (length(y) != nobs) 
    stop("x and y have different number of observations.")
  if (missing(kern)) {
    kern = rbfdot(sigma=sigest(x))
    cat("'kern' is missing: Gaussian kernel is used.\n")
  }
  ####################################################################
  # qval = as.double(qval)
  maxit = as.integer(maxit)
  eps = as.double(eps)
  gam = as.double(1e-7) # a tiny value to avoid matrix sigularity
  ####################################################################
  ## lambda setup
  if (missing(lambda)) {
    if (nobs <= 1e4) 
      lambda = 10 ^ seq(3, -3, len=100) else
      lambda = 10 ^ seq(0, -7, len=100)
    ulam = lambda
  } else {
    ulam = as.double(sort(lambda, decreasing=TRUE))
  }
  nlam = as.integer(length(lambda))
  fit = dwdpath(x, y, nobs, np, kern, qval, ulam, nlam, wt, eps, maxit, gam)
  fit.call = this.call
  class(fit) = c(class(fit), "kerndwd")
  fit
} 

dwdpath = function(x, y, nobs, np, kern, qval, ulam, nlam, wt, eps, maxit, gam) {
  #################################################################### 
  ## check index 
  if (qval <= 0) {
    warning("The parameter 'qval' must be positive; set to 1.")
    qval = 1
  }  
  ####################################################################
  if (missing(wt) || is.null(wt)) {
    if (class(kern)[[1]] == "vanillakernel" && nobs >= np) {
    ## linear DWD
      fit = .Fortran("ldwd", qval, as.double(x), 
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], np + 1, anlam) 
    } else {
    ## kernel DWD
      Kmat = kernelMatrix(kern, x)
      fit = .Fortran("kdwd", qval, as.double(Kmat), 
        nobs, as.double(y), nlam, ulam, eps, maxit, gam, 
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], nobs + 1, anlam) 
    }
  } else {
    if (length(wt) != nobs)
      stop("The length of the weight vector is not n.")
    if (any(wt < 0))
      stop("The weights must be nonnegative.")
    if (class(kern)[[1]] == "vanillakernel" && nobs > np) {
    ## weighted linear DWD
      fit = .Fortran("wldwd", qval, as.double(x), as.double(wt),
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], np + 1, anlam) 
    } else {
    ## weighted kernel DWD
      Kmat = kernelMatrix(kern, x)
      fit = .Fortran("wkdwd", qval, as.double(Kmat), as.double(wt),
        nobs, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], nobs + 1, anlam)  
    } 
  }
  ####################################################################
  ## wrap up output
  info = list(qval = qval, eps = eps, maxit = signif(maxit),
    kern = capture.output(show(kern)))
  if (!missing(wt)) info = c(info, wt=list(wt))
  outlist = list(alpha = alpha, lambda = ulam[seq(anlam)], 
    npass = fit$npass, jerr = fit$jerr, info = info)
  class(outlist) = c("dwdpath")
  outlist
}
