kerndwd = function(x, y, kern, qval=1, lambda, wt=NULL,
  eps=1e-05, maxit=1e+05) {
  ####################################################################
  ## data setup
  this.call = match.call()
  y = c(-1, 1)[as.factor(drop(y))]
  x = as.matrix(x)
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels.")
  if (qval <= 0) 
    stop("qval must be positive.")  
  nobs = as.integer(NROW(x))
  np = as.integer(NCOL(x))
  if (length(y) != nobs) 
    stop("x and y have different number of observations.")
  if (missing(kern)) {
    kern = rbfdot(sigma=1)
    cat("kern is missing: Gaussian kernel is used.\n")
  }  
  ####################################################################
  qval = as.double(qval)
  maxit = as.integer(maxit)
  eps = as.double(eps)
  gam = as.double(1e-8) # a tiny value to avoid matrix sigularity
  ####################################################################
  ## lambda setup
  if (is.null(lambda)) {
    stop("Users have to provide a lambda sequence.")
  } else {
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  ####################################################################
  ## call Fortran core
  if (is.null(wt)) {
    if (class(kern)[[1]] == "vanillakernel" && nobs > np) {
    ## linear DWD
      ftran = "ldwd"
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(x), 
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], 
        np + 1, anlam) 
    } else {
    ## kernel DWD
      ftran = "kdwd"
      Kmat = kernelMatrix(kern, x)
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(Kmat), 
        nobs, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], 
        nobs + 1, anlam) 
    }
  } else {
    if (length(wt) != nobs)
      stop("The length of the weight vector is not n.")
    if (class(kern)[[1]] == "vanillakernel" && nobs > np) {
    ## weighted linear DWD
      ftran = "wldwd"
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(x), as.double(wt),
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], 
        np + 1, anlam) 
    } else {
    ## weighted kernel DWD
      ftran = "wkdwd"
      Kmat = kernelMatrix(kern, x)
      if(abs(qval - as.integer(qval)) < .Machine$double.eps) {
        ftran = paste0(ftran, "int")
        qval = as.integer(qval)
      }
      fit = .Fortran(ftran, qval, as.double(Kmat), as.double(wt),
        nobs, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="kerndwd")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], 
        nobs + 1, anlam)  
    } 
  }
  ####################################################################
  ## wrap up output
  info = list(qval = qval, eps = eps, maxit = signif(maxit),
    kern = capture.output(show(kern)))
  if (!is.null(wt)) info = c(info, wt=list(wt))
  outlist = list(alpha = alpha, lambda = ulam[seq(anlam)], 
    npass = fit$npass, jerr = fit$jerr, info = info,
    call = this.call)
  class(outlist) = c("kerndwd")
  outlist
} 
