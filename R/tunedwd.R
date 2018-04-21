tunedwd = function(x, y, kern, lambda, qvals=1, eps=1e-5, maxit=1e+5, 
  nfolds=5, foldid=NULL) {
  ####################################################################
  ## data setup
  this.call = match.call()
  y = c(-1, 1)[as.factor(drop(y))]
  x = as.matrix(x)
  nobs = as.integer(NROW(x))
  np = as.integer(NCOL(x))
  if (length(y) != nobs) 
    stop("x and y have different number of observations.")
  if (missing(kern)) {
    kern = rbfdot(sigma=sigest(x))
    cat("'kern' is missing: Gaussian kernel is used.\n")
  }
  maxit = as.integer(maxit)
  eps = as.double(eps)
  gam = as.double(1e-6)
  if (missing(lambda)) {
    if (nobs <= 1e4) 
      lambda = 10 ^ seq(3, -3, len=100) else
      lambda = 10 ^ seq(0, -7, len=100)
    ulam = lambda
  } else {
    ulam = as.double(sort(lambda, decreasing=TRUE))
  } 
  nlam = as.integer(length(lambda))
  qlen = length(qvals)
  ####################################################################
  ## fit DWD 
  if (is.null(foldid)) 
    foldid = sample(rep(seq(nfolds), length=nobs)) else 
    nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be larger than 3; nfolds=5 recommended.")
  predmat = array(NA, c(nobs, nlam, qlen))
  is.linr = (class(kern)[[1]] == "vanillakernel" && nobs > np)
  nlams = matrix(NA, qlen, nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    x.train = x[!which, , drop=FALSE]
    y.train = y[!which]
    nobs.i = length(y.train)
    x.test = drop(as.matrix(x[which, , drop=FALSE]))
    if (is.null(dim(x.test))) dim(x.test) = c(1, np)
    if (is.linr) {
      fit = .Fortran("lqdwd", qvals, qlen, as.double(x.train), 
        nobs.i, np, as.double(y.train), nlam, ulam, eps, maxit, gam,
        anlam=integer(qlen), npass=integer(nlam*qlen), jerr=integer(qlen), 
        btmat=double((np+1)*nlam*qlen), PACKAGE="kerndwd")
      alpha = array(fit$btmat, c(np+1, nlam, qlen))
    } else {
      Kmat = kernelMatrix(kern, x.train)
      fit = .Fortran("kqdwd", qvals, qlen, as.double(Kmat), 
        nobs.i, as.double(y.train), nlam, ulam, eps, maxit, gam, 
        anlam=integer(qlen), npass=integer(nlam*qlen), jerr=integer(qlen), 
        alpmat=double((nobs.i+1)*nlam*qlen), PACKAGE="kerndwd")
      alpha = array(fit$alpmat, c(nobs.i+1, nlam, qlen))
      
    }
    nlams[,i] = fit$anlam
    for (qq in seq(qlen)) {
      if (is.linr) {
        alp = matrix(alpha[,seq(fit$anlam[qq]),qq], np+1, fit$anlam[qq])
        nfit = x.test %*% alp[-1, , drop=FALSE]
      } else {
        alp = matrix(alpha[,seq(fit$anlam[qq]),qq], nobs.i+1, fit$anlam[qq])
        nfit = kernelMult(kern, x.test, x.train, alp[-1,, drop=FALSE])
      }
      nfit = sweep(nfit, MARGIN=2, alp[1,,drop=FALSE], "+")
      predmat[which, seq(nlam), qq] = nfit 
    }
  }
  nlams = rep(nlam, nfolds)
  cvres = as.list(qlen)
  cvms = rep(NA, qlen)
  for (qq in seq(qlen)) {
    cvraw = (y != ifelse(predmat[,,qq] > 0, 1, -1))
    if (length(y)/nfolds >= 3) {
      cvob = cvcompute(cvraw, foldid, nlams)
      cvraw = cvob$cvraw
      cvn = cvob$cvn
    } else cvn = nobs - colSums(is.na(predmat[,,qq]))    
    cvm = colMeans(cvraw, na.rm=TRUE)
    cvsd = sqrt(colMeans(scale(cvraw, cvm, FALSE)^2, na.rm=TRUE)/(cvn - 1))
    cvres[[qq]] = list(qval=qvals[qq], lambda=ulam, cvm=cvm, cvsd=cvsd, 
      cvup=cvm+cvsd, cvlo=cvm-cvsd, predmat=predmat[,,qq])
    cvres[[qq]] = c(cvres[[qq]], as.list(getmin(ulam, cvm, cvsd)))
    cvms[qq] = cvres[[qq]]$cvm.1se
  }
  q.index = qlen + 1 - which.min(rev(cvms))
  #q.index = which.min(cvms)
  q.tune = qvals[q.index]
  lam.tune = cvres[[q.index]]$lambda.1se
  obj = list(lam.tune=lam.tune, q.tune=q.tune, cvres=cvres)
  class(obj) = "tunedwd.kerndwd"
  obj
} 
