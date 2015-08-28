cv.kerndwd = function(x, y, kern, lambda, qval=1,
  pred.loss=c("misclass", "loss"), nfolds=5, foldid, ...) {
  ####################################################################
  ## data setup
  if (missing(pred.loss)) 
    pred.loss = "misclass" else pred.loss = match.arg(pred.loss)
  if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
    warning("Only 'misclass' and 'loss' available for 
            DWD; 'misclass' used.")
    pred.loss = "misclass"
  }
  typenames = c(misclass = "mis-classification error",
                loss     = "DWD loss")
  y = c(-1, 1)[as.factor(drop(y))]
  x = as.matrix(x)
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels.")
  x.row = as.integer(NROW(x))
  if (length(y) != x.row) 
    stop("x and y have different number of observations.")  
  ####################################################################
  ## fit DWD nfold times and save them 
  if (missing(foldid)) 
    foldid = sample(rep(seq(nfolds), length=x.row)) else 
    nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be larger than 3; nfolds=5 recommended.")
  predmat = matrix(NA, x.row, length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = kerndwd(x=x[!which, , drop=FALSE], 
      y=y[!which], kern=kern, lambda=lambda, qval=qval, ...)
    preds = predict(fitobj, kern, x[!which, , drop=FALSE], 
      x[which, , drop=FALSE], type="link")
    nlami = length(fitobj$lambda)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  ## select the lambda according to predmat
  cvraw = switch(pred.loss, 
    loss     = dwdloss(y * predmat, qval),
    misclass = (y != ifelse(predmat > 0, 1, -1)))
  if (length(y)/nfolds >= 3) {
    cvob = cvcompute(cvraw, foldid, nlams)
    cvraw = cvob$cvraw
    cvn = cvob$N
  } else cvn = length(y) - colSums(is.na(predmat))    
  cvm = colMeans(cvraw, na.rm=TRUE)
  cvsd = sqrt(colMeans(scale(cvraw, cvm, FALSE)^2, 
    na.rm=TRUE)/(cvn - 1))
  ## wrap up output
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
    cvup = cvm + cvsd, cvlo = cvm - cvsd, 
    name = as.character(typenames[pred.loss]))
  obj = c(out, as.list(getmin(lambda, cvm, cvsd)))
  class(obj) = "cv.kerndwd"
  obj
} 

