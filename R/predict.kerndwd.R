predict.kerndwd = function(object, kern, x, newx,
    type=c("class", "link"), ...) {
  type = match.arg(type)
  if (class(kern)[[1]] == "vanillakernel" && NROW(x) > NCOL(x))
    nfit = newx %*% object$alpha[-1, ] else 
    nfit = kernelMult(kern, newx, x, object$alpha[-1, ])
  nfit = sweep(nfit, MARGIN = 2, object$alpha[1, ], "+")
  switch(type, link = nfit, class = ifelse(nfit > 0, 1, -1))
} 
