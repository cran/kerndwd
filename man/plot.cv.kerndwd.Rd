\name{plot.cv.kerndwd}
\alias{plot.cv.kerndwd}
\title{plot the cross-validation curve}
\description{
Plot cross-validation error curves with the upper and lower standard deviations versus log \code{lambda} values.}
\usage{
\method{plot}{cv.kerndwd}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{A fitted \code{\link{cv.kerndwd}} object.}
	\item{sign.lambda}{Against \code{log(lambda)} (default) or its negative if \code{sign.lambda=-1}.}
	\item{\dots}{Other graphical parameters being passed to \code{plot}.}
}
\details{This function plots the cross-validation error curves. This function is modified based on the \code{plot.cv} function of the \code{glmnet} package.
}

\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang@umn.edu}}
\references{
Wang, B. and Zou, H. (2015)
``Another Look at DWD: Thrifty Algorithm and Bayes Risk Consistency in RKHS". \cr
\url{http://arxiv.org/abs/1508.05913v1.pdf}\cr
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22\cr
\url{http://www.jstatsoft.org/v33/i01/paper}
}
 
\seealso{\code{\link{cv.kerndwd}}.}
\examples{
set.seed(1)
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
kern = rbfdot(sigma=1)
m.cv = cv.kerndwd(BUPA$X, BUPA$y, kern, family="DWD", 
  qval=1, lambda=lambda, eps=1e-5, maxit=1e5)
m.cv
}
