\name{tunedwd}
\alias{tunedwd}
\title{fast tune procedure for DWD}
\description{A fast implementaiton of cross-validation for \code{\link{kerndwd}} to find the optimal values of the tuning parameter \code{lambda}.}
\usage{
tunedwd(x, y, kern, lambda, qvals=1, eps=1e-5, maxit=1e+5, nfolds=5, foldid=NULL)
}
\arguments{
    \item{x}{A matrix of predictors, i.e., the matrix \code{x} used in \code{\link{kerndwd}}.}
    \item{y}{A vector of binary class labels, i.e., the \code{y} used in \code{\link{kerndwd}}. \code{y} has two levels.}
    \item{kern}{A kernel function.}
    \item{lambda}{A user specified \code{lambda} candidate sequence for cross-validation.}
    \item{qvals}{A vector containing the index of the generalized DWD. Default value is 1.}
    \item{eps}{The algorithm stops when (i.e. \eqn{\sum_j(\beta_j^{new}-\beta_j^{old})^2}{sum(j)(beta_new[j]-beta_old[j])^2} is less than \code{eps}, where \eqn{j=0,\ldots, p}. Default value is \code{1e-5}.}
    \item{maxit}{The maximum of iterations allowed. Default is 1e5.}
    \item{nfolds}{The number of folds. Default value is 5. The allowable range is from 3 to the sample size.}
    \item{foldid}{An optional vector with values between 1 and \code{nfold}, representing the fold indices for each observation. If supplied, \code{nfold} can be missing.}
}


\details{This function returns the best tuning parameters
\code{q} and \code{lambda} by cross-validation. An efficient tune method is employed to accelerate the algorithm.
}

\value{
A \code{tunedwd.kerndwd} object including the cross-validation results is return.
    \item{lam.tune}{The optimal \code{lambda} value.}
    \item{q.tune}{The optimal \code{q} value.}
}


\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang-wang@uiowa.edu}}
\references{
Wang, B. and Zou, H. (2018)
``Another Look at Distance Weighted Discrimination," 
\emph{Journal of Royal Statistical Society, Series B}, \bold{80}(1), 177--198. \cr
\url{https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12244}\cr
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22.\cr
\url{https://www.jstatsoft.org/v33/i01/paper}
}

\seealso{\code{\link{kerndwd}}.}
\examples{
set.seed(1)
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
kern = rbfdot(sigma=sigest(BUPA$X))
ret = tunedwd(BUPA$X, BUPA$y, kern, qvals=c(1,2,10), lambda=lambda, eps=1e-5, maxit=1e5)
ret
}

