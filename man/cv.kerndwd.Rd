\name{cv.kerndwd}
\alias{cv.kerndwd}
\title{cross-validation}
\description{Carry out a cross-validation for \code{\link{kerndwd}} to find optimal values of the tuning parameter \code{lambda}.}
\usage{
cv.kerndwd(x, y, kern, lambda, nfolds=5, foldid, wt, ...)
}
\arguments{
    \item{x}{A matrix of predictors, i.e., the matrix \code{x} used in \code{\link{kerndwd}}.}
    \item{y}{A vector of binary class labels, i.e., the \code{y} used in \code{\link{kerndwd}}. \code{y} has to be two levels.}
    \item{kern}{A kernel function.}
    \item{lambda}{A user specified \code{lambda} candidate sequence for cross-validation.}
    \item{nfolds}{The number of folds. Default value is 5. The allowable range is from 3 to the sample size.}
    \item{foldid}{An optional vector with values between 1 and \code{nfold}, representing the fold indices for each observation. If supplied, \code{nfold} can be missing.}
     \item{wt}{A vector of length \eqn{n} for weight factors. When \code{wt} is missing or \code{wt=NULL}, an unweighted DWD is fitted. } 
    \item{\dots}{Other arguments being passed to \code{\link{kerndwd}}.}
}

\details{This function computes the mean cross-validation error and the standard error by fitting \code{\link{kerndwd}} with every fold excluded alternatively. This function is modified based on the \code{cv} function from the \code{glmnet} package.
}
\value{A \code{\link{cv.kerndwd}} object including the cross-validation results is return..
    \item{lambda}{The \code{lambda} sequence used in \code{\link{kerndwd}}.}
    \item{cvm}{A vector of length \code{length(lambda)}: mean cross-validated error.}
    \item{cvsd}{A vector of length \code{length(lambda)}: estimates of standard error of \code{cvm}.}
    \item{cvupper}{The upper curve: \code{cvm + cvsd}.}
    \item{cvlower}{The lower curve: \code{cvm - cvsd}.}
    \item{lambda.min}{The \code{lambda} incurring the minimum cross validation error \code{cvm}.}
    \item{lambda.1se}{The largest value of \code{lambda} such that error is within one standard error of the minimum.}
    \item{cvm.min}{The cross-validation error corresponding to \code{lambda.min}, i.e., the least error.}
    \item{cvm.1se}{The cross-validation error corresponding to \code{lambda.1se}.}
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

\seealso{\code{\link{kerndwd}} and \code{\link{plot.cv.kerndwd}}}
\examples{
set.seed(1)
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(3, -3, length.out=10))
kern = rbfdot(sigma=sigest(BUPA$X))
m.cv = cv.kerndwd(BUPA$X, BUPA$y, kern, qval=1, lambda=lambda, eps=1e-5, maxit=1e5)
m.cv$lambda.min
}

