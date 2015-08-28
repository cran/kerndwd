\name{cv.kerndwd}
\alias{cv.kerndwd}
\title{cross-validation for DWD}
\description{Carry out a cross-validation for \code{\link{kerndwd}} to find the optimal value of the tuning parameter \code{lambda}.}
\usage{
cv.kerndwd(x, y, kern, lambda, qval=1, 
  pred.loss=c("misclass", "loss"), nfolds=5, foldid, ...)
}
\arguments{
		\item{x}{A matrix of predictors, i.e., the matrix \code{x} used in \code{\link{kerndwd}}.}
		\item{y}{A vector of binary class labels, i.e., the \code{y} used in \code{\link{kerndwd}}.}
		\item{kern}{A kernel function.}
    \item{lambda}{A user specified \code{lambda} candidate sequence for cross-validation.}
    \item{qval}{The index of the generalized DWD. Default value is 1.}
		\item{pred.loss}{If the value is "loss", the DWD loss is used as the cross-validation criterion; if the value is "misclass" or else, mis-classification rate is the criterion.}
    \item{nfolds}{The number of folds. Default value is 5. The allowable range is from 3 to the sample size.}
		\item{foldid}{An optional vector with values between 1 and \code{nfold}, representing the folder indices for each observation.  
If supplied, \code{nfold} can be missing.}
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
		\item{name}{``mis-classification error" or "DWD loss", depending on \code{pred.loss}, for plotting purposes.}
		\item{lambda.min}{The \code{lambda} incurring the minimum cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is within one standard error of the minimum.}
    \item{cvm.min}{The cross-validation error corresponding to \code{lambda.min}, i.e., the least error.}
		\item{cvm.1se}{The cross-validation error corresponding to \code{lambda.1se}.}
}

\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang@umn.edu}}
\references{
Wang, B. and Zou, H. (2015)
``Another Look at DWD: Thrifty Algorithm and Bayes Risk Consistency in RKHS". \cr
\url{http://arxiv.org/abs/1508.05913v1.pdf}\cr
\cr
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22\cr
\url{http://www.jstatsoft.org/v33/i01/paper}
\cr
}

\seealso{\code{\link{kerndwd}} and \code{\link{plot.cv.kerndwd}}}
\examples{
set.seed(1)
data(Haberman)
Haberman$X = scale(Haberman$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
kern = rbfdot(sigma=1)
m.cv = cv.kerndwd(Haberman$X, Haberman$y, kern, qval=1, 
  lambda=lambda, eps=1e-5, maxit=1e5)
m.cv$lambda.min
}

