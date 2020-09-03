\name{predict.kerndwd}
\alias{predict.kerndwd}
\title{predict class labels for new observations}
\description{
Predict the binary class labels or the fitted values of an \code{\link{kerndwd}} object.
}

\usage{
\method{predict}{kerndwd}(object, kern, x, newx, type=c("class", "link"), ...)
}

\arguments{
		\item{object}{A fitted \code{\link{kerndwd}} object.}
    \item{kern}{The kernel function used when fitting the \code{\link{kerndwd}} object.}
    \item{x}{The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{\link{kerndwd}} object.}
		\item{newx}{A matrix of new values for \code{x} at which predictions are to be made. We note that \code{newx} must be a matrix, \code{predict} function does not accept a vector or other formats of \code{newx}.}
		\item{type}{\code{"class"} or \code{"link"}? \code{"class"} produces the predicted binary class labels and \code{"link"} returns the fitted values. Default is \code{"class"}.}
    \item{\dots}{Not used. Other arguments to \code{predict}.} 
}
 
\details{If \code{"type"} is \code{"class"}, the function returns the predicted class labels. If \code{"type"} is \code{"link"}, the result is \eqn{\beta_0 + x_i'\beta} for the linear case and \eqn{\beta_0 + K_i'\alpha} for the kernel case.}

\value{Returns either the predicted class labels or the fitted values, depending on the choice of \code{type}.}

\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang-wang@uiowa.edu}}
\references{
Wang, B. and Zou, H. (2018)
``Another Look at Distance Weighted Discrimination," 
\emph{Journal of Royal Statistical Society, Series B}, \bold{80}(1), 177--198. \cr
\url{https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12244}\cr
}
 
\seealso{\code{\link{kerndwd}}}
\examples{
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
kern = rbfdot(sigma=sigest(BUPA$X))
m1 = kerndwd(BUPA$X, BUPA$y, kern,
  qval=1, lambda=lambda, eps=1e-5, maxit=1e5)
predict(m1, kern, BUPA$X, tail(BUPA$X))
}

