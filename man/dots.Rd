\name{kernel functions}
\alias{kern}
\alias{dots}
\alias{rbfdot}
\alias{polydot}
\alias{vanilladot}
\alias{laplacedot}
\alias{besseldot}
\alias{anovadot}
\alias{splinedot}
\alias{rbfkernel-class}
\alias{polykernel-class}
\alias{vanillakernel-class}
\alias{anovakernel-class}
\alias{besselkernel-class}
\alias{laplacekernel-class}
\alias{splinekernel-class}
\alias{sigest}
\title{Kernel Functions}
\description{
  Kernel functions provided in the R package \code{kernlab}. Details can be seen in the reference below.\cr
  The Gaussian RBF kernel \eqn{k(x,x') = \exp(-\sigma \|x - x'\|^2)}{k(x, x') = exp(-\sigma|x - x'|^2)} \cr
  The Polynomial kernel \eqn{k(x,x') = (scale <x, x'> + offset)^{degree}}{k(x,x') = (scale <x, x'> + offset)^degree}\cr
  The Linear kernel \eqn{k(x,x') = <x, x'>}{k(x, x') = <x, x'>}\cr
  The Laplacian kernel \eqn{k(x,x') = \exp(-\sigma \|x - x'\|)}{k(x, x') = exp(-\sigma|x - x'|)} \cr
  The Bessel kernel \eqn{k(x,x') = (- \mathrm{Bessel}_{(\nu+1)}^n \sigma \|x - x'\|^2)}{k(x, x') = (-Bessel{(\nu+1)}^n \sigma|x - x'|^2)} \cr
  The ANOVA RBF kernel \eqn{k(x,x') = \sum_{1\leq i_1 \ldots < i_D \leq N} 
  \prod_{d=1}^D k(x_{id}, {x'}_{id})}{k(x, x') = sum_{1 <= i_1, ..., < i_D <= N} 
  prod_{d=1}^D k(x_id, x'_id)} where k(x, x) is a Gaussian RBF kernel. \cr
  The Spline kernel \eqn{ \prod_{d=1}^D 1 + x_i x_j + x_i x_j \min(x_i,
    x_j)  - \frac{x_i + x_j}{2} \min(x_i,x_j)^2 +
    \frac{\min(x_i,x_j)^3}{3}}{prod_{d=1}^D 1 + x_i x_j + x_i x_j min(x_i,
    x_j)  - (x_i + x_j)/2 min(x_i,x_j)^2 +
    (min(x_i,x_j)^3)/3}.
  The parameter \code{sigma} used in \code{rbfdot} can be selected by \code{sigest()}.
}
\usage{
rbfdot(sigma = 1)
polydot(degree = 1, scale = 1, offset = 1)
vanilladot()
laplacedot(sigma = 1)
besseldot(sigma = 1, order = 1, degree = 1)
anovadot(sigma = 1, degree = 1)
splinedot()
sigest(x)
}

\arguments{
  \item{sigma}{The inverse kernel width used by the Gaussian, the
    Laplacian, the Bessel, and the ANOVA kernel.}
  \item{degree}{The degree of the polynomial, bessel or ANOVA
    kernel function. This has to be an positive integer.}
  \item{scale}{The scaling parameter of the polynomial kernel function.}
  \item{offset}{The offset used in a polynomial kernel.}
  \item{order}{The order of the Bessel function to be used as a kernel.}
  \item{x}{The design matrix used in \code{kerndwd} when \code{sigest} is called to estimate \code{sigma} in \code{rbfdot()}.}
}
\details{These R functions and descriptions are directly duplicated and/or adapted from the R package \code{kernlab}.}  

\value{
 Return an S4 object of class \code{kernel} which can be used as the argument of \code{kern} when fitting a \code{\link{kerndwd}} model.
}

\references{
Wang, B. and Zou, H. (2017)
``Another Look at Distance Weighted Discrimination," \emph{Journal of Royal Statistical Society, Series B}, to appear. \cr
\url{http://arxiv.org/abs/1508.05913}\cr

Karatzoglou, A., Smola, A., Hornik, K., and Zeileis, A. (2004)
``kernlab -- An S4 Package for Kernel Methods in R", 
\emph{Journal of Statistical Software}, \bold{11}(9), 1--20.\cr
\url{http://www.jstatsoft.org/v11/i09/paper}\cr
}
\examples{
data(BUPA)
# generate a linear kernel
kfun = vanilladot()

# generate a Laplacian kernel function with sigma = 1
kfun = laplacedot(sigma=1)

# generate a Gaussian kernel function with sigma estimated by sigest()
kfun = rbfdot(sigma=sigest(BUPA$X))

# set kern=kfun when fitting a kerndwd object
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
m1 = kerndwd(BUPA$X, BUPA$y, kern=kfun,
  qval=1, lambda=lambda, eps=1e-5, maxit=1e5)
}
\keyword{kernel methods}

