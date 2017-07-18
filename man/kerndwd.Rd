\name{kerndwd}
\alias{kerndwd}
\title{solve Linear DWD and Kernel DWD}
\description{Fit the linear distance weighted discrimination (DWD) model and the DWD on Reproducing kernel Hilbert space. The solution path is computed at a grid of values of tuning parameter \code{lambda}.}
\usage{
kerndwd(x, y, kern, lambda, qval=1, wt, eps=1e-05, maxit=1e+05)
}
\arguments{
	\item{x}{A numerical matrix with \eqn{N} rows and \eqn{p} columns for predictors.}
	\item{y}{A vector of length \eqn{N} for binary responses. The element of \code{y} is either -1 or 1.}
	\item{kern}{A kernel function; see \code{\link{dots}}.}
  \item{lambda}{A user supplied \code{lambda} sequence.}
  \item{qval}{The exponent index of the generalized DWD. Default value is 1.}
  \item{wt}{A vector of length \eqn{n} for weight factors. When \code{wt} is missing or \code{wt=NULL}, an unweighted DWD is fitted. } 
	\item{eps}{The algorithm stops when (i.e. \eqn{\sum_j(\beta_j^{new}-\beta_j^{old})^2}{sum(j)(beta_new[j]-beta_old[j])^2} is less than \code{eps}, where \eqn{j=0,\ldots, p}. Default value is \code{1e-5}.}
	\item{maxit}{The maximum of iterations allowed. Default is 1e5.}
}

\details{Suppose that the generalized DWD loss is \eqn{V_q(u)=1-u}{V_q(u) = 1 - u} if \eqn{u \le q/(q+1)}{u <= q/(q+1)} and \eqn{\frac{1}{u^q}\frac{q^q}{(q+1)^{(q+1)}}}{(1/u)^q * q^q/(q+1)^{(q+1)}} if \eqn{u > q/(q+1)}{u > q/(q+1)}. The value of \eqn{\lambda}, i.e., \code{lambda}, is user-specified. 

In the linear case (\code{kern} is the inner product and N > p), the \code{\link{kerndwd}} fits a linear DWD by minimizing the L2 penalized DWD loss function,
\deqn{\frac{1}{N}\sum_{i=1}^n V_q(y_i(\beta_0 + X_i'\beta)) + \lambda \beta' \beta.}{(1/N) * sum_i [V_q(y_i(\beta_0 + X_i'\beta))] + \lambda \beta' \beta.} 

In the kernel case, the \code{\link{kerndwd}} fits a kernel DWD by minimizing
\deqn{\frac{1}{N}\sum_{i=1}^n V_q(y_i(\beta_0 + K_i' \alpha)) + \lambda \alpha' K \alpha,}{(1/N) * sum_i [V_q(y_i(\beta_0 + K_i' \alpha))] + \lambda \alpha' K \alpha,}
where \eqn{K}{K} is the kernel matrix and \eqn{K_i}{K_i} is the ith row. 

The weighted linear DWD and the weighted kernel DWD are formulated as follows,
\deqn{\frac{1}{N}\sum_{i=1}^n w_i \cdot V_q(y_i(\beta_0 + X_i'\beta)) + \lambda \beta' \beta,}{(1/N) * sum_i [w_i * V_q(y_i(\beta_0 + X_i'\beta))] + \lambda \beta' \beta,}
\deqn{\frac{1}{N}\sum_{i=1}^n w_i \cdot V_q(y_i(\beta_0 + K_i' \alpha)) + \lambda \alpha' K \alpha,}{(1/N) * sum_i [w_i * V_q(y_i(\beta_0 + K_i' \alpha))] + \lambda \alpha' K \alpha,}
where \eqn{w_i} is the ith element of \code{wt}. The choice of weight factors can be seen in the reference below.
}

\value{
An object with S3 class \code{\link{kerndwd}}.
  \item{alpha}{A matrix of DWD coefficients at each \code{lambda} value. The dimension is \code{(p+1)*length(lambda)} in the linear case and \code{(N+1)*length(lambda)} in the kernel case.}
  \item{lambda}{The \code{lambda} sequence.}
  \item{npass}{Total number of MM iterations for all lambda values. }
  \item{jerr}{Warnings and errors; 0 if none.}
  \item{info}{A list including parameters of the loss function, \code{eps}, \code{maxit}, \code{kern}, and \code{wt} if a weight vector was used.}
  \item{call}{The call that produced this object.}
}

\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang@umn.edu}}
\references{
Wang, B. and Zou, H. (2017)
``Another Look at Distance Weighted Discrimination," \emph{Journal of Royal Statistical Society, Series B}, to appear. \cr
\url{http://arxiv.org/abs/1508.05913}\cr

Karatzoglou, A., Smola, A., Hornik, K., and Zeileis, A. (2004)
``kernlab -- An S4 Package for Kernel Methods in R", 
\emph{Journal of Statistical Software}, \bold{11}(9), 1--20.\cr
\url{http://www.jstatsoft.org/v11/i09/paper}\cr
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software}, \bold{33}(1), 1--22.\cr
\url{http://www.jstatsoft.org/v33/i01/paper}\cr
Marron, J.S., Todd, M.J., and Ahn, J. (2007)
``Distance-Weighted Discrimination"", 
\emph{Journal of the American Statistical Association}, \bold{102}(408), 1267--1271.\cr
\url{https://faculty.franklin.uga.edu/jyahn/sites/faculty.franklin.uga.edu.jyahn/files/DWD3.pdf}\cr
Qiao, X., Zhang, H., Liu, Y., Todd, M., Marron, J.S. (2010)
``Weighted distance weighted discrimination and its asymptotic properties", 
\emph{Journal of the American Statistical Association}, \bold{105}(489), 401--414.\cr
\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996856/}\cr
}

\seealso{\code{\link{predict.kerndwd}}, \code{\link{plot.kerndwd}}, and \code{\link{cv.kerndwd}}.}
\examples{
data(BUPA)
# standardize the predictors
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)

# a grid of tuning parameters
lambda = 10^(seq(3, -3, length.out=10))

# fit a linear DWD
kern = vanilladot()
DWD_linear = kerndwd(BUPA$X, BUPA$y, kern,
  qval=1, lambda=lambda, eps=1e-5, maxit=1e5)

# fit a DWD using Gaussian kernel
kern = rbfdot(sigma=1)
DWD_Gaussian = kerndwd(BUPA$X, BUPA$y, kern,
  qval=1, lambda=lambda, eps=1e-5, maxit=1e5)

# fit a weighted kernel DWD
kern = rbfdot(sigma=1)
weights = c(1, 2)[factor(BUPA$y)]
DWD_wtGaussian = kerndwd(BUPA$X, BUPA$y, kern,
  qval=1, lambda=lambda, wt = weights, eps=1e-5, maxit=1e5)
}
\keyword{DWD}
\keyword{kernel methods}
\keyword{classification}
