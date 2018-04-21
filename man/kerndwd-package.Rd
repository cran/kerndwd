\name{kerndwd-package}
\alias{kerndwd-package}
\docType{package}
\title{
Kernel Distance Weighted Discrimination
}
\description{
Extremely novel efficient procedures for solving linear generalized DWD and kernel generalized DWD in reproducing kernel Hilbert spaces for classification. The algorithm is based on the majorization-minimization (MM) principle to compute the entire solution path at a given fine grid of regularization parameters.\cr
}

\details{
Suppose \code{x} is predictor and \code{y} is a binary response. The package computes the entire solution path over a grid of \code{lambda} values.

The main functions of the package \code{\link{kerndwd}} include:\cr
\code{kerndwd}\cr
\code{cv.kerndwd}\cr
\code{tunedwd}\cr
\code{predict.kerndwd}\cr
\code{plot.kerndwd}\cr
\code{plot.cv.kerndwd}\cr
}

\author{Boxiang Wang and Hui Zou\cr
Maintainer: Boxiang Wang  \email{boxiang@umn.edu}}
\references{
Wang, B. and Zou, H. (2018)
``Another Look at Distance Weighted Discrimination," 
\emph{Journal of Royal Statistical Society, Series B}, \bold{80}(1), 177--198. \cr
\url{https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12244}\cr
Karatzoglou, A., Smola, A., Hornik, K., and Zeileis, A. (2004)
``kernlab -- An S4 Package for Kernel Methods in R", 
\emph{Journal of Statistical Software}, \bold{11}(9), 1--20.\cr
\url{http://www.jstatsoft.org/v11/i09/paper}\cr
Marron, J.S., Todd, M.J., Ahn, J. (2007)
``Distance-Weighted Discrimination"", 
\emph{Journal of the American Statistical Association}, \bold{102}(408), 1267--1271.\cr
\url{https://faculty.franklin.uga.edu/jyahn/sites/faculty.franklin.uga.edu.jyahn/files/DWD3.pdf}\cr
}
 
\keyword{DWD}
\keyword{kernel methods}
\keyword{classification}
\keyword{package}

