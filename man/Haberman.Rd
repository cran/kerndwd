\name{Haberman}
\alias{Haberman}

\title{Haberman's survival data set from Alon et al. (1999)}
\usage{
data(Haberman)
}

\description{
Haberman's survival data (3 predictors for 306 samples) from a study on the survival of patients who had taken breast cancer surgery. 
}

\details{This data set contains 305 observations and 3 predictors representing the survival status of 305 patients. The three predictors are age of patient when the surgery was conducted, years of operations, and the number of positive axillary nodes detected.  Among 305 subjects, 224 survived 5 years or longer and 81 died within 5 years.  
}

\value{
A list with the following elements:
  \item{X}{A numerical matrix for predictors: 305 rows and 3 columns; each row corresponds to a patient.}
  \item{y}{A numeric vector of length 305 representing the survival status (1 for died within 5 years; -1 for survived 5 years or longer).}
}

\source{The data set is available for download from UCI machine learning repository. \cr
\url{https://archive.ics.uci.edu/ml/datasets/Haberman's+Survival}}

\references{
Haberman, S. J. (1976).
``Generalized residuals for log-linear models'',
\emph{Proceedings of the 9th International Biometrics Conference}, Boston, 104--122.\cr
\cr
Lichman, M. (2013). ``UCI Machine Learning Repository", Irvine, CA: University of California, School of Information and Computer Science.\cr
\url{http://archive.ics.uci.edu/ml}
}

\examples{
# load kerndwd library
library(kerndwd)

# load data set
data(Haberman)

# the number of samples predictors
dim(Haberman$X)

# the number of samples for each class
sum(Haberman$y == -1) # survived
sum(Haberman$y == 1)  # died
}

\keyword{data set}
