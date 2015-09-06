\name{BUPA}
\alias{BUPA}

\title{BUPA's liver disorders data}
\usage{
data(BUPA)
}

\description{
BUPA's liver disorders data: 345 male individuals' blood test result and liver disorder status. 
}

\details{This data set consists of 345 observations and 6 predictors representing the blood test result liver disorder status of 345 patients. The three predictors are mean corpuscular volume (MCV), alkaline phosphotase (ALKPHOS), alamine aminotransferase (SGPT), aspartate aminotransferase (SGOT), gamma-glutamyl transpeptidase (GAMMAGT), and the number of alcoholic beverage drinks per day (DRINKS). 
}

\value{
A list with the following elements:
  \item{X}{A numerical matrix for predictors: 345 rows and 6 columns; each row corresponds to a patient.}
  \item{y}{A numeric vector of length 305 representing the liver disorder status.}
}

\source{The data set is available for download from UCI machine learning repository. \cr
\url{https://archive.ics.uci.edu/ml/datasets/Liver+Disorders}}

\examples{
# load data set
data(BUPA)

# the number of samples predictors
dim(BUPA$X)

# the number of samples for each class
sum(BUPA$y == -1) 
sum(BUPA$y == 1)
}

\keyword{data set}
