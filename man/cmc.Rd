\name{cmc}
\alias{cmc}
\docType{data}
\title{Contraceptive Choice in Indonesia
}
\description{
Data of married women in Indonesia who were not pregnant (or did not know they were pregnant) at the time of the survey. The dataset contains demographic and socio-economic characteristics of the women along with their preferred method of contraception (no use, long-term methods, short-term methods).
}
\usage{data(cmc)}
\format{
  A data frame containing 1,437 observations on the following 10 variables.
  \describe{
    \item{\code{W_AGE}}{wife's age in years.}
    \item{\code{W_EDU}}{ordered factor indicating wife's education, with levels \code{"low"}, \code{"2"}, \code{"3"} and \code{"high".}}
    \item{\code{H_EDU}}{ordered factor indicating wife's education, with levels \code{"low"}, \code{"2"}, \code{"3"} and \code{"high".}}
    \item{\code{NCHILD}}{number of children.}
    \item{\code{W_REL}}{factor indicating wife's religion, with levels \code{"non-Islam"} and \code{"Islam"}.}
     \item{\code{W_WORK}}{factor indicating if the wife is working.}
     \item{\code{H_OCC}}{factor indicating husband's occupation, with levels \code{"1"}, \code{"2"}, \code{"3"} and \code{"4"}. The labels are not known.}
      \item{\code{SOL}}{ordered factor indicating the standard of living index with levels \code{"low"}, \code{"2"}, \code{"3"} and \code{"high"}.}
       \item{\code{MEDEXP}}{factor indicating media exposure, with levels \code{"good"} and \code{"not good"}.}
        \item{\code{CM}}{factor indicating the contraceptive method used, with levels \code{"no-use"}, \code{"long-term"} and \code{"short-term"}.}
  }
}
\source{
This dataset is part of the 1987 National Indonesia Contraceptive Prevalence Survey and was created by Tjen-Sien Lim. It has been taken from the UCI Machine Learning Repository at \url{http://archive.ics.uci.edu/ml/}.
}
\references{
Lim, T.-S., Loh, W.-Y. & Shih, Y.-S. (1999). A Comparison of Prediction Accuracy, Complexity, and Training Time of Thirty-three Old and New Classification Algorithms. \emph{Machine Learning}, \emph{40}(3), 203-228.
}
\examples{
data(cmc)
}
\keyword{datasets}
