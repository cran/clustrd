\name{mybond}
\alias{mybond}
\docType{data}
\title{James Bond films}
\description{
The data set refers to 26 James Bond films produced up to 2021, based on 10 film characteristics: 7 continuous (year of release, production budget, box office gross in the USA and worldwide, running time, IMDB average rating, Rotten Tomatoes rating) and 3 categorical (Bond actor, native country of the actor playing the villain, native country of the actor playing the Bond girl). All figures in USD are adjusted for inflation. Most of the data was compiled from the Wikipedia page: \url{https://en.wikipedia.org/wiki/List_of_James_Bond_films}.}

\usage{data("mybond")}
\format{
   A data frame with 26 observations on the following 10 variables.
  \describe{
    \item{\code{year}}{Year of release}
    \item{\code{budget}}{Official production budget (in million USD)}
    \item{\code{grossusa}}{Box office gross in the USA (in million USD)}
    \item{\code{grosswrld}}{Box office gross worldwide (in million USD)}
    \item{\code{rtime}}{Running time in minutes}
    \item{\code{IMDB}}{IMDB rating}
    \item{\code{rottentomatoes}}{Rotten Tomatoes rating}
    \item{\code{actor}}{Bond actor}
    \item{\code{villaincnt}}{Native country of the actor playing the villain}
    \item{\code{bondgirlcnt}}{Native country of the actor playing the Bond girl}
  }
}


\examples{
data(mybond)
}
\keyword{datasets}