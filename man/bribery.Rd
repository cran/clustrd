\name{bribery}
\alias{bribery}
\docType{data}
\title{Bribery cases in Russia}
\description{
The data set refers to a collection of 55 articles on Bribery cases from central Russian newspapers 1999-2000 (Mirkin, 2005). The variables reflect the following fivefold structure of bribery situations: two interacting sides - the office and the client, their interaction, the corrupt service rendered, and the environment in which it all occurs. These structural aspects can be characterized by eleven variables that have been manually recovered from the newspaper articles.}

\usage{data("bribery")}
\format{
   A data frame with 55 observations on 11 categorical variables).
  \describe{
    \item{\code{Of}}{Type of Office}
    \item{\code{Cl}}{Level of Client}
    \item{\code{Serv}}{Type of service: obstruction of justice, favours, cover-up, change of category, extortion of money for rendering free services}
    \item{\code{Occ}}{Frequency of occurrence}
    \item{\code{Init}}{Who initiated the bribery act}
    \item{\code{Brib}}{Bribe Level in $}
    \item{\code{Typ}}{Type of corruption}
    \item{\code{Net}}{Corruption network}
    \item{\code{Con}}{Condition of corruption}
    \item{\code{Bran}}{Branch at which the corrupt service occurred}
    \item{\code{Pun}}{Punishment}
  }
}
\references{
Mirkin, B. (2005). \emph{Clustering for data mining: a data recovery approach}. Chapman and Hall/CRC.
}

\examples{
data(bribery)
}
\keyword{datasets}