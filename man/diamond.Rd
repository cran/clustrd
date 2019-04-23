\name{diamond}
\alias{diamond}
\docType{data}
\title{Diamond Stone Pricing}
  \description{Data on 308 diamond stones sold in Singapore. The main attributes are diamond weight, colour, clarity, certification body and price in Singapore $. The weight of a diamond stone is indicated in terms of carat units. Since stones may be divided into 3 clusters due to their size, namely small (less than 0.5 carats), medium (0.5 to less than 1 carat) and large (1 carat and over), following Chu (2001), three binary variables have been built representing the three caratage ranges, and three quantitative variables (denoted Small, Medium, Large) have been derived by multiplying such binary variables by carats. So, the "Small" variable has nonzero values (i.e., the carat values) only for the smallest diamonds (less than 0.5 carats), and likewise for the other two variables. Thus, these variables are weighted binary variables. The colour of a diamond is graded from D (completely colourless), E, F, G, ..., to I (almost colorless). Clarity refers to the diamond's internal and external imperfections. Clarity is graded on a scale from IF (internally flawless), to very very slightly imperfect (VVS1 or VVS2), and very slightly imperfect, VS1 or VS2. Three certification bodies were used: New York based Gemmological Institute of America (GIA), Antwerp based International Gemmological Institute (IGI) and Hoge Raad Voor Diamant (HRD).  
}
\usage{data(diamond)}
\format{
  A data frame with 308 observations on the following 7 variables.
  \describe{
    \item{\code{Small}}{weighted binary variable with nonzero values (i.e., the carat values) for diamonds with less than 0.5 carats.}
    \item{\code{Medium}}{weighted binary variable with nonzero values (i.e., the carat values) for diamonds from 0.5 to less than 1 carat.}
    \item{\code{Large}}{weighted binary variable with nonzero values (i.e., the carat values) for diamonds from 1 carat and over.}
    \item{\code{Colour}}{the color of the diamond with a factor with levels (D, E, F, G, H, I).}
    \item{\code{Clarity}}{the clarity of the diamond with a factor with levels (IF, VVS1, VVS2, VS1, VS2).}
    \item{\code{Certification}}{the certification body with a factor with levels (GIA, IGI, HRD).}
    \item{\code{Price}}{the price of a diamond in Singapore $.}
  }
}

\references{Chu, S. (2001). Pricing the C's of Diamond Stones, \emph{Journal of Statistics Education}, \emph{9}(2).}
\keyword{datasets}
