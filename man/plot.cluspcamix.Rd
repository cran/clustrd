\name{plot.cluspcamix}
\alias{plot.cluspcamix}

\title{
Plotting function for \code{cluspcamix()} output.
}
\description{
Plotting function that creates a scatterplot of the objects, a correlation circle of the variables or a biplot of both objects and variables. Optionally, for metric variables, it returns a parallel coordinate plot showing cluster means and for categorical variables, a series of barplots showing the standardized residuals per attribute for each cluster.}

\usage{
\method{plot}{cluspcamix}(x, dims = c(1, 2), cludesc = FALSE, 
topstdres = 20, objlabs = FALSE, attlabs = NULL, attcatlabs = NULL, 
subplot = FALSE, what = c(TRUE,TRUE), max.overlaps = 10, \dots)
}

\arguments{
  \item{x}{Object returned by \code{cluspcamix()}}
  \item{dims}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and
second dimension vertical}
  \item{what}{Vector of two logical values specifying the contents of the plots. First entry indicates whether a scatterplot of the objects and cluster centroids is displayed and the second entry whether a correlation circle of the variables is displayed. The default is \code{c(TRUE, TRUE)} and the resultant plot is a biplot of both objects and variables}
  \item{cludesc}{A logical value indicating if a parallel coordinate plot showing cluster means is produced (default = \code{FALSE)}}
    \item{topstdres}{Number of largest standardized residuals used to describe each cluster (default = 20). Works only in combination with \code{cludesc = TRUE}}
      \item{subplot}{A logical value indicating whether a subplot with the full distribution of the standardized residuals will appear at the bottom left corner of the corresponding plots. Works only in combination with \code{cludesc = TRUE}}
      \item{objlabs}{A logical value indicating whether object labels will be plotted; if \code{TRUE} row names of the data matrix are used (default = \code{FALSE}). Warning: when \code{TRUE}, execution time of the plotting function will increase dramatically as the number of objects gets larger}
    \item{attlabs}{Vector of custom labels of continuous attributes; if not provided, default labeling is applied}
    \item{attcatlabs}{Vector of custom labels of categorical attributes (categories); if not provided, default labeling is applied}
      \item{max.overlaps}{Maximum number of text labels allowed to overlap. Defaults to 10}
  \item{\dots}{Further arguments to be transferred to \code{cluspcamix()}}
}
\value{
The function returns a ggplot2 scatterplot of the solution obtained via \code{cluspcamix()} that can be further customized using the \pkg{ggplot2} package. When \code{cludesc = TRUE}, for metric variables, the function also returns a ggplot2 parallel coordinate plot and for categorical variables, a series of ggplot2 barplots showing the largest (or all) standardized residuals per attribute for each cluster.
}
\references{
van de Velden, M., Iodice D'Enza, A., & Markos, A. (2019). Distance-based clustering of mixed data. \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, e1456.

Vichi, M., Vicari, D., & Kiers, H. A. L. (2019). Clustering and dimension reduction for mixed variables. \emph{Behaviormetrika}. doi:10.1007/s41237-018-0068-6.
}
\seealso{
\code{\link{plot.clusmca}, \link{plot.cluspca}}
}
\examples{
data(diamond)
#Mixed Reduced K-means solution with 3 clusters in 2 dimensions 
#after 10 random starts
outmixedRKM = cluspcamix(diamond, 3, 2, method = "mixedRKM", nstart = 10)
#Scatterplot (dimensions 1 and 2)
plot(outmixedRKM, cludesc = TRUE)
}