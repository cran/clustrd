\name{plot.cluspca}
\alias{plot.cluspca}

\title{
Plotting function for \code{cluspca()} output.
}
\description{
Plotting function that creates a scatterplot of the objects, a correlation circle of the variables or a biplot of both objects and variables. Optionally, it returns a parallel coordinate plot showing cluster means.}

\usage{
\method{plot}{cluspca}(x, dims = c(1, 2), cludesc = FALSE, what = c(TRUE,TRUE), attlabs, \dots)
}

\arguments{
  \item{x}{Object returned by \code{cluspca()}}
  \item{dims}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and
second dimension vertical}
  \item{what}{Vector of two logical values specifying the contents of the plots. First entry indicates whether a scatterplot of the objects and cluster centroids is displayed and the second entry whether a correlation circle of the variables is displayed. The default is \code{c(TRUE, TRUE)} and the resultant plot is a biplot of both objects and variables}
  \item{cludesc}{A logical value indicating if a parallel coordinate plot showing cluster means is produced (default = \code{FALSE)}}
    \item{attlabs}{Vector of custom attribute labels; if not provided, default labeling is applied}
  \item{\dots}{Further arguments to be transferred to \code{cluspca()}}
}
\value{
The function returns a ggplot2 scatterplot of the solution obtained via \code{cluspca()} that can be further customized using the \pkg{ggplot2} package. When \code{cludesc = TRUE} 
the function also returns a ggplot2 parallel coordinate plot.
}
\references{
De Soete, G., and Carroll, J. D. (1994). K-means clustering in a low-dimensional Euclidean space. In Diday E. et al. (Eds.), \emph{New Approaches in Classification and Data Analysis}, Heidelberg: Springer, 212-219.

Vichi, M., and Kiers, H.A.L. (2001). Factorial K-means analysis for two-way data. \emph{Computational Statistics and Data Analysis}, 37, 49-64.
}
\seealso{
\code{\link{plot.clusmca}, \link{plot.cluspcamix}}
}
\examples{
data("macro")
#Factorial K-means (3 clusters in 2 dimensions) after 100 random starts
outFKM = cluspca(macro, 3, 2, method = "FKM", rotation = "varimax")
#Scatterplot (dimensions 1 and 2) and cluster description plot
plot(outFKM, cludesc = TRUE)

data("iris", package = "datasets")
#Compromise solution between PCA and Reduced K-means
#on the iris dataset (3 clusters in 2 dimensions) after 100 random starts
outclusPCA = cluspca(iris[,-5], 3, 2, alpha = 0.3, rotation = "varimax")
table(outclusPCA$cluster,iris[,5])
#Save the ggplot2 scatterplot
map = plot(outclusPCA)$map
#Customization (adding titles)
map + ggtitle(paste("A compromise solution between RKM and FKM on the iris: 
3 clusters of sizes ", paste(outclusPCA$size, 
collapse = ", "),sep = "")) + xlab("Dimension 1") + ylab("Dimension 2") + 
theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}

