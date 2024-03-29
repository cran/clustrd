\name{plot.clusmca}
\alias{plot.clusmca}
\title{
Plotting function for \code{clusmca()} output.
}
\description{
Plotting function that creates a scatterplot of the object scores and/or the attribute scores and the cluster centroids. Optionally, the function returns a series of barplots showing the standardized residuals per attribute for each cluster. 
}
\usage{
\method{plot}{clusmca}(x, dims = c(1,2), what = c(TRUE,TRUE),
cludesc = FALSE, topstdres = 20, objlabs = FALSE, attlabs = NULL, 
subplot = FALSE, max.overlaps=10, \dots)
}

\arguments{
  \item{x}{Object returned by \code{clusmca()}}
  \item{dims}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and
second dimension vertical}
  \item{what}{Vector of two logical values specifying the contents of the plots. First entry indicates whether a scatterplot of the objects is displayed in principal coordinates. Second entry indicates whether a scatterplot of the attribute categories is displayed in principal coordinates. Cluster centroids are always displayed. The default is \code{c(TRUE, TRUE)} and the resultant plot is a biplot of both objects and attribute categories with gamma-based scaling (see van de Velden et al., 2017)}
  \item{cludesc}{A logical value indicating whether a series of barplots is produced showing the largest (in absolute value) standardized residuals per attribute for each cluster (default = \code{FALSE})}
  \item{topstdres}{Number of largest standardized residuals used to describe each cluster (default = 20). Works only in combination with \code{cludesc = TRUE}}
  \item{objlabs}{A logical value indicating whether object labels will be plotted; if \code{TRUE} row names of the data matrix are used (default = \code{FALSE}). Warning: when \code{TRUE}, execution time of the plotting function will increase dramatically as the number of objects gets larger}
  \item{attlabs}{Vector of custom attribute labels; if not provided, default labeling is applied}
  \item{subplot}{A logical value indicating whether a subplot with the full distribution of the standardized residuals will appear at the bottom left corner of the corresponding plots. Works only in combination with \code{cludesc = TRUE}}
  \item{max.overlaps}{Maximum number of text labels allowed to overlap. Defaults to 10}
  \item{\dots}{Further arguments to be transferred to \code{clusmca()}}
}
\value{
The function returns a ggplot2 scatterplot of the solution obtained via \code{clusmca()} that can be further customized using the \pkg{ggplot2} package. When \code{cludesc = TRUE} the function also returns a series of ggplot2 barplots showing the largest (or all) standardized residuals per attribute for each cluster.
}

\references{
Hwang, H., Dillon, W. R., and Takane, Y. (2006). An extension of multiple correspondence analysis for identifying heterogenous subgroups of respondents. \emph{Psychometrika}, 71, 161-171.

Iodice D'Enza, A., and Palumbo, F. (2013). Iterative factor clustering of binary data. \emph{Computational Statistics}, \emph{28}(2), 789-807.

van de Velden M., Iodice D'Enza, A., and Palumbo, F. (2017).  Cluster correspondence analysis. \emph{Psychometrika}, \emph{82}(1), 158-185.
}
\seealso{
\code{\link{plot.cluspca}, \link{plot.cluspcamix}}
}
\examples{
data("mybond")
#Cluster Correspondence Analysis with 3 clusters in 2 dimensions after 10 random starts
outclusCA = clusmca(mybond[,8:10], 3, 2, nstart = 100, seed = 234)
#Save the ggplot2 scatterplot  
map = plot(outclusCA, max.overlaps = 40)$map
#Customization (adding titles)
map + ggtitle(paste("Cluster CA plot of the James bond categorical data: 3 clusters of sizes ", 
                    paste(outclusCA$size, collapse = ", "),sep = "")) + 
    xlab("Dim. 1") + ylab("Dim. 2") + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
    
data("mybond")
#i-FCB with 3 clusters in 2 dimensions after 10 random starts
outclusCA = clusmca(mybond[,8:10], 3, 2, method = "iFCB", nstart= 10)
#Scatterlot with the observations only (dimensions 1 and 2) 
#and cluster description plots showing the 20 largest std. residuals 
#(with the full distribution showing in subplots)
plot(outclusCA, dim = c(1,2), what = c(TRUE, FALSE), cludesc = TRUE, 
subplot = TRUE)
}

