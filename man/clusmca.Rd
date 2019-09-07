\name{clusmca}
\alias{clusmca}
\alias{print.clusmca}
\alias{summary.clusmca}
\alias{fitted.clusmca}

\title{
Joint  dimension  reduction  and  clustering  of  categorical  data.
}

\description{
This function implements MCA K-means (Hwang, Dillon and Takane, 2006), i-FCB (Iodice D' Enza and Palumbo, 2013) and Cluster Correspondence Analysis (van de Velden, Iodice D' Enza and Palumbo, 2017). The methods combine variants of Correspondence Analysis for dimension reduction with K-means for clustering.
}

\usage{
clusmca(data, nclus, ndim, method=c("clusCA","iFCB","MCAk"),
alphak = .5, nstart = 100, smartStart = NULL, gamma = TRUE, 
inboot = FALSE, seed = NULL)

\method{print}{clusmca}(x, \dots)

\method{summary}{clusmca}(object, \dots)

\method{fitted}{clusmca}(object, mth = c("centers", "classes"), \dots)
}

\arguments{
  \item{data}{Dataset with categorical variables}
  \item{nclus}{Number of clusters (nclus = 1 returns the MCA solution; see Details)}
  \item{ndim}{Dimensionality of the solution}
  \item{method}{Specifies the method. Options are MCAk for MCA K-means, iFCB for Iterative Factorial Clustering of Binary variables and clusCA for Cluster Correspondence Analysis (default = \code{"clusCA")}}
  \item{alphak}{Non-negative scalar to adjust for the relative importance of MCA (\code{alphak = 1}) and K-means (\code{alphak = 0}) in the solution (default = .5). Works only in combination with \code{method = "MCAk"}}
  \item{nstart}{Number of random starts (default = 100)}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{gamma}{Scaling parameter that leads to similar spread in the object and variable scores (default = \code{TRUE})}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when \code{smartStart = NULL}. The default value is NULL.}
  \item{inboot}{Used internally in the bootstrap functions to perform bootstrapping on the indicator matrix.}
  \item{x}{For the \code{print} method, a class of \code{clusmca}}
  \item{object}{For the \code{summary} method, a class of \code{clusmca}}
  \item{mth}{For the \code{fitted} method, a character string that specifies the type of fitted value to return: \code{"centers"} for the observations center vector, or \code{"class"} for the observations cluster membership value}
  \item{\dots}{Not used}
}

\details{
For the K-means part, the algorithm of Hartigan-Wong is used by default.

The hidden \code{print} and \code{summary} methods print out some key components of an object of class \code{clusmca}. 

The hidden \code{fitted} method returns cluster fitted values. If method is \code{"classes"}, this is a vector of cluster membership (the cluster component of the "clusmca" object). If method is \code{"centers"}, this is a matrix where each row is the cluster center for the observation. The rownames of the matrix are the cluster membership values.

When \code{nclus} = 1 the function returns the MCA solution with objects in principal and variables in standard coordinates. \code{plot(object)} shows the corresponding asymmetric biplot.
}

\value{
\item{obscoord}{Object scores}
\item{attcoord}{Variable scores}
\item{centroid}{Cluster centroids}
\item{cluster}{Cluster membership}
\item{criterion}{Optimal value of the objective criterion}
\item{size}{The number of objects in each cluster}
\item{nstart}{A copy of \code{nstart} in the return object}
\item{odata}{A copy of \code{data} in the return object}
}
\references{
Hwang, H., Dillon, W. R., and Takane, Y. (2006). An extension of multiple correspondence analysis for identifying heterogenous subgroups of respondents. \emph{Psychometrika}, 71, 161-171.

Iodice D'Enza, A., and Palumbo, F. (2013). Iterative factor clustering of binary data. \emph{Computational Statistics}, \emph{28}(2), 789-807.

van de Velden M., Iodice D' Enza, A., and Palumbo, F. (2017).  Cluster correspondence analysis. \emph{Psychometrika}, \emph{82}(1), 158-185.
}

\seealso{
\code{\link{cluspca}}, \code{\link{cluspcamix}}, \code{\link{tuneclus}}
}
\examples{
data(cmc)
# Preprocessing: values of wife's age and number of children were categorized 
# into three groups based on quartiles
cmc$W_AGE = ordered(cut(cmc$W_AGE, c(16,26,39,49), include.lowest = TRUE))
levels(cmc$W_AGE) = c("16-26","27-39","40-49") 
cmc$NCHILD = ordered(cut(cmc$NCHILD, c(0,1,4,17), right = FALSE))
levels(cmc$NCHILD) = c("0","1-4","5 and above")

#Cluster Correspondence Analysis solution with 3 clusters in 2 dimensions 
#after 10 random starts
outclusCA = clusmca(cmc, 3, 2, method = "clusCA", nstart = 10, seed = 1234)
outclusCA
#Scatterplot (dimensions 1 and 2)
plot(outclusCA)

#MCA K-means solution with 3 clusters in 2 dimensions after 10 random starts
outMCAk = clusmca(cmc, 3, 2, method = "MCAk", nstart = 10, seed = 1234)
outMCAk
#Scatterplot (dimensions 1 and 2)
plot(outMCAk)

#nclus = 1 just gives the MCA solution
#outMCA = clusmca(cmc, 1, 2)
#outMCA
#Scatterplot (dimensions 1 and 2) 
#asymmetric biplot with scaling gamma = TRUE
#plot(outMCA)
}
