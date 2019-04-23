\name{tuneclus}
\alias{tuneclus}
\alias{print.tuneclus}
\alias{summary.tuneclus}
\alias{fitted.tuneclus}

\title{
Cluster quality assessment for a range of clusters and dimensions.
}

\description{
This function facilitates the selection of the appropriate number of clusters and dimensions for joint dimension reduction and clustering methods. 
}

\usage{
tuneclus(data, nclusrange = 3:4, ndimrange = 2:3, 
method = c("RKM","FKM","mixedRKM","mixedFKM","clusCA","iFCB","MCAk"), 
criterion = "asw", dst = "full", alpha = NULL, alphak = NULL, 
center = TRUE, scale = TRUE, rotation = "none", nstart = 100, 
smartStart = NULL, seed = NULL)

\method{print}{tuneclus}(x, \dots)

\method{summary}{tuneclus}(object, \dots)

\method{fitted}{tuneclus}(object, mth = c("centers", "classes"), \dots)

}

\arguments{
  \item{data}{Continuous, Categorical ot Mixed data set}
  \item{nclusrange}{An integer vector with the range of numbers of clusters which are to be compared by the cluster validity criteria. Note: the number of clusters should be greater than one}
  \item{ndimrange}{An integer vector with the range of dimensions which are to be compared by the cluster validity criteria}
  \item{method}{Specifies the method. Options are \code{RKM} for reduced K-means, \code{FKM} for factorial K-means, \code{mixedRKM} for mixed reduced K-means, \code{mixedFKM} for mixed factorial K-means, \code{MCAk} for MCA K-means, \code{iFCB} for Iterative Factorial Clustering of Binary variables and \code{clusCA} for Cluster Correspondence Analysis}
  \item{criterion}{One of \code{asw}, \code{ch} or \code{crit}. Determines whether average silhouette width, Calinski-Harabasz index or objective value of the selected method is used (default = \code{"asw")}}
  \item{dst}{Specifies the data used to compute the distances between objects. Options are \code{full} for the original data (after possible scaling) and \code{low} for the object scores in the low-dimensional space (default = \code{"full")}}
    \item{alpha}{Adjusts for the relative importance of (mixed) RKM and FKM in the objective function; \code{alpha = 1} reduces to PCA, \code{alpha = 0.5} to (mixed) reduced K-means, and \code{alpha = 0} to (mixed) factorial K-means}
    \item{alphak}{Non-negative scalar to adjust for the relative importance of MCA (\code{alphak = 1}) and K-means (\code{alphak = 0}) in the solution (default = .5). Works only in combination with \code{method = "MCAk"}}
  \item{center}{A logical value indicating whether the variables should be shifted to be zero centered (default = \code{TRUE)}}
  \item{scale}{A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place (default = \code{TRUE)}}
  \item{rotation}{Specifies the method used to rotate the factors. Options are none for no rotation, varimax for varimax rotation with Kaiser normalization and promax for promax rotation (default = \code{"none")}}
  \item{nstart}{Number of starts (default = 100)}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when smartStart = NULL. The default value is NULL.}
  \item{x}{For the \code{print} method, a class of \code{clusmca}}
  \item{object}{For the \code{summary} method, a class of \code{clusmca}}
  \item{mth}{For the \code{fitted} method, a character string that specifies the type of fitted value to return: \code{"centers"} for the observations center vector, or \code{"class"} for the observations cluster membership value}
  \item{\dots}{Not used}
}

\details{
For the K-means part, the algorithm of Hartigan-Wong is used by default.

The hidden \code{print} and \code{summary} methods print out some key components of an object of class \code{tuneclus}. 

The hidden \code{fitted} method returns cluster fitted values. If method is \code{"classes"}, this is a vector of cluster membership (the cluster component of the "tuneclus" object). If method is \code{"centers"}, this is a matrix where each row is the cluster center for the observation. The rownames of the matrix are the cluster membership values.}

\value{
\item{clusobjbest}{The output of the optimal run of \code{cluspca()} or \code{clusmca()}}
\item{nclusbest}{The optimal number of clusters}
\item{ndimbest}{The optimal number of dimensions}
\item{critbest}{The optimal criterion value for \code{nclusbest} clusters and \code{ndimbest} dimensions}
\item{critgrid}{Matrix of size \code{nclusrange x ndimrange} with the criterion values for the specified ranges of clusters and dimensions (values are calculated only when the number of clusters is greater than the number of dimensions; otherwise values in the grid are left blank)}}

\references{Calinski, R.B., and Harabasz, J., (1974). A dendrite method for cluster analysis. \emph{Communications in Statistics}, 3, 1-27.

Kaufman, L., and Rousseeuw, P.J., (1990). \emph{Finding Groups in Data: An Introduction to Cluster Analysis}. Wiley, New York.}

\seealso{
\code{\link{cluspca}}, \code{\link{clusmca}}
}

\examples{
# Reduced K-means for a range of clusters and dimensions
data(macro)
# Cluster quality assessment based on the average silhouette width in the low dimensional space
# nstart = 1 for speed in example
# use more for real applications
bestRKM = tuneclus(macro, 3:4, 2:3, method = "RKM", 
criterion = "asw", dst = "low", nstart = 1, seed = 1234)
bestRKM
#plot(bestRKM)

# Cluster Correspondence Analysis for a range of clusters and dimensions
data(bribery)
# Cluster quality assessment based on the Callinski-Harabasz index in the full dimensional space
bestclusCA = tuneclus(bribery, 4:5, 3:4, method = "clusCA",
criterion = "ch", nstart = 20, seed = 1234)
bestclusCA
#plot(bestclusCA, cludesc = TRUE)

# Mixed reduced K-means for a range of clusters and dimensions
data(diamond)
# Cluster quality assessment based on the average silhouette width in the low dimensional space
# nstart = 5 for speed in example
# use more for real applications
bestmixedRKM = tuneclus(diamond[,-7], 3:4, 2:3, 
method = "mixedRKM", criterion = "asw", dst = "low", 
nstart = 5, seed = 1234)
bestmixedRKM
#plot(bestmixedRKM)
}
