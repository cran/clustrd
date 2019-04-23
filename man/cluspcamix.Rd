\name{cluspcamix}
\alias{cluspcamix}
\alias{print.cluspcamix}
\alias{summary.cluspcamix}
\alias{fitted.cluspcamix}

\title{
Joint  dimension  reduction  and  clustering  of mixed-type data.
}

\description{
This function implements clustering and dimension reduction for mixed-type variables, i.e., categorical and metric (see, Yamamoto & Hwang, 2014; van de Velden, Iodice D'Enza, & Markos 2019; Vichi, Vicari, & Kiers, 2019). This framework includes Mixed Reduced K-means and Mixed Factorial K-means, as well as a compromise of these two methods. The methods combine Principal Component Analysis for dimension reduction with K-means for clustering.
}

\usage{
cluspcamix(data, nclus, ndim, method=c("mixedRKM","mixedFKM"), 
center = TRUE, scale = TRUE, alpha=NULL, rotation="none", 
nstart = 100, smartStart=NULL, seed=NULL, binary = FALSE)

\method{print}{cluspcamix}(x, \dots)

\method{summary}{cluspcamix}(object, \dots)

\method{fitted}{cluspcamix}(object, mth = c("centers", "classes"), \dots)
}

\arguments{
  \item{data}{Dataset with categorical and metric variables}
  \item{nclus}{Number of clusters (nclus = 1 returns the PCAMIX solution; see Details)}
  \item{ndim}{Dimensionality of the solution}
    \item{method}{Specifies the method. Options are mixedRKM for mixed reduced K-means and mixedFKM for mixed factorial K-means (default = \code{"mixedRKM"})}
 \item{center}{A logical value indicating whether the variables should be shifted to be zero centered (default = \code{TRUE)}}
  \item{scale}{A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place (default = \code{TRUE)}}
  \item{alpha}{Adjusts for the relative importance of Mixed RKM and Mixed FKM in the objective function; \code{alpha} = 0.5 leads to mixed reduced K-means, \code{alpha} = 0 to mixed factorial K-means, and \code{alpha} = 1 reduces to the tandem approach (PCAMIX followed by K-means)}
    \item{rotation}{Specifies the method used to rotate the factors. Options are \code{none} for no rotation, \code{varimax} for varimax rotation with Kaiser normalization and \code{promax} for promax rotation (default = \code{"none"})}
  \item{nstart}{Number of random starts (default = 100)}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when \code{smartStart = NULL}. The default value is NULL.}
  \item{binary}{If \code{TRUE} then all categorical variables are 0-1 (dummy) variables.}
  \item{x}{For the \code{print} method, a class of \code{cluspcamix}}
  \item{object}{For the \code{summary} method, a class of \code{cluspcamix}}
  \item{mth}{For the \code{fitted} method, a character string that specifies the type of fitted value to return: \code{"centers"} for the observations center vector, or \code{"class"} for the observations cluster membership value}
  \item{\dots}{Not used}
}

\details{
For the K-means part, the algorithm of Hartigan-Wong is used by default.

The hidden \code{print} and \code{summary} methods print out some key components of an object of class \code{cluspcamix}. 

The hidden \code{fitted} method returns cluster fitted values. If method is \code{"classes"}, this is a vector of cluster membership (the cluster component of the "cluspcamix" object). If method is \code{"centers"}, this is a matrix where each row is the cluster center for the observation. The rownames of the matrix are the cluster membership values.

When \code{nclus} = 1 the function returns the solution of PCAMIX (Kiers, 1991) and \code{plot(object)} shows the corresponding biplot.
}

\value{
\item{obscoord}{Object scores}
\item{attcoord}{Variable scores}
\item{centroid}{Cluster centroids}
\item{cluster}{Cluster membership}
\item{criterion}{Optimal value of the objective criterion}
\item{size}{The number of objects in each cluster}
\item{scale}{A copy of \code{scale} in the return object}
\item{center}{A copy of \code{center} in the return object}
\item{nstart}{A copy of \code{nstart} in the return object}
\item{odata}{A copy of \code{data} in the return object}
}
\references{
van de Velden, M., Iodice D'Enza, A., & Markos, A. (2019).
Distance-based clustering of mixed data. \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, e1456.

Vichi, M., Vicari, D., & Kiers, H.A.L. (2019). Clustering and
dimension reduction for mixed variables. \emph{Behaviormetrika}. doi:10.1007/s41237-018-0068-6.

Yamamoto, M., & Hwang, H. (2014). A general formulation of
cluster analysis with dimension reduction and subspace 
separation. \emph{Behaviormetrika}, \emph{41}, 115-129.
}

\seealso{
\code{\link{cluspca}}, \code{\link{clusmca}}, \code{\link{tuneclus}}
}
\examples{
data(diamond)
#Mixed Reduced K-means solution with 3 clusters in 2 dimensions 
#after 10 random starts
outmixedRKM = cluspcamix(diamond, 3, 2, method = "mixedRKM", nstart = 10, seed = 1234)
outmixedRKM 
#Scatterplot (dimensions 1 and 2)
plot(outmixedRKM)

#Tandem analysis: PCAMIX followed by K-means solution 
#with 3 clusters in 2 dimensions after 10 random starts 
outTandem = cluspcamix(diamond, 3, 2, alpha = 1, nstart = 10, seed = 1234)
outTandem
#Scatterplot (dimensions 1 and 2)
plot(outTandem)

#nclus = 1 just gives the PCAMIX solution
#outPCAMIX = cluspcamix(diamond, 1, 2)
#outPCAMIX
#Biplot (dimensions 1 and 2) 
#plot(outPCAMIX)
}
