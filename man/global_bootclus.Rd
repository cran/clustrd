\name{global_bootclus}
\alias{global_bootclus}

\title{
Global stabiliy assessment of Joint Dimension Reduction and Clustering methods by bootstrapping.
}

\description{
Runs joint dimension and clustering algorithms repeatedly for different numbers of clusters on bootstrap replica of the original data and returns corresponding cluster assignments, and cluster agreement indices comparing pairs of partitions.
}

\usage{
global_bootclus(data, nclusrange = 3:4, ndim = NULL, 
method = c("RKM","FKM","mixedRKM","mixedFKM","clusCA","MCAk","iFCB"), 
nboot = 10, alpha = NULL, alphak = NULL, center = TRUE, 
scale = TRUE, nstart = 100, smartStart = NULL, seed = NULL)
}

\arguments{
  \item{data}{Continuous, Categorical ot Mixed data set}
  \item{nclusrange}{An integer or an integer vector with the number of clusters or a range of numbers of clusters (should be greater than one)}
  \item{ndim}{Dimensionality of the solution; if NULL it is set to nclus - 1}
  \item{method}{Specifies the method. Options are RKM for Reduced K-means, FKM for Factorial K-means, mixedRKM for Mixed Reduced K-means, mixedFKM for Mixed Factorial K-means, MCAk for MCA K-means, iFCB for Iterative Factorial Clustering of Binary variables and clusCA for Cluster Correspondence Analysis.}
  \item{nboot}{Number of bootstrap pairs of partitions}
  \item{alpha}{Adjusts for the relative importance of (mixed) RKM and FKM in the objective function; \code{alpha = 1} reduces to PCA/PCAMIX, \code{alpha = 0.5} to (mixed) reduced K-means, and \code{alpha = 0} to (mixed) factorial K-means}  
  \item{alphak}{Non-negative scalar to adjust for the relative importance of MCA (\code{alphak = 1}) and K-means (\code{alphak = 0}) in the solution (default = .5). Works only in combination with \code{method = "MCAk"}}
\item{center}{A logical value indicating whether the metric variables should be shifted to be zero centered (default = \code{TRUE)}}
\item{scale}{A logical value indicating whether the metric variables should be scaled to have unit variance before the analysis takes place (default = \code{TRUE)}}
  \item{nstart}{Number of random starts (default = 100)}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when \code{smartStart = NULL}. The default value is NULL.}
  \item{\dots}{Not used}
}

\details{
The algorithm for assessing global cluster stability is similar to that in Dolnicar and Leisch (2010) and can be summarized in three steps:

\emph{Step 1. Resampling:} Draw bootstrap samples S_i and T_i of size \emph{n} from the data and use the original data, X, as evaluation set E_i = X. Apply the clustering method of choice to S_i and T_i and obtain C^{S_i} and C^{T_i}.

\emph{Step 2. Mapping:} Assign each observation x_i to the closest centers of C^{S_i} and C^{T_i} using Euclidean distance, resulting in partitions C^{XS_i} and C^{XT_i}, where C^{XS_i} is the partition of the original data, X, predicted from clustering bootstrap sample S_i (same for T_i and C^{XT_i}).

\emph{Step 3. Evaluation:} Use the Adjusted Rand Index (ARI, Hubert & Arabie, 1985) or the Measure of Concordance (MOC, Pfitzner 2008) as measure of agreement and stability.

Inspect the distributions of ARI/MOC to assess the global reproducibility of the clustering solutions.

While nboot = 100 is recommended, smaller run numbers could give quite informative results as well, if computation times become too high.

Note that the stability of a clustering solution is assessed, but stability is not the only important validity criterion - clustering solutions obtained by very inflexible clustering methods may be stable but not valid, as discussed in Hennig (2007).
}

\value{
\item{nclusrange}{An integer or an integer vector with the number of clusters or a range of numbers of clusters}
\item{clust1}{Partitions, C^{XS_i} of the original data, X, predicted from clustering bootstrap sample S_i (see Details)}
\item{clust2}{Partitions, C^{XT_i} of the original data, X, predicted from clustering bootstrap sample T_i (see Details)}
\item{index1}{Indices of the original data rows in bootstrap sample S_i}
\item{index2}{Indices of the original data rows in bootstrap sample T_i}
\item{rand}{Adjusted Rand Index values}
\item{moc}{Measure of Concordance values}
}

\references{
Hennig, C. (2007). Cluster-wise assessment of cluster stability. \emph{Computational Statistics and Data Analysis}, \emph{52}, 258-271.

Dolnicar, S., & Leisch, F. (2010). Evaluation of structure and reproducibility of cluster solutions using the bootstrap. \emph{Marketing Letters}, \emph{21}(1), 83-101.
}

\seealso{
\code{\link{local_bootclus}}
}
\examples{
## 5 bootstrap replicates and nstart = 1 for speed in example,
## use more for real applications
data(diamond)
boot_mixedRKM = global_bootclus(diamond[,-7], nclusrange = 3:5,
method = "mixedRKM", nboot = 5, nstart = 1, seed = 1234)

boxplot(boot_mixedRKM$rand, xlab = "Number of clusters", ylab =
"adjusted Rand Index")

## 5 bootstrap replicates and nstart = 10 for speed in example,
## use more for real applications
#data(macro)
#boot_RKM = global_bootclus(macro, nclusrange = 2:5,
#method = "RKM", nboot = 5, nstart = 10, seed = 1234)

#boxplot(boot_RKM$rand, xlab = "Number of clusters", ylab =
#"adjusted Rand Index")

## 5 bootstrap replicates and nstart = 1 for speed in example,
## use more for real applications
#data(bribery)
#boot_cluCA = global_bootclus(bribery, nclusrange = 2:5, 
#method = "clusCA", nboot = 5, nstart = 1, seed = 1234)

#boxplot(boot_cluCA$rand, xlab = "Number of clusters", ylab =
#"adjusted Rand Index")
}
