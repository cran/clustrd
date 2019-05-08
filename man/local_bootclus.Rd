\name{local_bootclus}
\alias{local_bootclus}

\title{
Cluster-wise stability assessment of Joint Dimension Reduction and Clustering methods by bootstrapping.
}

\description{
Assessment of the cluster-wise stability of a joint dimension and clustering method. The data is resampled using bootstrapping and the Jaccard similarities of the original clusters to the most similar clusters in the resampled data are computed. The mean over these similarities is used as an index of the stability of a cluster. The method is similar to the one described in Hennig (2007). 
}

\usage{
local_bootclus(data, nclus, ndim = NULL, 
method = c("RKM","FKM","mixedRKM","mixedFKM","clusCA","MCAk","iFCB"), 
scale = TRUE, center= TRUE, alpha = NULL, nstart=100, 
nboot=10, alphak = .5, seed = NULL)
}

\arguments{
  \item{data}{Continuous, Categorical ot Mixed data set}
  \item{nclus}{Number of clusters}
  \item{ndim}{Dimensionality of the solution}
  \item{method}{Specifies the method. Options are RKM for Reduced K-means, FKM for Factorial K-means, mixedRKM for Mixed Reduced K-means, mixedFKM for Mixed Factorial K-means, MCAk for MCA K-means, iFCB for Iterative Factorial Clustering of Binary variables and clusCA for Cluster Correspondence Analysis.}
  \item{scale}{A logical value indicating whether the metric variables should be scaled to have unit variance before the analysis takes place (default = \code{TRUE)}}
  \item{center}{A logical value indicating whether the metric variables should be shifted to be zero centered (default = \code{TRUE)}}
  \item{alpha}{Adjusts for the relative importance of (mixed) RKM and FKM in the objective function; \code{alpha = 1} reduces to PCA/PCAMIX, \code{alpha = 0.5} to (mixed) reduced K-means, and \code{alpha = 0} to (mixed) factorial K-means}
  \item{nstart}{Number of random starts (default = 100)}
  \item{nboot}{Number of bootstrap pairs of partitions}
  \item{alphak}{Non-negative scalar to adjust for the relative importance of MCA (\code{alphak = 1}) and K-means (\code{alphak = 0}) in the solution (default = .5). Works only in combination with \code{method = "MCAk"}}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when \code{smartStart = NULL}. The default value is NULL.}
}


\details{
The algorithm for assessing local cluster stability is similar to that in Hennig (2007) and can be summarized in three steps:

\emph{Step 1. Resampling:} Draw bootstrap samples S_i and T_i of size n from the data and use the original data as evaluation set E_i = X. Apply a joint dimension reduction and clustering method to S_i and T_i and obtain C^{S_i} and C^{T_i}.

\emph{Step 2. Mapping}: Assign each observation x_i to the closest centers of C^{S_i} and C^{T_i} using Euclidean distance, resulting in partitions C^{XS_i} and C^{XT_i}.

\emph{Step 3. Evaluation}: Obtain the maximum Jaccard agreement between each original cluster C_k and each one of the two bootstrap clusters, C_^k'{XS_i} and C_^k'{XT_i} as measure of agreement and stability, and take the average of each pair.

Inspect the distributions of the maximum Jaccard coefficients to assess the cluster level (local) stability of the solution.

Here are some guidelines for interpretation. Generally, a valid, stable cluster should yield a mean Jaccard similarity value of 0.75 or more. Between 0.6 and 0.75, clusters may be considered as indicating patterns in the data, but which points exactly should belong to these clusters is highly doubtful. Below average Jaccard values of 0.6, clusters should not be trusted. "Highly stable" clusters should yield average Jaccard similarities of 0.85 and above. 

While B = 100 is recommended, smaller run numbers could give quite informative results as well, if computation times become too high.

Note that the stability of a cluster is assessed, but stability is not the only important validity criterion - clusters obtained by very inflexible clustering methods may be stable but not valid, as discussed in Hennig (2007).
}

\value{
\item{nclus}{An integer with the number of clusters}
\item{clust1}{Partitions, C^{XS_i} of the original data, X, predicted from clustering bootstrap sample S_i (see Details)}
\item{clust2}{Partitions, C^{XT_i} of the original data, X, predicted from clustering bootstrap sample T_i (see Details)}
\item{index1}{Indices of the original data rows in bootstrap sample S_i}
\item{index2}{Indices of the original data rows in bootstrap sample T_i}
\item{Jaccard}{Mean Jaccard similarity values}
}
\references{
Hennig, C. (2007). Cluster-wise assessment of cluster stability. \emph{Computational Statistics and Data Analysis}, \emph{52}, 258-271.
}

\seealso{
\code{\link{global_bootclus}}
}

\examples{
## 5 bootstrap replicates and nstart = 10 for speed in example,
## use more for real applications
data(iris)
bootres = local_bootclus(iris[,-5], nclus = 3, ndim = 2,
method = "RKM", nboot = 5, nstart = 1, seed = 1234)

boxplot(bootres$Jaccard, xlab = "cluster number", ylab =
"Jaccard similarity")

## 5 bootstrap replicates and nstart = 5 for speed in example,
## use more for real applications
#data(diamond)
#bootres = local_bootclus(diamond[,-7], nclus = 4, ndim = 3,
#method = "mixedRKM", nboot = 5, nstart = 10, seed = 1234)

#boxplot(bootres$Jaccard, xlab = "cluster number", ylab =
#"Jaccard similarity")

## 5 bootstrap replicates and nstart = 1 for speed in example,
## use more for real applications
#data(bribery)
#bootres = local_bootclus(bribery, nclus = 5, ndim = 4,
#method = "clusCA", nboot = 10, nstart = 1, seed = 1234)

#boxplot(bootres$Jaccard, xlab = "cluster number", ylab =
#"Jaccard similarity")
}
