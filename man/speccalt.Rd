\name{speccalt}
\alias{speccalt}
\title{
speccalt
}
\description{
Spectral clustering algorithm, with optional automatic estimation of the optimal number of clusters.
}
\usage{
speccalt(kern, k=NA, maxk=20)
}
\arguments{
  \item{kern}{
A kernel matrix object, as returned e.g. by local.rbfdot.
  }
  \item{k}{
Parameterized number of clusters. If NA, the automatic estimation procedure (bartlett) is used.
	}
	\item{maxk}{
If k is NA, an upper bound for the automatic estimation. Defaults to 20.
}
}
\value{
A numeric vector of cluster labels. Clusters are indexed from 1 to k.
}
\seealso{
local.rbfdot
bartlett
}
\author{
Pierrick Bruneau
}
\examples{
kern <- local.rbfdot(synth2)
clust <- speccalt(kern) # with automatic estimation
clust2 <- speccalt(kern, 4) # without automation
}

