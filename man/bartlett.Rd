\name{bartlett}
\alias{bartlett}
\title{
bartlett
}
\description{
Adaptation of the Bartlett statistical test for equal variances, to a procedure of detection of the optimal number of clusters, in the context of spectral clustering.
}
\usage{
bartlett(eigvals, thres=0.95, maxk=20)
}
\arguments{
  \item{eigvals}{
A vector of eigenvalues, as returned by the eigendecomposition of the normalized Laplacian (see examples), i.e. in decreasing order.
	}
	\item{thres}{
1 - alpha significance threshold, defaults to 95%.
	}
	\item{maxk}{
upper bound to the returned value.
	}
}
\value{
The computed optimal number of clusters.
}
\author{
Pierrick Bruneau
}
\examples{
	# normalized Laplacian computation
	kern <- local.rbfdot(synth1)
	nelts <- dim(kern)[1]
	diag(kern) <- 0
	deg <- sapply(1:nelts, function(i) {
		return(sum(kern[i,]))
	})
	deg <- diag(deg)
	L <- diag(nelts) - solve(deg) %*% kern
	eig <- eigen(L)

	optK <- bartlett(eig$values)$k
}

