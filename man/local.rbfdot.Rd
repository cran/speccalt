\name{local.rbfdot}
\alias{local.rbfdot}
\title{
local.rbfdot
}
\description{
Locally adapted RBF kernel function. Taken and adapted from kernlab.
}
\usage{
local.rbfdot(data, n=5)
}
\arguments{
  \item{data}{
A data set matrix, with N elements as rows.
}
	\item{n}{
Parameter for median computation, defaults to 5.
}
}
\value{
A NxN matric of RBF kernel values.
}
\author{
Pierrick Bruneau
}
\examples{
kern <- local.rbfdot(synth1)
}

