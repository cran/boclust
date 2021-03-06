% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clust_Method.R
\name{OverlapClust}
\alias{OverlapClust}
\title{Overlap Clustering}
\usage{
OverlapClust(data.simi, p = c(0.9, 0.75, 0.5), lin = 0.25,
  pro.show = FALSE)
}
\arguments{
\item{data.simi}{The similarity matrix of Bossa scores obtained by \code{\link{BossaSimi}}}

\item{p}{A set of quantiles(90%, 75% and median) of the positive values of
similarity matrix to form clusters at different levels of within-cluster similarity.}

\item{lin}{A tuning parameter to control the size of each overlap cluster before
merging, smaller lin leads to larger cluster size.}

\item{pro.show}{A logical indicator whether show the details of the process.}
}
\value{
A list including two data.frame: overlap sub-clusters and cluster center for each.
}
\description{
Do overlap clustering with Bossa similarity in different levels of \code{p}.
}
\examples{
{
data(bo.simu.data)

## calculate the similarity matrix
bossa.simi <- BossaSimi(bo.simu.data)$bossa.simi

overlap.clust <- OverlapClust(bossa.simi)
}
}
