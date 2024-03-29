% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{clusteringSummary}
\alias{clusteringSummary}
\title{Matrix of cluster characteristics}
\usage{
clusteringSummary(
  gg,
  att = c("lec", "wt", "fc", "infomap", "louvain", "sgG1", "sgG2", "sgG5", "spectral")
)
}
\arguments{
\item{gg}{graph to analyse}

\item{att}{vector of attribute names that contains membership data}
}
\value{
matrix of clustering characteristics
}
\description{
Matrix of cluster characteristics
}
