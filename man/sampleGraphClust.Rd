% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{sampleGraphClust}
\alias{sampleGraphClust}
\title{Perturbe graph and calculate its clustering}
\usage{
sampleGraphClust(gg, mask, alg, type, reclust = FALSE, Cnmin = -1, Cnmax = 10)
}
\arguments{
\item{gg}{graph}

\item{mask}{percentage of elements to perturbe}

\item{alg}{clustering alg.}

\item{type}{edges=>1 or nodes=>2  to mask}

\item{reclust}{logical to decide wether to invoke reclustering}

\item{Cnmin}{Cn min for Spectral algorithm}

\item{Cnmax}{Cn max for reclustering algorithms}
}
\description{
Perturbe graph and calculate its clustering
}
