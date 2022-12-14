% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAplotInference.R
\name{CAplotInference}
\alias{CAplotInference}
\title{CAplotInference run a CA
(with the \code{ExPosition} package)
and generates the inference graphs and tables.
Note: \emph{Still Under Development}.}
\usage{
CAplotInference(
  resCA,
  data,
  DESIGN = NULL,
  make_design_nominal = TRUE,
  k = 0,
  graphs = 12,
  printGraphs = FALSE,
  col4I = NULL,
  col4J = NULL,
  niter.boot = 100,
  niter.perm = 100,
  fast = FALSE,
  save2pptx = FALSE,
  title4pptx = "CA Inference Results"
)
}
\arguments{
\item{resCA}{Output from epCA}

\item{data}{A data frame or a matrix with
data suitable for a CA.}

\item{DESIGN}{Default: NULL.
A design vector (could be factor or character)
or (Boolean) matrix used to assign rows
to groups.}

\item{make_design_nominal}{if TRUE (Default) transform
the vector from \code{DESIGN} into
a Boolean matrix.
Passed to
\code{InPosition::epCA.inference.battery}.}

\item{k}{number
of factor to keep; when equql to
0  (Default), all factors are kept.
Passed to
\code{InPosition::epCA.inference.battery}.}

\item{graphs}{do we want graphs? At the moment, this parameter is ignored.
Current Default is \code{12} which indicates that
the graphs are generated for the first
2 components. Note that current version
is creating output only for the first two
components.}

\item{printGraphs}{(Default: FALSE)
do we want to print the graphics as \code{.png}?}

\item{col4I}{a color vector for
plotting the rows (if \code{NULL}
Default) use colors from \code{ExPosition::epCA}.}

\item{col4J}{a color vector for
plotting the columns (if \code{NULL}
Default) use colors from \code{ExPosition::epCA}.}

\item{niter.boot}{Number of bootstrap samples to run. Default: 100}

\item{niter.perm}{Number of permutations to run. Default: 100}

\item{fast}{Use the fast version of the perm and bootstrap from data4PCCAR}

\item{save2pptx}{Default: FALSE}

\item{title4pptx}{Title of the PPTX, Default:
'CA Inference Results'.}
}
\value{
A list made of two lists
}
\description{
Generates the inference graphs and tables for CA.
}
\details{
Work in Progress
}
\author{
Luke Moraglia
}
