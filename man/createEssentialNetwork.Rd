% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createEssentialNetwork.R
\name{createEssentialNetwork}
\alias{createEssentialNetwork}
\title{Create a synthetic lethal network}
\usage{
createEssentialNetwork(
  datapath = ".",
  outpath = ".",
  method = "ISLE",
  n = 0,
  tissue = "all"
)
}
\arguments{
\item{datapath}{Path to directory containing Dependency Map data}

\item{outpath}{Output path to save .rds object, optionally including filename.}

\item{method}{Character string describing which method to use. One of c("ISLE")}

\item{n}{The number of genes to compute associations for, useful for debugging and testing.
Default is 0, which computes the entire network.}

\item{tissue}{Character string defining which tissue to use. Default is "all"}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{    A dataframe containing a matrix of effect sizes (effectSizes), a matrix of p-values
}\if{html}{\out{</div>}}

(pvals), and a list of network parameters.
}
\description{
createEssentialNetwork takes as arguments a path to Cancer Dependency Map essentiality data
and parameters defining the network properties. It then computes a gene essentiality
network. All the variations produce a test statistic and a p-value.
}
