% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bam_classes.R
\docType{class}
\name{bioindex-class}
\alias{bioindex-class}
\alias{bioindex}
\alias{biodiversity_index}
\title{Class  \code{bioindex}}
\value{
An object of class bioindex
}
\description{
Class  \code{bioindex}
}
\section{Slots}{

\describe{
\item{\code{alpha}}{A matrix with the richness of species per site}

\item{\code{omega}}{A matrix with the range size of every species}

\item{\code{wBeta}}{A numeric value with Whittaker’s multiplicative beta index}

\item{\code{laBeta}}{A numeric value with Lande’s additive beta index}

\item{\code{leBeta}}{A numeric value with Legendre’s beta index}

\item{\code{nestedness}}{A numeric value with Wright and Reeves' nestedness}

\item{\code{dispersion_field}}{A matrix with the set of ranges of all species that
occur in at each locality}

\item{\code{richness_field}}{A matrix with the number of  shared species in
each site}
}}

\author{
Luis Osorio-Olvera & Jorge Soberón
}
