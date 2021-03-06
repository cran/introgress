\name{LociDataSim1}
\docType{data}
\alias{LociDataSim1}
\title{Marker Information for Simulation 1}
\description{
This data object provides the marker information associated with
\code{AdmixSimData1}. Rows correspond to markers (110 markers
spaced evenly across 10 linkage groups at 10 cM intervals). The first
column gives the locus name,  (c1.f is the marker at the location of the
genomic region that determined each
indvidual's fitness, see \code{AdmixSimData1}). The second column (type)
indicates that all of the markers are co-dominant (C). The third and
fourth columns specify each markers linkage group (between 1 and 10) and
relative location on each linkage group, respectively.
}
\usage{data(LociDataSim1)}
\format{A matrix with 110 rows and 4 columns.}
\source{From simulations conducted by CAB and ZG.}
\references{
Gompert Z. and Buerkle C. A. (2009) A powerful regression-based
method for admixture mapping of isolation across the genome of
hybrids. \emph{Molecular Ecology}, \bold{18}, 1207-1224.

Gompert Z. and Buerkle C. A. (2009) \emph{introgress}: a software
package for mapping components of isolation in hybrids. \emph{Molecular
  Ecology Resources}, \emph{in preparation}.  
} 

\seealso{
\code{\link{AdmixDataSim1}}
}

\keyword{datasets}

