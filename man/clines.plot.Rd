\name{clines.plot}
\alias{clines.plot}
\title{Clines Plot}
\description{
This function produces graphical plots of genomic clines using the
output from the function \code{genomic.clines}. 

}
\usage{
clines.plot(cline.data=NULL,marker.order=NULL,rplots=3,cplots=3,pdf=TRUE,
            out.file="clines.pdf")
}
\arguments{
  
  \item{cline.data}{a list that is the product of the
    \code{genomic.clines} function.}
  
  \item{marker.order}{an optional numeric or character vector specifying
    the order to plot marker results in; if \code{marker.order = NULL}
    markers are plotted in the order they were originally provided.}
  
  \item{rplots}{numerical value specifying the number of plots per row
    in the output.}
  
  \item{cplots}{numerical value specifying the number of plots per
    column in the output.}
  
  \item{pdf}{logical specifying whether to print the plots to a pdf file;
    if \code{pdf=FALSE} plots are printed to the current graphical
    device.}
  
  \item{out.file}{a character string for the filename for the
    output if \code{pdf=TRUE}.}
}  

\details{
  
This function produces graphical plots based on the data object produced
by the \code{genomic.clines} function. A separate plot is produced for
each marker. Plots depict the probability of a given genotype as a
function of hybrid index. Plots for the observed data are shown with
solid (homozygous population 1 genotype for co-dominant markers, or
population 1 allele for dominant or haploid markers) and dashed lines
(inter-specific heterozygotes for co-dominant markers, not shown for
dominant or haploid markers). If significance testing was conducted with
the \code{genomic.clines} function, the corresponding 95\% confidence
intervals are shown in dark and light green. Circles indicate
individuals that are homozygous for alleles from population 1
(co-dominant) or population 1 alleles (dominant or haploid) on the top
line, interspecific heterozygotes (co-dominant only) in the middle, and
homozygotes for alleles from population 2 (co-dominant) or population 2
alleles (dominant or haploid) on the bottom line.  These number of
individiuals with each genotype is printed on the right vertical
axis. The locus name and \emph{P}-value are printed in each plot.
  
See Gompert and Buerkle (2009a, 2009b) for additional details and examples.
}

\value{
A plot is produced, but a data object is not returned.  
}
  
\references{
Gompert Z. and Buerkle C. A. (2009) A powerful regression-based
method for admixture mapping of isolation across the genome of
hybrids. \emph{Molecular Ecology}, \emph{in press}.

Gompert Z. and Buerkle C. A. (2009) \emph{introgress}: a software
package for mapping components of isolation in hybrids. \emph{Molecular
  Ecology Resources}, \emph{in preparation}.  

 } 

\author{
Zachariah Gompert \email{zgompert@uwyo.edu}, C. Alex Buerkle \email{buerkle@uwyo.edu}
}  

\seealso{
\code{\link{genomic.clines}}
 } 
 
\examples{
## this code assumes the data object clines.out2 has been produced using
## the sample code for the "genomic.clines" function

## produce plots of genomic clines 
\dontrun{
clines.plot(cline.data=clines.out2)
}
} 
