\name{plot_theta}
\alias{plot_theta}
\title{Allele Summary Plot}
\description{
 Creates an Allele Summary plot (allele B vs. allele A) for each SNP specified in snpsfilename. The points in the plot are the (theta A, theta B) ordered
 pairs for all the samples of the SNP. If a plotfilename is specified, 
 it will save the plot as a .ps file, otherwise the plot is shown on screen.
}
\usage{
plot_theta(genotypefile = "Xba.rlmm",
	   thetafile = "Xba.theta",
	   Pick.Obj = "FALSE",
	   plotfile = "plots.ps",
	   snpsfile = "snps.lst")
}
\arguments{
  \item{genotypefile}{Character string specifying the directory AND name of the .rlmm file created by \code{Classify} (optional)}
  \item{thetafile}{Character string specifying the directory AND name of the .rlmm file created by \code{Create_Thetafile} (optional)}
  \item{Pick.Obj}{At this point, it should always be left as the default FALSE, ie: it is for development purposes only (optional)}
  \item{plotfile}{The name where to store the plot as a .ps file, if blank such as "", it will display on screen instead (optional)}
  \item{snpsfile}{A list of SNPs to plot, with one SNP following another on a newline (optional)}
}

\author{Nusrat Rabbee <nrabbee@post.harvard.edu>, Gary Wong
  <wongg62@gmail.com>
}

\keyword{hplot}
