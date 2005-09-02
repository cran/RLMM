\name{create_Thetafile}
\alias{create_Thetafile}
\title{Calculating Parameter Estimates}
\description{
Assuming that the *.norm files are created, this step of the data analysis will calculate estimates of theta A and theta B values for each SNP and chip based on normalized probe intensity data from the *.norm files. The theta values are produced from fitting a probe-level additive model to the log2 A probe
intensities and the B intensities separately.
}
\usage{
create_Thetafile(probefiledir = getwd(),
                 start = 1,
                 end = -1,
                 thetafile = "")
}
\arguments{
\item{probefiledir}{Character string specifying the directory with the *.norm files (optional)}
\item{start}{An integer value specifying which SNP number we should start at when calculating the theta values (optional)}
\item{end}{An integer value specifying which SNP number we should stop at when calculating the theta values (optional)}
\item{thetafile}{A character string specifying the name the theta file will be saved as (optional)}
}

\author{Nusrat Rabbee <nrabbee@post.harvard.edu>, Gary Wong
  <wongg62@gmail.com>
}

\keyword{methods}
