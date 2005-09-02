\name{normalize_Rawfiles}
\alias{normalize_Rawfiles}
\title{Normalize PM Intensity values}
\description{
Given a directory with *.raw files, it will normalize the PMA and PMB intensities in 
each file using Xba.CQV (composite quantile vector) and return the normalized 
values written to *.norm files corresponding to its *.raw files. 
EG: If two *.raw files are used, two *.norm files will be returned.
This normalization simply puts the probe data on the same scale
as the training data.
}
\usage{
normalize_Rawfiles(cqvfile = "",
                   probefiledir = getwd())
}
\arguments{
\item{cqvfile}{Character string specifying the CQV filename (e.g., Xba.CQV) (required)}
\item{probefiledir}{Character string specifying location of the *.raw files and *.norm files (optional)}
}

\author{Nusrat Rabbee <nrabbee@post.harvard.edu>, Gary Wong
  <wongg62@gmail.com>
}

\keyword{manip}
