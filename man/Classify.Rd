\name{Classify}
\alias{Classify}
\title{Classification of SNPs based on theta estimates}
\description{
This function entails classification of SNPs based on the theta estimates (thetafile), genotype information (A regions file), and some internal files.
Currently, this algorithm works for the Affymetrix 100K - Xba dataset.
}
\usage{
Classify(genotypefile = "",
         regionsfile = "",
	 thetafile = "",
	 callrate = 100)
	}
\arguments{
\item{genotypefile}{Name of the classified SNPs with the genotypes (required)}
\item{regionsfile}{Character string specifying the directory AND name of regionsfile - e.g., "Xba.regions" (required)}
\item{thetafile}{Character string specifying the directory AND name of thetafile (required)}
\item{callrate}{Call Rate percentage; The user can specify any number from the list: 80,82,84,86,88,90,92,94,96,98,100. Default is 100\%(optional)}
}

\details{
For each SNP, Mahalanobis distances from each chip's (theta A, theta B) ordered pair to the genotype centers is calculated. Each chip is assigned the genotype of the cluster which it is closest to (ie: AA, AB, BB).
}

\author{Nusrat Rabbee <nrabbee@post.harvard.edu>, Gary Wong
  <wongg62@gmail.com>
}

\keyword{classif}
