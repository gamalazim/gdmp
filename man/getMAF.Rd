% File man/getMAF.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{getMAF}
\alias{getMAF}
\title{Return frequency of the minor allele of a SNP}
\description{
Given individual genotypes of a single SNP, the function returns allele frequency
of the minor allele.
}

\usage{
getMAF(snpG)
}

\arguments{
  \item{snpG}{is a column vector in the genotypes array, created by \code{toArray}
  and converted to integer genotypes by \code{snpRecode}.
  The column represents genotypes of a single SNP for all individuals in data.}
}

\details{
Allele frequency is calculated for one of the two alleles and is returned if it is
below 0.5, otherwise \sQuote{1 - allele frequency} is returned. The function retuns
a frequency of 0 for non polymorphic SNPs.
}

\value{
Allele frequency of the minor allele.
}

\seealso{
\code{\link{snpSelect}}, \code{\link{snpRecode}}, \code{\link{toArray}}
}

\references{
  Falconer and Mackay (1996).
  Introduction to Quantitative Genetics (4th Edition).
  \emph{Pearson Education Limited, Edinburgh, England}.
}

\examples{
set.seed(002016)
snpG.1 <- c(rep(0, 100), rep(1, 80), rep(2, 9))
snpG.2 <- c(rep(0, 100), rep(1, 80))
snpG.3 <- c(rep(0, 100), rep(2, 9))
snpG.4 <- c(rep(0, 100))

getMAF(snpG.1)
#        0
#0.2592593

getMAF(snpG.2)
#        1
#0.2222222

getMAF(snpG.3)
#         0
#0.08256881

getMAF(snpG.3)

}
