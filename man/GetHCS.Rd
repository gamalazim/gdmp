% File man/GetHCS.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{GetHCS}
\alias{GetHCS}
\title{Check SNP genotypes and exclude redundant SNPs}
\description{
Given individual genotypes of a set of SNPs, the function checks for the existence
of redundant SNPs to exclude. A margin of error of 0.5\% is allowed by default.
}

\usage{
GetHCS(ga.r, Exclude=1:ncol(ga.r), allow = .005)
}

\arguments{
  \item{ga.r}{Matrix of genotypes created by \code{toArray} and converted to integer
  genotypes matrix by \code{snpRecode}.}
  \item{Exclude}{Numeric indices or column names of \code{ga.r} to check.}
  \item{allow}{allowed margin of error, default is 0.005.}
}

\details{
Test for similar SNP genotypes across a set of individuals. SNPs are considered
identical if the number of different genotypes in the population tested remains below
an allowed error margin of 0.5\%. Say, \code{Exclude <- 1:100} with SNP #1 similar to #25,
then \code{Exclude[25]} will be flagged for exclusion, whereas \code{Exclude[1]}
will not be flagged for exclusion.

In addition to identical SNPs, the function flaggs SNP genotypes that are entirely
opposite within error margin as redundant as well. Thus, SNPs are declared highly correlated
if the genotypes are all the same (0-0, 1-1, and 2-2) or all opposite (0-2, 1-1, 2-0) within
the error margin specified.
}

\value{
If \code{Exclude} contains SNP names, a character
vector of excluded SNPs is returned, and if it contains integer values, a numeric vector
of excluded SNPs is returned.
}

\seealso{
\code{\link{is.identical}}, \code{\link{snpRecode}}, \code{\link{toArray}}
}

\examples{

## Simulate random allele designations for 100 bi-allelic SNPs
set.seed(2016)
desig <- array(sample(c('A','C','G','T'), size = 200, repl = TRUE), dim=c(100, 2))

## Simulate random SNP genotypes for 20 individuals - put them in array format
## '-' indicates an unknown base
ga <- array(0, dim=c(20, 100))
for(i in 1:20)
  for(j in 1:100)
    ga[i, j] <- paste(sample(c(desig[j,],"-"), 2, prob=c(.47, .47, .06), repl=TRUE), collapse='')

## Recode the matrix, place recoded genotypes in ga.r
desig <- data.frame(AlleleA_Forward = factor(desig[,1]), AlleleB_Forward = factor(desig[,2]))
ga.r <- array(5, dim=c(20, 100))
for(i in 1:100) ga.r[,i] <- snpRecode(ga[,i], desig[i,])

## Check all SNP genotypes in ga.r for similarity across individuals
## Allow for a margin of error of 0.5%
GetHCS(ga.r)
#[1] 42 91  # SNPs 42 & 91 are similar to earlier SNPs in the vector, 'Exclude'

## Check SNP genotypes from 1 to 50 for similarity across individuals
GetHCS(ga.r, Exclude=1:50)
#[1] 42

}

