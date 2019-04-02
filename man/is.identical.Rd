% File man/is.identical.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{is.identical}
\alias{is.identical}
\title{Check two SNPs for near identity across genotyped samples}
\description{
Given individual genotypes of two SNPs, the function checks if the two sets of
genotypes are identical or completely opposite. A margin of error of 0.5\% is allowed.
}


\usage{
is.identical(x, y, allow = .005)
}

\arguments{
  \item{x, y}{two column vectors in the genotypes array, created by \code{toArray}
  and converted to integer genotypes by \code{snpRecode}.
  The two columns are genotypes of two SNPs for all individuals in data.}
  \item{allow}{allowed margin of error}
}

\details{
Test if two SNPs are identical in genotypes across a set of individuals. SNPs are considered
identical if the number of different genotypes in the population tested remains below
an allowed error margin of 0.5\%.

In addition to identical SNPs, the function considers SNP genotypes that are entirely
opposite within error margin as redundant. Thus, two SNPs are declared highly correlated
if the genotypes are all the same (0-0, 1-1, and 2-2) or are all opposite (0-2, 1-1, 2-0) within
the error margin specified.
}

\value{
A logical value, of TRUE for identical SNPs or FALSE for different SNPs, is returned.
}

\seealso{
\code{\link{GetHCS}}, \code{\link{snpRecode}}, \code{\link{toArray}}
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

## Check the first 2 SNPs for being identical based
## on a minimum allowed margin of error of 0.5%
is.identical(ga.r[,1], ga.r[,2], allow = .005)
# [1] FALSE

## Create an instance of exactly identical SNP genotypes
ga.r <- cbind(ga.r, ga.r[,1])  # SNP #1 and #101 are exactly identical
is.identical(ga.r[,1], ga.r[,101], allow = 0)
# [1] TRUE

## Create an instance of identical SNP genotypes with a 5% error
ga.r <- cbind(ga.r, ga.r[,1])  # SNP #1 and #101 are 100% identical
ga.r[20,101] <- 2  # a different genotype, to make SNP #1 & #101 only 95% identical
is.identical(ga.r[,1], ga.r[,101]) # use default allow of .005
# [1] FALSE

is.identical(ga.r[,1], ga.r[,101], allow = .05) # allow for a 5% marging of error
# [1] TRUE
}

