% File man/snpRecode.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{snpRecode}
\alias{snpRecode}
\title{Recode a matrix of SNP genotypes as 0, 1, and 2}
\description{
\code{snpRecode} is a function to convert SNP genotypes to 0, 1, and 2 for the homozygous,
  heterozygous, and other homozygous genotype, respectively.}

\usage{
snpRecode(snpG, designat)
}

\arguments{
\item{snpG}{is a column vector in the genotypes array, created by \code{toArray}.
  The column represents genotypes of a single SNP for all or a subset of individuals in data.}
\item{designat}{is the 2-base allele designations for each SNP. This is sometimes
  called allele report data, where the specefic bases of alleles A and B are reported. Formated as
  data frame with two factors for alleles A and B. See \sQuote{Examples}.}
}

\details{
Recode snp genotypes by counting the number of copies of allele A in an element of \code{snpG}
which is a column vector in the genotypes array, \code{ga}, where
\itemize{
  \item \code{snpG} is a column vector in the genotypes array,
  \item \code{ga} is the genotypes array created by \code{toArray}. It contains elements such as "AA", "AG", "GA", "-A", "- -".
}
Unknown genotypes are those with non A/G/C/T bases, those are coded as 5.
}

\value{
A column vector of the integers 0, 1, and 2 is created based on the number of copies of allele A
in each element of the supplied vector of genotypes. A value of 5 is used to indicate an unknown
genotype.
}

\seealso{
\code{\link{toArray}}
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
    ga[i, j] <- paste(sample(c(desig[j,],"-"), 2, prob=c(.46, .46, .08), repl=TRUE), collapse='')

## Recode the matrix, place recoded genotypes in ga.r
desig <- data.frame(AlleleA_Forward = factor(desig[,1]), AlleleB_Forward = factor(desig[,2]))
ga.r <- array(5, dim=c(20, 100))
for(i in 1:100) ga.r[,i] <- snpRecode(ga[,i], desig[i,])

## Tabulate recoded genotypes in the matrix ga.r
table(ga.r)
#   0   1   2   5
# 326 632 701 341

}

