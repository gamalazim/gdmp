% File man/is.hwEq.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{is.hwEq}
\alias{is.hwEq}
\title{Check a SNP for Hardy-Weinberg equilibrium}
\description{
Given individual genotypes of a single SNP, the function checks a population of
individuals for Hardy-Weinberg equilibrium.
}


\usage{
is.hwEq(snpG, diff)
}

\arguments{
  \item{snpG}{is a column vector in the genotypes array, created by \code{toArray}
  and converted to integer genotypes by \code{snpRecode}.
  The column represents genotypes of a single SNP for all individuals in data.}
  \item{diff}{heterozygosity difference used in HW equilibrium, see \sQuote{Details}.}
}

\details{
A logical function to check for HW equilibrium, the expected frequency of the heterozygous
genotype is estimated based on the two homozygous genetypes by \sQuote{sqrt(4*p^2*q^2)}.
The absolute difference between the expected and observed frequency of the heterozygous
genotype of the SNP needs to be smaller than the minimum difference of \code{diff}.
}

\value{
A logical TRUE for H-W equilibrium or FALSE for H-W disequilibrium is returned.
}

\seealso{
\code{\link{snpSelect}}, \code{\link{snpRecode}}, \code{\link{toArray}}
}

\references{
  Falconer and Mackay (1996).
  Introduction to Quantitative Genetics (4th Edition).
  \emph{Pearson Education Limited, Edinburgh, England}.

  Wiggans et al. (2009).
  Selection of single-nucleotide polymorphisms and quality of genotypes
  used in genomic evaluation of dairy cattle in the United States and Canada.
  \emph{Journal of Dairy Science}, \bold{92}, 3431-3436.
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
    ga[i, j] <- paste(sample(c(desig[j,],"-"), 2, prob=c(.47, .47, .06), repl = TRUE), collapse='')

## Recode the matrix, place recoded genotypes in ga.r
desig <- data.frame(AlleleA_Forward = factor(desig[,1]), AlleleB_Forward = factor(desig[,2]))
ga.r <- array(5, dim=c(20, 100))
for(i in 1:100) ga.r[,i] <- snpRecode(ga[,i], desig[i,])

## Check the first 10 SNPs for H-W equilibrium based on a minimum
## allowed difference of 0.15 between observed and expected heterozygosity
apply(ga.r[,1:10], 2, is.hwEq, diff=.15)
# [1]  TRUE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE  TRUE FALSE
}

