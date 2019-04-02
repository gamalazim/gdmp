% File man/snpSelect.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{snpSelect}
\alias{snpSelect}
\title{Select genotyped individuals and SNPs based on call rate, heterozygosity, HW equilibrium,
  or minor allele frequency.}
\description{
\code{snpSelect}
is a function to exclude individuals and SNPs that do not satisfy one of the following criteria,
\itemize{
 \item SNP genotypes are above the provided call rate threshold, default = 0.85
 \item SNP genotypes are above the provided heterozygosity rate, default = 0.01
 \item SNP loci are not far from Hardy-Weinberg equilibrium
 \item SNP genotypes are above the provided minimum for minor allele frequency, default = 0.005
  }

}

\usage{
snpSelect(ga.r, select.method=c("call.rate", "heterozygosity", "HW.Eq", "MAF"),
  call.rate.min = .85, hz.min = .01, hz.diff = .15, MAF.min = .005)
}

\arguments{
  \item{ga.r}{Matrix of SNP genotypes as prepared by \code{toArray} and recoded by \code{snpRecode}.}
  \item{select.method}{The SNP criterion to use in selecting individuals to keep, default is 'call.rate'.}
  \item{call.rate.min}{Call rate threshold provided or the minimum allowed call rate.}
  \item{hz.min}{Heterozygosity threshold provided, less polymorphic SNPs are excluded.}
  \item{hz.diff}{Heterozygosity difference used in HW equilibrium, see \sQuote{Details}.}
  \item{MAF.min}{Minor allele frequency threshold provided, SNPs with lower frequencies
    for the minor allele are excluded.}
}

\details{
This is a selection tool to extract high quality individuals and SNPs. Based on the
selection criterion provided, \code{snpSelect} keeps individuals and SNPs that meet the criterion.
Providing multiple selection criteria to \code{select.method} does not work. If a cleanup
requires extraction of individuals and SNPs based on multiple criteria, the function may
be used iteratively.

In the case of \code{select.method} = \sQuote{call.rate}, individuals not SNPs are excluded.
In all other cases, SNPs are excluded.

For selection based on HW equilibrium, the expected frequency of the heterozygous
genotype is estimated based on the two homozygous genetypes by \sQuote{sqrt(4*p^2*q^2)}.
The absolute difference between the expected and observed frequency of the heterozygous
genotype of each SNP needs to be smaller than the minimum difference of \code{hz.diff}.
The logical function \code{is.hwEq(x, diff)} with a TRUE or FALSE output can be used to
test SNP genotypes in the integer vector, \code{x}, based of the difference, \code{diff}.
}

\value{
A matrix with individuals exceeding the \code{call.rate.min} threshold or a matrix with
SNPs exceeding \code{hz.min}. A minimum heterozygosity instead of 0 is recommended to leave
room for genotying errors. The function can also be used to exclude SNPs not in
Hardy-Weinberg equilibrium if the difference between expected and observed genotype
frequencies is greater than \code{hz.min}. Finally, \code{snpSelect} may be used to exclude SNPs
with low minor allele frequencies. If \code{MAF.min} is not supplied, the default of 0.5\% is
used.
}

\seealso{
\code{\link{getMAF}}, \code{\link{is.hwEq}}, \code{\link{snpRecode}}, \code{\link{toArray}}
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
    ga[i, j] <- paste(sample(c(desig[j,],"-"), 2, prob=c(.47, .47, .06), repl=TRUE), collapse='')

## Recode the matrix, place recoded genotypes in ga.r
desig <- data.frame(AlleleA_Forward = factor(desig[,1]), AlleleB_Forward = factor(desig[,2]))
ga.r <- array(5, dim=c(20, 100))
for(i in 1:100) ga.r[,i] <- snpRecode(ga[,i], desig[i,])

## Exclude individuals with call rates below 85%
dim(snpSelect(ga.r, select.method = "call.rate", call.rate.min = 0.85))
#[1]  15 100

## Exclude SNPs with heterozygosity <= 1%
dim(snpSelect(ga.r, select.method = "heterozygosity", hz.min = 0.01))
#[1] 20 79

## If the difference between expected and observed genotype frequencies
## is > 0.15, exclude the SNP.
dim(snpSelect(ga.r, select.method = "HW.Eq", hz.diff = 0.15))
#[1] 20 29

## Select SNPs with minor allele frequencies greater than 0.5%
dim(snpSelect(ga.r, select.method = "MAF", MAF.min = 0.005))
#[1] 20 79
}

