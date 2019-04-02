% File man/arrayAppend.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{arrayAppend}
\alias{arrayAppend}
\title{Append two arrays by adding new data to top array}
\description{
\code{arrayAppend} is a function to append two arrays by adding new data to array whose
name is placed in the first argument of the function.}

\usage{
arrayAppend(topA, bottomA, missingVal = 5)
}

\arguments{
\item{topA}{Old data.}
\item{bottomA}{New data.}
\item{missingVal}{Value to use for missing SNPs in new data, default = 5.}
}

\details{
  Two matrices created by the function, \code{\link{toArray}}, can be appended together using
  \code{arrayAppend}. Appending is more efficient after recoding each of the two matrices by using
  \code{\link{snpRecode}}. The first input argument \code{topA} holds old data, i.e., data to hold on to its
  SNP list. The second argument, \code{bottomA}, holds new data that needs to be fit to old SNP list.
  Therefore, missing SNPs in new data will be created anew but will be given the value of \sQuote{missingVal},
  with a default of 5. This also assumes that you will be imputing these missing values later on.

  \code{arrayAppend} goes by the SNP list of the old data in the first input argument, so if the
  new data matrix has more SNPs, they will be ignored. If you need to consider the extra SNPs,
  switch the order of the input matrices of the function.
}

\value{
  An object of class matrix with \sQuote{number of rows = number of individuals in both matrices,
  topA and bottomA} and \sQuote{number of columns = number of SNPs as listed in
  old data}.
}

\seealso{
\code{\link{toArray}}
}

\examples{

ga.old <- matrix(sample(c(0,1,2,5), size = 30, repl = TRUE), nrow=3, ncol=10,
  dimnames = list(paste("Individual", 1:3, sep="."), paste("SNP", 1:10, sep=".")))

ga <- matrix(sample(c(0,1,2,5), size = 24, repl = TRUE), nrow=3, ncol=8,
dimnames = list(paste("Individual", 4:6, sep="."), paste("SNP", 1:8, sep=".")))

arrayAppend(ga.old, ga)
arrayAppend(ga.old, ga)

#
#             SNP.1 SNP.2 SNP.3 SNP.4 SNP.5 SNP.6 SNP.7 SNP.8 SNP.9 SNP.10
#Individual.1     5     5     0     1     2     1     2     5     1      2
#Individual.2     2     2     2     0     1     1     1     1     2      0
#Individual.3     1     1     5     0     5     5     0     5     5      5
#Individual.4     1     2     1     1     0     2     0     0     5      5
#Individual.5     2     1     2     1     0     1     2     1     5      5
#Individual.6     0     2     0     2     5     0     1     1     5      5
#
#Note that SNP.9 and SNP.10 were added to new data but were assigned the default
#code, 5, for missing values. You should impute all missing values afterwards.
#

#If the order of input matrices is switched, the extra SNPs in the second argument
#will be ignored. This function goes by the SNP list of the first arguement.

arrayAppend(ga, ga.old)

#             SNP.1 SNP.2 SNP.3 SNP.4 SNP.5 SNP.6 SNP.7 SNP.8
#Individual.4     1     2     1     1     0     2     0     0
#Individual.5     2     1     2     1     0     1     2     1
#Individual.6     0     2     0     2     5     0     1     1
#Individual.1     5     5     0     1     2     1     2     5
#Individual.2     2     2     2     0     1     1     1     1
#Individual.3     1     1     5     0     5     5     0     5

}

