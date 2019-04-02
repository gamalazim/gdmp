% File man/toArray.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{toArray}
\alias{toArray}
\title{Turn Illumina's Genome Studio \sQuote{Final Report} file into an array}
\description{
\code{toArray} is a function to turn genotyping data into an array.
  This function does not read stored data in the \sQuote{Final Report} file, however,
  it converts a data fram with the same structure as \sQuote{Final Report} into an object of class matrix.
  The specifications of the genotyping data frame are given under \sQuote{Details}.
}

\usage{
toArray(finalRep)
}
\arguments{
\item{finalRep}{A data frame with the same structure as \sQuote{Final Report}, containing
     genotyping data. Specifications are given under \sQuote{Details}.}
}

\details{
Used to turn Illumina's Genome Studio \sQuote{Final Report} file into an array.
Specifications of the input data fram, \sQuote{finalRep}, are:
\itemize{
  \item Input is a data frame with its first 4 columns listed in the following order:
     \enumerate{
       \item SNP names, as factor, with equal number of SNPs per individual.
       \item Identification codes, as factor, for genotyped individuals.
       \item Allele 1 (one character: A, C, G, T, or -).
       \item Allele 2 (one character: A, C, G, T, or -).
     }
  \item All SNPs of individual 1 are listed first followed by SNPs of individual 2, and so on.
  \item SNPs are listed for each individual in the same order.
}

Note that it is easy to read the \sQuote{Final Report} file into a data frame which is then
used as an input to \code{toArray}. See \sQuote{Examples}.
}

\value{

  An object of class matrix with \sQuote{number of rows = number of individuals} and
  \sQuote{number of columns = number of SNPs}. Each element of the matrix consists of a
  two-character string for the two DNA bases of a single SNP locus.
}

\seealso{
\code{\link{arrayAppend}}
}

\examples{

## Read file './Final.Report', located in the current working directory,
## and place the first 4 columns in a data frame.

#d <- read.table("./Final.Report", skip=10)[,1:4]

## Use toArray to turn data read into a matrix

#ga <- toArray(finalRep = d)
#ga[1:6, 4000:4002]

#
#               SNP.4000           SNP.4001           SNP.4002
#Individual.1       "GG"               "CC"               "CC"
#Individual.2       "TG"               "AC"               "CC"
#Individual.3       "TG"               "AA"               "CC"
#Individual.4       "GG"               "AC"               "TC"
#Individual.5       "GG"               "AC"               "CC"
#Individual.6       "GG"               "AA"               "CC"
#
}

