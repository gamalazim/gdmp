% File man/read.findhap.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{read.findhap}
\alias{read.findhap}
\title{Read the output of Paul Vanraden's findhap, \sQuote{genotypes.filled}}
\description{
\code{read.findhap} is a tool that utilizes external C code to read quickly the output of
  findhap program stored by default in the \sQuote{genotypes.filled} file. \sQuote{findhap}
  is a standalone software program written by Paul Vanraden to impute missing genotypes. If you
  are not familiar with \sQuote{findhap}, you probably don't need this function.
}

\usage{
read.findhap(Nanim, Nmark, file="./genotypes.filled")
}
\arguments{

\item{Nanim}{Total number of individuals to read from imputed genotypes file.}
\item{Nmark}{Total number of SNP genotypes in the imputed-to chip.}
\item{file}{Findhap's output file for imputed genotypes. Default is \sQuote{genotypes.filled}.}

}


\details{
This function reads into R environment imputed SNP genotypes as formated by findhap program. This tool
utilizes external C code for quick read. The standard output of this function is an object of class matrix.

}

\value{

  An object of class matrix with \sQuote{number of rows = Nanim} and
  \sQuote{number of columns = Nmark + 3}. Refer to findhap manual for
  more details about the structure of \sQuote{genotypes.filled}.
}

\references{
  Paul Vanraden.
  Find haplotypes and impute genotypes using multiple chip sets and sequence data.
  http://aipl.arsusda.gov/software/findhap/.
}

