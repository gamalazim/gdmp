% File man/BayesAB.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{BayesAB}
\alias{BayesAB}
\title{Estimte SNP Effects by MCMC - Bayes (A and B) - beta}
\description{
\code{BayesAB} is an implementation of Bayes A and B for estimating direct SNP effects in
high dimensional data problems (\emph{p >> N}). \code{BayesAB} utilizes the C Function \code{cBaysABF} for speed ..
}

\usage{
BayesAB(ga, numiter = 5000, numMHIter = 10, Pi = .9, y)
}
\arguments{
\item{ga}{A matrix of genotypes with a number of rows identical to the number of genotyped individuals
          and a number of columns identical to the number of SNPs. Values in the matrix are 0, 1, 2, & 5
          for homozygous, heterozygous, other homozygous, & unknown genotypes, respectively.}
\item{numiter}{Number of iterations}
\item{numMHIter}{Number of Metropolis-Hastings iterations. Set to a value of 1 with Bayes A.}
\item{Pi}{Proportion of SNP loci with 0 effect. Set to 0 to run Bayes A instead of B.}
\item{y}{Trait phenotypes or conventional breeding values}
}


\details{
This function runs Bayes A and B to estimate direct SNP effects based on a matrix of genotypes,
\code{ga} and a vector of adjusted phenotypes, \code{y}, (Meuwissen et al., 2001; Genetics 157:1819-1829).
This is essential in high dimensional data problems with highly overparameterized models (\emph{p >> N}).
To run Bayes A, set \code{Pi} to 0 and \code{numMHIter} to 1.

Other data management functions in \code{gdmp} can be used to construct the integer matrix of genotypes,
\code{ga}, to use as input to \code{BayesAB}.

}


\value{

  A list object with a vector of SNP estimates, \code{meanb}, and a vector of genomic values for individuals,
  \code{aHat}, are returned. It is also possible to extract the number of SNP loci in \code{nLoci}.
}

\seealso{
\code{\link{BayesCpi}}
}
\references{
  Meuwissen et al. (2001).
  Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps.
  \emph{Genetics}, \bold{157}, 1819--1829.
}

