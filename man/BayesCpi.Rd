% File man/BayesCpi.Rd
% gdmp R package
% Copyright 2016 Gamal Abdel-Azim
% Distributed under GPL 2 or later

\name{BayesCpi}
\alias{BayesCpi}
\title{Estimte SNP Effects by MCMC - Bayes C and Cpi - beta}
\description{
\code{BayesCpi} is an implementation of Bayes Cpi to extend Bayes A and B for estimating direct SNP effects in
high dimensional data problems (\emph{p >> N}). \code{BayesCpi} treats the prior probability, \code{pi} = \emph{P(SNP has zero effect)}, as unknown.
The C Function \code{cBaysCpi} is utilized for for speed ..

}


\usage{
BayesCpi(ga, numiter = 5000, Pi = .9, y)
}
\arguments{
\item{ga}{A matrix of genotypes with a number of rows identical to the number of genotyped individuals
          and a number of columns identical to the number of SNPs. Values in the matrix are 0, 1, 2, & 5
          for homozygous, heterozygous, other homozygous, & unknown genotypes, respectively.}
\item{numiter}{Number of iterations}
\item{Pi}{Proportion of SNP loci with 0 effect for Bayes C}
\item{y}{Trait phenotypes or conventional breeding values}
}


\details{
This function runs Bayes C and Cpi to estimate direct SNP effects and the
proportion of loci with nonzero effects based on a matrix of genotypes,
\code{ga} and a vector of adjusted phenotypes, \code{y}, (Habier et al., 2011; BMC Bioinformatics 12:186).
As in other bayesian alphabet, Bayes Cpi is essential in high dimensional data problems with highly
overparameterized models (\emph{p >> N}). It extends Bayes A and B to estimate the proportion of loci
with nonzero effect.

Other data management functions in \code{gdmp} can be used to construct the integer matrix of genotypes,
\code{ga}, to use as input to \code{BayesCpi}.

}


\value{

  A list object with a vector of SNP estimates \code{meanb} and a vector of genomic values for individuals,
  \code{aHat} are returned. It is also possible to extract the estimated number of SNP loci in \code{nLoci}.
}

\seealso{
\code{\link{BayesAB}}
}
\references{
  Habier et al. (2011).
  Extension of the bayesian alphabet for genomic selection.
  \emph{BMC Bioinformatics}, \bold{12}, 186.
}

