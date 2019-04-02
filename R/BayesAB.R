##
## BayesAB
##
##
## R function based on RL Fernando's Bayes A and B
## BayesAB utilizes the C Function cBaysABF for speed ..
##

#'@export
BayesAB <- function(ga, numiter=5000, numMHIter=10, Pi=.9, y) {

# Parameters
nmarkers    = dim(ga)[2];    # number of markers
logPi     = log(Pi);
logPiComp = log(1-Pi);

nrecords = dim(ga)[1];
x = cbind(1,ga);

# inital values
vara=1;
mean2pq=.5;
scaleb   =  0.5*vara/(nmarkers*(1-Pi)*mean2pq);
consta   = 2.0*log(scaleb*2.0);

b = array(0.0,nmarkers+1);
meanb   = b;
b[1]   = mean(y);
Var    = array(0.0,nmarkers);
ppa    = array(0.0,nmarkers);
nLoci  = 0;

# adjust y
ycorr = y - x%*%b;

##
##
## The C Part -- MCMC in C
##

out <- .C("cBaysABF", as.integer(x), as.integer(numiter), as.integer(numMHIter),
        as.integer(nrecords),
        as.integer(nmarkers), as.double(b), as.double(Var), as.double(ycorr),
        as.double(logPiComp), as.double(logPi), as.double(scaleb),

        meanb=as.double(meanb),
        ppa=as.double(ppa),
        nLoci = as.integer(nLoci),
        PACKAGE="gdmp")


  list(
  meanb   = out$meanb/numiter,
  ppa     = out$ppa/numiter,
  aHat    = x %*% (out$meanb/numiter),
  nLoci   = out$nLoci
  )
}
