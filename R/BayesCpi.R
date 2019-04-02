##
## BayesCpi
##
##
## R function based on Habier's Bayes C Pi
## The function utilizes the C Function cBaysCPi for speed ..
#'@export
BayesCpi <- function(ga, numiter=5000, Pi=.9, y) {

# Parameters
nmarkers    = dim(ga)[2];    # number of markers
logPi     = log(Pi);
logPiComp = log(1-Pi);

nrecords = dim(ga)[1];
x = cbind(1,ga);

storePi = array(0.0,numiter);

# inital values
vara = 1;
mean2pq = .5;
nua = 4;
varEffects  = vara/(nmarkers*(1-Pi)*mean2pq);
scalec      = varEffects*(nua-2)/nua;

b = array(0.0,nmarkers+1);
meanb   = b;
b[1]   = mean(y);
Var    = array(0.0,nmarkers);
ppa    = array(0.0,nmarkers);
piMean = 0.0;
meanVar = 0.0;
nLoci = 0;

# adjust y
ycorr = y - x%*%b;



##
##
## The C Part -- MCMC in C
##

out <- .C("cBaysCPi", as.integer(x), as.integer(numiter), as.integer(nrecords),
        as.integer(nmarkers), as.double(b), as.double(Var), as.double(ycorr),
        as.double(logPiComp), as.double(logPi), as.double(varEffects),
        as.double(scalec), as.integer(nua), as.double(vara), as.double(mean2pq),

        meanb=as.double(meanb),
        ppa=as.double(ppa),
        meanVar=as.double(meanVar),
        storePi=as.double(storePi),
        piMean=as.double(piMean),
        nLoci = as.integer(nLoci),
        PACKAGE="gdmp")



cat("No. Loci in model = ", out$nLoci, "\n");
list(

  meanb   = out$meanb/numiter,
  ppa     = out$ppa/numiter,
  piMean  = out$piMean/numiter,
  meanVar = out$meanVar/numiter,
  storePi = out$storePi,
  aHat    = x %*% (out$meanb/numiter),
  nLoci   = out$nLoci
  )
}

