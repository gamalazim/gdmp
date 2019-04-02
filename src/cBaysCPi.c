///
/// R CMD SHLIB -o cBaysCPi.so -lm cBaysCPi.c gt.c ranlib.c com.c

#include "gt.h"
#include "ranlib.h"

// mcmc sampling


void cBaysCPi(int *data, int *pnumiter, int *pnrecords, int *pnmarkers,
          double *b, double *Var, double *ycorr, double *plogPiComp,
          double *plogPi, double *pvarEffects, double *pscalec, int *pnua,
          double *pvara, double *pmean2pq,
          double *meanb, double *ppa, double *pmeanVar, double *storePi,
          double *ppiMean, int *pnLoci ) {

  int i, j, k, locus, countLoci, iter, **x, nLoci=*pnLoci;
  int numiter= *pnumiter, nrecords= *pnrecords, nmarkers= *pnmarkers, nua= *pnua;
  double piMean = *ppiMean, logPiComp = *plogPiComp, logPi = *plogPi, vara= *pvara;
  double varEffects= *pvarEffects, scalec= *pscalec, meanVar= *pmeanVar, mean2pq= *pmean2pq;
  double vare, Sum, v0, v1, logDelta0, logDelta1, probDelta1, u;
  double lhs, rhs, invLhs, mean, xpx, aa, bb, Pi;

  // Stuff needed for sampling from statistical distributions
  long is1, is2;
	char phrase[1024];
	time_t rawtime;
	struct tm * timeinfo;

  // Use current date and time to set seeds
	time( &rawtime );
	timeinfo = localtime( &rawtime );
	strcpy(phrase, asctime (timeinfo) );
  phrtsd(phrase, &is1, &is2);
  setall(is1, is2);

  x = imatrix(0, nrecords-1, 0, nmarkers); // nrecords X (nmarkers + 1)
  k = 0;
  for(j=0; j<=nmarkers; j++)
    for(i=0; i<nrecords; i++) {
      x[i][j] = data[k]; k++; }


  for (iter=1; iter<=numiter; iter++){

  // sample vare
	// vare = ( t(ycorr)%*%ycorr )/rchisq(1,nrecords + 3);
	Sum = 0;
  for(i=0; i<nrecords; i++) Sum += (ycorr[i]*ycorr[i]);
  // vare = Sum/rchisq(1,nrecords + 3);
  vare = gennch((nrecords+3.0), 0.0);
	vare = Sum/vare;

  // sample intercept
  Sum = 0;
 	//ycorr = ycorr + x[,1]*b[1];
	for(i=0; i<nrecords; i++) { ycorr[i] += b[0]; Sum += ycorr[i]; }
	//rhs    = sum(ycorr)/vare;
	rhs    = Sum/vare;

  //invLhs = 1.0/(nrecords/vare);
	invLhs = vare/nrecords;
	mean = rhs*invLhs;

	// b[1] = rnorm(1,mean,sqrt(invLhs));
	b[0] = gennor( mean, sqrt(invLhs) );

  //ycorr = ycorr - x[,1]*b[1]
	for(i=0; i<nrecords; i++) ycorr[i] -= b[0];

	meanb[0] = meanb[0] + b[0];

  // sample delta (slide 48)  and effect for each locus
  nLoci = 0;
	for (locus=1; locus<=nmarkers; locus++){ // jump over Mu
    rhs = xpx = 0.0;
    //ycorr = ycorr + x[,1+locus]*b[1+locus];
		for(i=0; i<nrecords; i++) {
      ycorr[i] += (x[i][locus]*b[locus]);
      rhs += (x[i][locus]*ycorr[i]);
      xpx += (x[i][locus]*x[i][locus]);
    }
		//rhs = t(x[,1+locus])%*%ycorr;
		//xpx = t(x[,1+locus])%*%x[,1+locus];

    v0  =  xpx*vare;
		v1  =  (square(xpx)*varEffects + xpx*vare);
		logDelta0 = -0.5*(log(v0) + square(rhs)/v0) + logPi;
		logDelta1 = -0.5*(log(v1) + square(rhs)/v1) + logPiComp;
		probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1));
		u = ranf();  // unif(0, 1)
		if(u < probDelta1) {
			nLoci = nLoci + 1;
			lhs = xpx/vare + 1.0/varEffects;
			invLhs = 1.0/lhs;
			mean = invLhs*rhs/vare;
			b[locus]= gennor( mean, sqrt(invLhs) ); // b[locus+1]= rnorm(1,mean,sqrt(invLhs));

      // ycorr = ycorr - x[,1+locus]*b[1+locus];
      for(i=0; i<nrecords; i++) ycorr[i] -= (x[i][locus]*b[locus]);

			meanb[locus] = meanb[locus] + b[locus];
			ppa[locus-1] = ppa[locus-1] + 1;
			Var[locus-1] = varEffects;
		}
		else {
			b[locus] = 0.0;
			Var[locus-1] = 0.0;
		}

	}

  // sample common variance
	countLoci = 0;
	Sum = 0.0;
	for(locus=0; locus<nmarkers; locus++){
		if(Var[locus]>0.0){
			countLoci++;
			Sum += square(b[1+locus]);
		}
	}
	// varEffects = (scalec*nua + Sum)/rchisq(1,nua+countLoci); # slide 50
  varEffects = (scalec*nua + Sum);
  varEffects /= gennch( (nua+countLoci), 0.0 );
	meanVar = meanVar + varEffects;

  // sample Pi
	aa = nmarkers-countLoci + 1;
	bb = countLoci + 1;
	//Pi = rbeta(1, aa, bb);  # slide 52
	Pi = genbet(aa, bb);
	storePi[iter-1] = Pi;
	piMean = piMean + Pi;
	scalec =  (nua-2)/nua * vara/((1-Pi)*nmarkers*mean2pq);
	logPi     = log(Pi);
	logPiComp = log(1-Pi);
  }

  *pnLoci = nLoci;
  *ppiMean = piMean;
  *pmeanVar= meanVar;

  free_imatrix(x, 0, nrecords-1, 0, nmarkers);

}

