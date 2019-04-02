#include "gt.h"
#include "ranlib.h"

// MCMC Sampling

void cBaysABF(int *data, int *pnumiter, int *pnumMHIter, int *pnrecords, int *pnmarkers,
  double *b, double *Var, double *ycorr, double *plogPiComp, double *plogPi,
  double *pscaleb, double *meanb, double *ppa, int *pnLoci ) {

  int i, j, k, locus, iter, mhiter, **x, nLoci=*pnLoci, numMHIter = *pnumMHIter;
  int numiter= *pnumiter, nrecords= *pnrecords, nmarkers= *pnmarkers;
  double logPiComp = *plogPiComp, logPi = *plogPi;
  double scaleb= *pscaleb;
  double vare, Sum, v1, v2, u;
  double lhs, rhs, invLhs, mean, xpx;

  double consta   = 2.0*log(scaleb*2.0);

  double logDataNullModel, logDataOld, logDataNew;
  double logPriorOld, logPriorNew, logPosteriorOld, logPosteriorNew;
  double varCandidate, logProposalNew, logProposalOld, sv, acceptProb;

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

  // sample variance and effect for each locus
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

    v1  =  (square(xpx)*Var[locus-1] + xpx*vare);
		v2  =  xpx*vare;

		logDataNullModel = -0.5*(log(v2) + square(rhs)/v2);
		if (Var[locus-1] > 0.0){
			logDataOld  = -0.5*(log(v1) + square(rhs)/v1);
			sv = scaleb;
			logPriorOld = consta -3.0*log(Var[locus-1]) - sv*4/(2.0*Var[locus-1]) + logPiComp;
		}
		else {
			logDataOld  = logDataNullModel;
			logPriorOld = logPi;
		}
		logPosteriorOld = logDataOld + logPriorOld;
		for (mhiter=1; mhiter<=numMHIter; mhiter++){
			u = ranf(); // runif(1);
			varCandidate = 0;
			if (u > 0.5){ // include in the model and check if it was there before
				if (Var[locus-1] > 0){
					varCandidate = Var[locus-1]*2;
					varCandidate /= gennch(4, 0);
				}
				else {
					varCandidate = scaleb*4;
					varCandidate /= gennch(4, 0);
				}
			}
			if (varCandidate > 0.0){
				v1  =  (square(xpx)*varCandidate + xpx*vare);
				logDataNew =  -0.5*(log(v1) + square(rhs)/v1);
				sv = scaleb;
				logPriorNew = consta -3.0*log(varCandidate) - sv*4/(2.0*varCandidate) + logPiComp;
				logPosteriorNew = logDataNew + logPriorNew;
				if(Var[locus-1]>0){
					sv = varCandidate*0.5;
					logProposalOld = 2.0*log(sv*2.0) -3.0*log(Var[locus-1]) - sv*4/(2.0*Var[locus-1]);
					sv = Var[locus-1]*0.5;
					logProposalNew = 2.0*log(sv*2.0) -3.0*log(varCandidate) - sv*4/(2.0*varCandidate);
				}
				else {
					logProposalOld = 0.0;
					sv = scaleb;
					logProposalNew = consta -3.0*log(varCandidate) - sv*4/(2.0*varCandidate);
				}
			}
			else {
				logDataNew = logDataNullModel;
				logPriorNew = logPi;
				logPosteriorNew = logDataNew + logPriorNew;
				if (Var[locus-1]>0){
					sv = scaleb;
					logProposalOld = consta -3.0*log(Var[locus-1]) - sv*4/(2.0*Var[locus-1]);
					logProposalNew = 0.0;
				}
				else {
					logProposalOld = 0.0;
					logProposalNew = 0.0;
				}
			}
			acceptProb = exp(logPosteriorNew+logProposalOld-logPosteriorOld-logProposalNew);
			u = ranf(); // runif(1);
// printf("locus= %d, Var[%d]= %f, varCandidate= %f, u= %f, acceptProb= %f \n", locus, locus,Var[locus-1],varCandidate,u,acceptProb);
			if(u < acceptProb) {
				Var[locus-1] = varCandidate;
				logPosteriorOld = logPosteriorNew;
			}

		}
		if(Var[locus-1]) {
			nLoci = nLoci + 1;
			lhs = xpx/vare + 1.0/Var[locus-1];
			invLhs = 1.0/lhs;
			mean = invLhs*rhs/vare;

      //b[1+locus]= rnorm(1,mean,sqrt(invLhs));
			b[locus]= gennor( mean, sqrt(invLhs) );

      //ycorr = ycorr - x[,1+locus]*b[1+locus];
      for(i=0; i<nrecords; i++) ycorr[i] -= (x[i][locus]*b[locus]);

			meanb[locus] = meanb[locus] + b[locus];
			ppa[locus-1] = ppa[locus-1] + 1;
		}
		else b[locus] = 0.0;
	}

  }

  *pnLoci = nLoci;

  free_imatrix(x, 0, nrecords-1, 0, nmarkers);

}

