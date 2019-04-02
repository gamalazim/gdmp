#include "gt.h"

int isIdentical(int *x, int *y, int k, int n) {
// A function called from cGetHCS to compare two SNPs

  int i, xyNot = 0, xyNotOp = 0;

  for(i=0; i<n; i++) if(x[i] != y[i]) xyNot++;
  for(i=0; i<n; i++) if(x[i] != abs(y[i]-2)) xyNotOp++;
  if(xyNotOp < xyNot) xyNot = xyNotOp;
  if(xyNot <= k) return 1;
  else return 0;
}

void cGetHCS(int *data, int *Exclude, int *ll, int *nn, int *rAllow) {

  int i, j, n = *nn, l = *ll, k = 0, **gar;
  int allow = (*rAllow);

  // reversed dimensions of R's ga.r so you can send a SNP as a row
  gar = imatrix(0, n, 0, l);
  for(i=0; i<n; i++) for(j=0; j<l; j++) { gar[i][j] = data[k]; k++; }

  for(i=0; i<(n-1); i++) {
    if(Exclude[i]) {
      for(j=(i+1); j<n; j++) {
	     if(Exclude[j]) {
	       if(isIdentical(gar[i], gar[j], allow, l)) Exclude[j] = 0;
	     }
      }
    }
  }
  free_imatrix(gar, 0, n, 0, l);
}
