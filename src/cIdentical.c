
#include <R.h>
void cIdentical(int *x, int *y, int *k, int *n, int *result) {
  int nn = (*n);
  int kk = (*k);
  int i, xyNot = 0, xyNotOp = 0;

  for(i=0; i<nn; i++) if(x[i] != y[i]) xyNot++;
  for(i=0; i<nn; i++) if(x[i] != abs(y[i]-2)) xyNotOp++;
  if(xyNotOp < xyNot) xyNot = xyNotOp;
  if(xyNot <= kk) (*result) = 1;
  else (*result) = 0;
}

