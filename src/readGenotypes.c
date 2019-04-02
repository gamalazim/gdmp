// read genotypes.filled -- output of findhap

#include <R.h>
#include "gt.h"

void cReadGenotypes(int *r, int *c, int *h, int *g, char **path) {
//
// Used this because 'read.fwf' of R is extrtemely slow.
// The R function, 'read.fwf' could simpley take hours to read
// 'genotypes.filled' of a single chromosome.
//
// This was not easy to write because of the fixed format of the FORTRAN output
// of findahap. The only way was to create the 2-character string 'a' with '\0' as
// its 2nd character and line[j] as its 1st character, then use 'atoi(a)'. Note that
// atoi is declared as 'int atoi(const char *nptr)' which means that a string not a
// character is the argument of atoi().
//
// Note that 'line' must be allocated to avoid crashes if declared as line[size].
// Also 'line' is declared and allocated so that it acomodates any larger future
// widths caused by greater chips or even sequences of millions of loci.
//

  char *line; char a[2];
  int i, j, k, l;
  int rr = (*r);
  int cc = (*c)+31;
  int size = cc + 1024; // added to read the whole line (need to add 3 at least)
  FILE *fp;

  if((fp=fopen(path[0],"r"))==NULL)
    Error("cReadGenotypes(): Can't read file 'genotypes.filled'");

  line = cvector(0, size);
  a[1] = '\0';
  k = 0;
  l = 0;
  for(i=1; i<=rr; i++) {
    if(fgets(line, size, fp)) 
      sscanf(line, "%d %d %d", &h[l], &h[l+1], &h[l+2]);
    l+=3;
    for(j=31; j < cc; j++) { a[0] = line[j]; g[k++] = atoi(a); }
  }

  fclose(fp);
  free_cvector(line, 0, size);
}
