#include<R.h>

void mystack(int *pn, int *pnfus, int *pnfuvars, int *pnbasevars, double *basevars, double *fuvars, double *out)
{
    int i, j, k, l, n, nfus, nfuvars, nbasevars;
    
    n = *pn;
    nfus = *pnfus;
    nfuvars = *pnfuvars;
    nbasevars = *pnbasevars;

    for(i=0;i<n;i++){
      for(j=0;j<nfus;j++){
        for(k=0;k<nbasevars;k++) *(out + n*nfus*k + nfus*i + j) = *(basevars + n*k + i);
        *(out + n*nfus*nbasevars + nfus*i + j) = j;
	for(k=0;k<nfuvars;k++) *(out + n*nfus*(nbasevars+1+k) + nfus*i + j) = *(fuvars + nfus*n*k + n*j + i);
      }
    }
}
