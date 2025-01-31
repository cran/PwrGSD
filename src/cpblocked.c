#include<R.h>
#include<Rmath.h>
#include "PwrGSD_mem.h"

#define MIN(x,y) (x < y ? x : y)

#include "PwrGSD.h"

typedef struct{
  int index;
  double time;
  int event;
  int arm;
} itea;

void cpblocked(itea *Yord, int *pn, double *time, int *nrisk, int *nevent, int *pntimes, int *pnevtypes, int *pnblocks)
{
  int n,i,j,k,l,nevty,nb,isev,ntimes,cont;
  int *dr,*dnev,*nr;
  double yhold,yn,yo;

  n = *pn;
  ntimes = *pntimes;
  nevty = *pnevtypes;
  nb = *pnblocks;
  f = &compitea;
  qsort(Yord, n, sizeof(itea), f);

  dr =   Calloc(nb, int);
  nr =   Calloc(nb, int);
  dnev = Calloc(nb*nevty, int);

  i=n-1;
  l=0;
  for(j=0;j<nb;j++) *(nr + j)=0;
  while((i>=0) && (l < ntimes)){
    isev = 0;
    for(j=0;j<nb;j++) *(dr + j) = 0;
    for(j=0;j<nb;j++) for(k=0;k<nevty;k++) *(dnev + nb*k + j) = 0;
    yn = (Yord+i)->time;
    yhold = yn;
    cont=1;
    while((yn == yhold) && cont){
      yo = yn;
      isev = isev || ((Yord+i)->event > 0);
      for(j=0;j<nb;j++)
      {
        for(k=0;k<nevty;k++)
          *(dnev + nb*k + j) = *(dnev + nb*k + j) + ((Yord+i)->arm==j) * ((Yord+i)->event==(k+1));
        *(dr + j) = *(dr + j) + ((Yord+i)->arm==j);
      }
      i--;
      cont=(i>=0);
      if(cont) yn=(Yord+i)->time;
    }
    for(j=0;j<nb;j++) *(nr + j) = *(nr + j) + *(dr+j);
    if(isev){
      for(j=0;j<nb;j++){
        *(nrisk + nb*(ntimes - 1 - l) + j) = *(nr + j);
        *(time + ntimes - 1 - l) = yo;
        for(k=0;k<nevty;k++)
          *(nevent + nevty*nb*(ntimes - 1 - l) + nevty*j + k) = *(dnev + k*nb + j);
      }
      l++;
    }
  }
  Free(dr);
  Free(dnev);
  Free(nr);
}

int compitea(const void *x, const void *y)
{
  itea *xx, *yy;
  xx = (itea *) x;
  yy = (itea *) y;
  return (1*(xx->time > yy->time) - 1*(xx->time < yy->time));
}
