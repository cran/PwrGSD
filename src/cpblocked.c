#include<R.h>
#include<Rmath.h>
#define MIN(x,y) (x < y ? x : y)

typedef struct{
  int index;
  double time;
  int event;
  int arm;
} itea;

typedef int CmprFun(const void *x, const void *y);
CmprFun compitea, *f;

void cpblocked(itea *Yord, int *pn, double *time, int *nrisk, int *nevent, int *pntimes, int *pnblocks)
{
  int n,i,j,l,nb,isev,ntimes,cont;
  int *dr,*dnev,*nr;
  double yhold,yn,yo;

  n = *pn;
  ntimes = *pntimes;
  nb = *pnblocks;
  f = &compitea;
  qsort(Yord, n, sizeof(itea), f);

  dr = (int *)Calloc(nb, int);
  dnev = (int *)Calloc(nb, int);
  nr = (int *)Calloc(nb, int);

  i=n-1;
  l=0;
  for(j=0;j<nb;j++) *(nr + j)=0;
  while(i>=0 && l < ntimes){
    isev = 0;
    for(j=0;j<nb;j++) *(dr + j) = 0;
    for(j=0;j<nb;j++) *(dnev +j) = 0;
    yn = (Yord+i)->time;
    yhold = yn;
    cont=1;
    while(yn == yhold && cont){
      yo = yn;
      isev = isev || 1*((Yord+i)->event > 0);
      for(j=0;j<nb;j++){
        *(dnev + j) = *(dnev + j) + ((Yord+i)->arm==j) * (Yord+i)->event;
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
        *(nevent + nb*(ntimes - 1 - l) + j) = *(dnev + j);
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
