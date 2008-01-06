#include<R.h>
#include<Rmath.h>
#define MIN(x,y) (x < y ? x : y)

typedef struct{
  int index;
  double time;
  int event;
  int arm;
} itea;

void cpblocked(itea *Yord, int *pn, double *time, int *nrisk, int *nevent, int *pntimes, int *pnblocks);

void mysurvfit(double *TOS, int *Event, int *Arm, int *pn, double *time, int *nrisk, int *nevent, 
               int *pntimes, int *pnblocks)
{
  int i,j,n;
  itea *YY;

  n = *pn;
  
  YY = (itea *) Calloc(n, itea);

  for (i=0;i<n;i++){
    (YY+i)->index = i;
    (YY+i)->time = *(TOS+i);
    (YY+i)->event = *(Event+i);
    (YY+i)->arm = *(Arm+i);
  }

  cpblocked(YY, pn, time, nrisk, nevent, pntimes, pnblocks);

  Free(YY);
}
