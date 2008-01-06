#include<R.h>
#include<Rmath.h>
#define MIN(x,y) (x < y ? x : y)
#define EPS 1.0e-10

typedef struct{
  int index;
  double time;
  int event;
  int arm;
} itea;

void cpblocked(itea *Yord, int *pn, double *time, int *nrisk, int *nevent, int *pntimes, int *pnblocks);

typedef void WtFun(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt);
WtFun flemhar, sflemhar, ramp, *wtfun;

void ISDstat(double *time, int *nrisk, int *nevent, int *pntimes, double *wt, double *stat, 
	     double *var);

void IntSurvDiff(double *TOS, int *Event, int *Arm, int *pn, int *wttyp, double *par, 
		double *time, int *nrisk, int *nevent, int *pntimes, double *stat, 
		double *var, double *wt)
{
  int i,n,ntimes;
  int *pnblocks;
  pnblocks = (int *)Calloc(1, int);
  *pnblocks = 2;

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

  if(*wttyp==0) wtfun = *flemhar;
  if(*wttyp==1) wtfun = *sflemhar;
  if(*wttyp==2) wtfun = *ramp;
  (*wtfun)(time, nrisk, nevent, pntimes, par, wt);
  ISDstat(time, nrisk, nevent, pntimes, wt, stat, var);

  Free(pnblocks);
  Free(YY);  
}

void ISDstat(double *time, int *nrisk, int *nevent, int *pntimes, double *wt, double *stat, 
	     double *var)
{
  int nt,ntm1,i;
  double xev0, xev1, xri0, xri1, q, stat_, var_, S0_, S1_, t_im1, dt, cs0, cs1, V;
  double *dV0, *dV1, *qS0dt, *qS1dt;
  nt = *pntimes;

  dV0=Calloc(nt,double);
  dV1=Calloc(nt,double);
  qS0dt=Calloc(nt,double);
  qS1dt=Calloc(nt,double);

  ntm1 = nt-1;
  stat_ = 0.0;
  var_ = 0.0;
  S0_=1.0;
  S1_=1.0;
  t_im1 = 0.0;
  for(i=0;i<nt;i++){
    q = *(wt+i);
    xev0 = (double) *(nevent + 2*i);
    xev1 = (double) *(nevent + 2*i + 1);
    xri0 = (double) *(nrisk + 2*i);
    xri1 = (double) *(nrisk + 2*i + 1);
    S0_ *= (xri0 > EPS ? (1.0 - xev0/xri0) : 1.0);
    S1_ *= (xri1 > EPS ? (1.0 - xev1/xri1) : 1.0);
    dt = *(time+i) - t_im1;
    t_im1 = *(time+i);
    *(qS0dt+i) = q * S0_ * dt;
    *(qS1dt+i) = q * S1_ * dt;
    stat_ += *(qS1dt+i) - *(qS0dt+i);
    *(dV0+i) = (xri0 > EPS ? xev0/(xri0*xri0) : 0.0);
    *(dV1+i) = (xri1 > EPS ? xev1/(xri1*xri1) : 0.0);
  }
  *stat = stat_;
  V = cs0 = cs1 = 0.0;
  for(i=ntm1;i>=0;i--){
    cs0 += *(qS0dt+i);
    cs1 += *(qS1dt+i);
    V += cs0 * cs0 * (*(dV0 + i)) + cs1 * cs1 * (*(dV1 + i));
  }
  *var = V;
  Free(dV0); 
  Free(dV1);
  Free(qS0dt); 
  Free(qS1dt);
}
