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

typedef void WtFun(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt);
WtFun flemhar, sflemhar, ramp, *wtfun;

void wlrstat(double *time, int *nrisk, int *nevent, double *wt, int *pntimes, double *UQ, 
	     double *varQ, double *UQt, double *varQt, double *var1t);

void WtdLogRank(double *TOS, int *Event, int *Arm, int *pn, int *wttyp, double *par,
		double *time, int *nrisk, int *nevent, double *wt, int *pntimes, double *UQ, double *varQ, 
                double *UQt, double *varQt, double *var1t)
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

  if(*wttyp==0) wtfun = &flemhar;
  if(*wttyp==1) wtfun = &sflemhar;
  if(*wttyp==2) wtfun = &ramp;
  (*wtfun)(time, nrisk, nevent, pntimes, par, wt);
  wlrstat(time, nrisk, nevent, wt, pntimes, UQ, varQ, UQt, varQt, var1t);

  Free(pnblocks);
  Free(YY);  
}

void wlrstat(double *time, int *nrisk, int *nevent, double *wt, int *pntimes, double *UQ, 
	     double *varQ, double *UQt, double *varQt, double *var1t)
{
  int nt,i;
  double xev, xev1, xri, xri1, Q, e, duQ, uQ, dv1, v1, dvQ, vQ;
  nt = *pntimes;
  duQ = 0.0;
  uQ = 0.0;
  dv1 = 0.0;
  v1 = 0.0;
  dvQ = 0.0;
  vQ = 0.0;
  for(i=0;i<nt;i++){
    xev = (double) (*(nevent + 2*i) + *(nevent + 2*i + 1));
    xev1 = (double) (*(nevent + 2*i + 1));
    xri = (double) (*(nrisk + 2*i) + *(nrisk + 2*i + 1));
    xri1 = (double) (*(nrisk + 2*i + 1));
    e = xri1/xri;
    Q = *(wt+i);
    duQ = Q * (xev1 - e * xev);
    uQ += duQ;
    *(UQt + i) = uQ;
    dv1 = e * (1.0 - e) * xev;
    v1 += dv1;
    *(var1t + i) = v1;
    dvQ = Q * Q * dv1;
    vQ += dvQ;
    *(varQt + i) = vQ;
  }
  *UQ = uQ;
  *varQ = vQ;
}

/*----------------------------------------------------------------------------------------
   Weighting Functions for Log-Rank Statistics:
----------------------------------------------------------------------------------------*/

/* 1.  Flemming Harington "G-rho,gamma" class:                                          */
void flemhar(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt)
{
  double s,xev,xri;
  int i,nt;
  s = 1.0;
  nt = *pntimes;
  for (i=0;i<nt;i++){
    xev = (double) (*(nevent + 2*i) + *(nevent + 2*i + 1));
    xri = (double) (*(nrisk + 2*i) + *(nrisk + 2*i + 1));
    s = s * (1.0-xev/xri);
    *(wt+i)=pow(s,*par)*pow(1.0-s,*(par+1));
  }
}

/* 2.  Flemming Harington "G-rho,gamma" class stopped at a specified time t:            */
void sflemhar(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt)
{
  double s,xev,xri,ti,val;
  int i,nt;
  s = 1.0;
  val = 1.0;
  nt = *pntimes;
  for (i=0;i<nt;i++){
    xev = (double) (*(nevent + 2*i) + *(nevent + 2*i + 1));
    xri = (double) (*(nrisk + 2*i) + *(nrisk + 2*i + 1));
    ti = *(time + i);
    s = s * (1.0-xev/xri);
    val = (ti <= *(par+2) ? *(wt+i)=pow(s,*par)*pow(1.0-s,*(par+1)) : val);
    *(wt+i)=val;
  }
}

/* 3. The ramp-plateau: increases linearly until a specified time, t, and then constant */
void ramp(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt)
{
  int i,nt;
  nt = *pntimes;
  for(i=0;i<nt;i++)
    *(wt + i) = MIN((*(time + i))/(*par), 1.0);
}


/*----------------------------------------------------------------------------------------
   add definition for alternate weighting function here...
   don't forget to make corresponding changes in this program and the
   R function interface to expand the menu of choices, and don't forget
   to add the function name to the globally declared type 'WtFun'...see above.
----------------------------------------------------------------------------------------*/
