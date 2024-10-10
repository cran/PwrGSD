#include<R.h>
#include<Rmath.h>
#include "PwrGSD_mem.h"

#define MIN(x,y) (x<y ? x : y)
#define MAX(x,y) (x>y ? x : y)
#define SQ(x) ((x)*(x))

typedef struct{
  int index;
  double time;
  int event;
  int arm;
} itea;

typedef struct{
  int index;
  double value;
} IND_VAL;

typedef int CmprFun(const void *x, const void *y);
CmprFun cmpr_IND_VAL;

void printmat(double *a, int nr, int nc, char *nm);
// void print_itea(itea *YY, int n, char *label);
// void print_tre(double *time, int *risk, int *event, int n, char *label);
// void print_pair(double *x, double *y, int n, char *x_lab, char *y_lab);
void makeRR(double *TI, double *Toth, double *T_R, double *tlook, int *pntot, double *pRR);
void makeYY(double *TI, double *Toth, double *u, double *tlook, int *pnn, int *pntimes, itea *YY);
void cpblocked(itea *Yord, int *pn, double *time, int *nrisk, int *nevent, int *pntimes, 
               int *pnevtypes, int *pnblocks);
void myapprox(double *x, double *y, int nxy, double *xout, double *yout, int nout, int method, 
              double yleft, double yright, double f);
void invrt3by3(double *A, double *Ainv);
void Hproj(double *ti, double *H, int *pnti, double *coef2, double *t_proj, double *H_proj, int *pn_proj);
void ls_quad(double *ti, double *H, int *pnti, double *coef3);


void project_end(double *T_R, double *T0, double *T1, double *Tc0, double *Tc1, int *pn, double *t_proj, 
                 double *v_Tend_proj, double *m_Tend_proj)
{
  int n,ntot,ntimes,one=1,three=3,i; 
  int *int_buff,*pntot,*pntimes,*nrisk,*nevent;
  double s,t_Q,tau_m_ter,H_tQ,H_tau_m_ter,H_tau,H_min,H_max,RR,norm,I_1,I_2,I_3,I_4,J_1,J_2;
  double *dbl_buff,*tlook,*time,*TI,*Toth,*H,*coef,*H_proj,*pRR;
  itea *YY;
  /*
        *t_proj=*par;              a.k.a. t_Q, the time the wt function levels off
        *(t_proj+1)=t_end-*accru;  a.k.a. T_end - t_er, the time of live removal
        *(t_proj+2)=t_end;         time of trial close-out.
  */
  n        = *pn;
  ntot     = 2*n;
  int_buff = Calloc(2*ntot+2, int);
  dbl_buff = Calloc(ntot+7, double);
  TI       = Calloc(ntot, double);
  Toth     = Calloc(ntot, double);

  pntot    = int_buff;
  pntimes  = int_buff + 1;
  nrisk    = int_buff + 2;
  nevent   = int_buff + 2 + ntot;
  tlook    = dbl_buff;
  coef     = dbl_buff + 1;
  H_proj   = dbl_buff + 3;
  pRR      = dbl_buff + 6;
  time     = dbl_buff + 7;

  YY       = Calloc(ntot, itea);

  *tlook   = *(t_proj+1);
  *pntot   = ntot;

  for(i=0;i<n;i++)
  {
    *(TI   + i)     = *(T0  + i);
    *(TI   + n + i) = *(T1  + i);
    *(Toth + i)     = *(Tc0 + i);
    *(Toth + n + i) = *(Tc1 + i);
  }

  makeYY(TI,Toth,T_R,tlook,pntot,pntimes,YY);
  cpblocked(YY, pntot, time, nrisk, nevent, pntimes, &one, &one);
  ntimes=*pntimes;

  H = Calloc(ntimes, double);
  s=0.0;
  for(i=0;i<ntimes;i++)
  {
    s += ((double)*(nevent + i))/((double)(*(nrisk + i)));
    *(H + i) = s;
  }

  ls_quad(time, H, pntimes, coef);
  Hproj(time, H, pntimes, coef+1, t_proj, H_proj, &three);

  makeRR(TI,Toth,T_R,tlook,pntot,pRR);
  RR = *pRR;

  t_Q         = *t_proj;
  tau_m_ter   = *(t_proj + 1);
  H_tQ        = *H_proj;
  H_tau_m_ter = *(H_proj + 1);
  H_tau       = *(H_proj + 2);

  norm = 1.0;
  if(H_tQ>0) norm = 1.0-exp(-H_tQ);
  /*
   compute the end of trial variance
   under Stopped Fleming-Harrington (0, 1) weighting (we assume that the
   actual weighting function which is deterministic linear rise from 0 to 1
   between 0 < t < t_c  is approximately the same as Stopped Fleming-Harrington weighting.
  */
  H_min = MIN(H_tQ, H_tau_m_ter);
  I_1   =  0.25*( (1.0-exp(-(RR+1.0)*H_min))/(RR+1.0) - 2.0*(1.0-exp(-(RR+2.0)*H_min))/(RR+2.0) + 
		  (1.0-exp(-(RR+3.0)*H_min))/(RR+3.0) );


  /* t_k is tau, t_min is tau_m_ter    t_min <- pmin(t.k, tau.m.ter) */
  /* H_min is H_tau_m_ter              H_min <- pmin(H.tk, H.tau.m.ter) */
  /* t_c is t_Q */
  I_2 = ((double)(t_Q < tau_m_ter))*SQ(1.0-exp(-H_tQ)) * (exp(-(RR+1.0)*H_tQ) - exp(-(RR+1.0)*H_tau_m_ter))/(4.0*(RR+1.0));

  /* t_k is tau, t_c is t_Q, t_min is t_Q   t.min <- pmin(t.k, t.c) */
  /* H_min is H_tQ */
  I_3 = ((double)(tau_m_ter < t_Q))/(4.0*(H_tau-H_tau_m_ter)) *
    
           ((exp(-(RR+1.0)*H_tau_m_ter)/(RR+1.0)- 2.0* exp(-(RR+2.0)*H_tau_m_ter)/(RR+2.0) +
             exp(-(RR+3.0)*H_tau_m_ter)/(RR+3.0)) * (H_tau - H_tau_m_ter) -

            (exp(-(RR+1.0)*H_tQ)/(RR+1.0) - 2*exp(-(RR+2.0)*H_tQ)/(RR+2.0) +
             exp(-(RR+3.0)*H_tQ)/(RR+3.0)) *(H_tau - H_tQ) - 

           ((exp(-(RR+1.0)*H_tau_m_ter)-exp(-(RR+1.0)*H_tQ))/SQ(RR+1.0) - 
          2.0*(exp(-(RR+2.0)*H_tau_m_ter)-exp(-(RR+2.0)*H_tQ))/SQ(RR+2.0) +

	    (exp(-(RR+3.0)*H_tau_m_ter)-exp(-(RR+3.0)*H_tQ))/SQ(RR+3.0)));
  
  /* t_k is tau, which is greater than t_max */
  H_max = MAX(H_tau_m_ter, H_tQ);
  I_4 = SQ(1.0-exp(-H_tQ))/(4.0*(H_tau-H_tau_m_ter)) * 
           (exp(-(RR+1.0)*H_max)/(RR+1.0) * (H_tau -H_max) - (exp(-(RR+1.0)*H_max)-exp(-(RR+1.0)*H_tau))/SQ(RR+1.0));

  *v_Tend_proj = (I_1 + I_2 + I_3 + I_4)/SQ(norm);
  
  /*
   compute drift at analysis j=1,2,...,k
  
   E U(t_j)/V(\tau)^{1/2} = n^{1/2} \beta_{\star} m(t_j, q_0)/m(\tau, q_0) m(\tau, 1)/V(\tau)^{1/2}
  
   using q_0(\xi) = Q(\xi) = 1 - exp(-H(\xi) /\ H(t_c)) 
  
   don't forget to define H.tk, H.tc = ???
  */
  J_1 = 0.25 * ((1.0-exp(-(RR+1.0)*H_tau_m_ter))/(RR+1.0) - (1.0-exp(-(RR+2.0)*H_tau_m_ter))/(RR+2.0));

  J_2 = 0.25 * (exp(-(RR+1.0)*H_tau_m_ter)/(RR+1.0) - exp(-(RR+2.0)*H_tau_m_ter)/(RR+2.0)) - 
               ((exp(-(RR+1.0)*H_tau_m_ter) - exp(-(RR+1.0)*H_tau))/SQ(RR+1.0) - 
                (exp(-(RR+2.0)*H_tau_m_ter) - exp(-(RR+2.0)*H_tau))/SQ(RR+2.0))/(4.0*(H_tau -H_tau_m_ter));

  *m_Tend_proj = (J_1 + J_2)/norm;

  Free(int_buff);
  Free(dbl_buff);
  Free(TI);
  Free(Toth);
  Free(YY);
  Free(H);
}

void makeYY(double *TI, double *Toth, double *u, double *tlook, int *pnn, int *pntimes, itea *YY)
{
  int ndths,j,i,nn=*pnn;
  double tout,tmpp;

  ndths = 0;
  j=0;
  for(i=0;i<nn;i++){
    tout = (*tlook - *(u+i) > 0.0 ? *tlook - *(u+i) : 0.0);
    if(*(u+i)<*tlook){
      tmpp = (*(Toth+i) < tout ? *(Toth+i) : tout);
      (YY+j)->index = j;
      (YY+j)->time = (*(TI+i) <= tmpp ? *(TI+i) : tmpp);
      (YY+j)->event = (*(TI+i) <= tmpp ? 1 : 0);
      ndths += (YY+j)->event;
      (YY+j)->arm = 0;
      j++;
    }
  }
  *pntimes = ndths;
}

void  makeRR(double *TI, double *Toth, double *T_R, double *tlook, int *pntot, double *pRR)
{
  int i,ntot=*pntot;
  double NCaM, NOthM;

  NCaM=0.0;
  NOthM=0.0;
  for(i=0;i<ntot;i++)
  {
    NCaM  += (double) ((*(TI   + i)) < MIN((*(Toth + i)),(*tlook) - (*(T_R + i))));
    NOthM += (double) ((*(Toth + i)) < MIN((*(TI   + i)),(*tlook) - (*(T_R + i))));
  }
  *pRR = NOthM/NCaM;
}

void Hproj(double *ti, double *H, int *pnti, double *coef2, double *t_proj, double *H_proj, int *pn_proj)
{
  int n_proj, nti, i, l, n_less;
  double max_ti, max_H;
  double *t_proj_s, *H_proj_s;
  IND_VAL *ind_val;
  CmprFun *f;

  n_proj = *pn_proj;
  nti = *pnti;

  t_proj_s = Calloc(n_proj, double);
  H_proj_s = Calloc(n_proj, double);
  ind_val = Calloc(n_proj, IND_VAL);

  max_ti = *(ti + nti -1);
  max_H = *(H + nti -1);

  for(i=0;i<n_proj;i++)
  {
    (ind_val + i)->index = i;
    (ind_val + i)->value= *(t_proj + i);
  }  
  
  f = &cmpr_IND_VAL;

  qsort(ind_val, n_proj, sizeof(IND_VAL), f);

  l=0;
  for(i=0;i<n_proj;i++)
  {
    *(t_proj_s + i) = (ind_val + i)->value;
    if(*(t_proj_s + i) <= max_ti) l++;
  }
  n_less=l;
  if(n_less > 0) myapprox(ti,H,nti,t_proj_s,H_proj_s,n_less,2,0,max_H,0);
  for(i=n_less;i<n_proj;i++) *(H_proj_s + i) = max_H + (*coef2 + 2*(*(coef2+1))*max_ti)*(*(t_proj + i) - max_ti);

  for(i=0;i<n_proj;i++) *(H_proj + (ind_val+i)->index) = *(H_proj_s + i);

  Free(ind_val);
  Free(H_proj_s);
  Free(t_proj_s);
}


void ls_quad(double *ti, double *H, int *pnti, double *coef3)
{
  int i,k,nti=*pnti,p=3;
  double s_ti, s_ti2, s_ti3, s_ti4, s_H, s_tiH, s_ti2H, s;
  double *XX, *XX_i, *XH;

  XX   = Calloc(p*p,double);
  XX_i = Calloc(p*p,double);
  XH   = Calloc(p,double);

  s_ti=0.0;
  s_ti2=0.0;
  s_ti3=0.0;
  s_ti4=0.0;
  s_H = 0.0;
  s_tiH = 0.0;
  s_ti2H = 0.0;
  for(i=0;i<nti;i++)
  {
    s_ti  += (*(ti + i));
    s_ti2 += (*(ti + i)) * (*(ti+i));
    s_ti3 += (*(ti + i)) * (*(ti+i)) * (*(ti+i));
    s_ti4 += (*(ti + i)) * (*(ti+i)) * (*(ti+i)) * (*(ti+i));

    s_H    += (*(H + i));
    s_tiH  += (*(ti + i)) * (*(H + i));
    s_ti2H += (*(ti + i)) * (*(ti + i)) * (*(H + i));
  }

  *XX = (double)nti;
  *(XX + 1) = *(XX + 3) = s_ti;
  *(XX + 2) = *(XX + 4) = *(XX + 6) = s_ti2;
  *(XX + 5) = *(XX + 7) = s_ti3;
  *(XX + 8) = s_ti4;

  *XH = s_H;
  *(XH + 1) = s_tiH;
  *(XH + 2) = s_ti2H;

  invrt3by3(XX, XX_i);

  for(i=0;i<p;i++)
  {
    s=0.0;
    for(k=0;k<p;k++) s += (*(XX_i + p*k + i))*(*(XH + k));
    *(coef3 + i) = s;
  }
}

void invrt3by3(double *A, double *Ainv)
{
  int i,j,k,p=3;
  double s,a,b,c,d,e,f;
  double *buff, *sqrA;
  buff = Calloc(p*p, double);
  sqrA = Calloc(p*p,double);

  *sqrA = pow(*A, 0.5);

  j=0;
  while((j>=0) && (j < p))
  {
    s = 0.0;
    k=0;
    while((k>=0) && (k <j))
    {
      s += (*(sqrA + p*k + j)) * (*(sqrA + p*k + j));
      k++;
    }
    *(sqrA + p*j + j) = pow(*(A + p*j + j) - s, 0.5);
    i=j+1;
    while((i>j) && (i<p))
    {
      s=0.0;
      k=0;
      while((k>=0) && (k <j))
      {
        s += (*(sqrA + p*k + i))*(*(sqrA + p*k + j));
        k++;
      }
      *(sqrA + p*j + i) = ((*(A + j*p + i)) - s)/(*(sqrA + p*j + j));
      i++;
    }
    j++;
  }

  a=*sqrA;
  b=*(sqrA+1);
  c=*(sqrA+2);
  d=*(sqrA+4);
  e=*(sqrA+5);
  f=*(sqrA+8);

  *buff = 1.0/a;
  *(buff+1)=-b/(a*d);
  *(buff+2)=(b*e-c*d)/(a*d*f);
  *(buff+4)=1.0/d;
  *(buff+5)=-e/(d*f);
  *(buff+8)=1.0/f;

  for(i=0;i<p;i++)
  for(j=0;j<p;j++)
  {
    s=0.0;
    for(k=0;k<p;k++) s+= (*(buff + p*i + k)) * (*(buff + p*j + k));
    *(Ainv + p*j + i) = s;
  }

  /* 
  printmat(A, p, p, "A");
  printmat(sqrA, p, p, "sqrA");
  printmat(Ainv, p, p, "Ainv");
  */
  
  Free(sqrA);
}

// void print_itea(itea *YY, int n, char *label)
// {
//   int i;
//   char cc;
//   Rprintf("%s:\n",label);
//   for(i=0;i<n;i++)
//   {
//     if((n%31)) Rprintf("index\t time\t event\t arm\n");
//     if((n%32)) scanf("continue? %c?\n",&cc);
//     Rprintf("%d\t %g\t %d\t %d\n", (YY+i)->index,(YY+i)->time,(YY+i)->event,(YY+i)->arm);
//   }
// }

// void print_tre(double *time, int *risk, int *event, int n, char *label)
// {
//   int i;
//   char cc;
//   Rprintf("%s:\n",label);
//   for(i=0;i<n;i++)
//   {
//     if((n%31)) Rprintf("time\t nrisk\t nevent\n");
//     if((n%32)) scanf("continue? %c?\n",&cc);
//     Rprintf("%g\t %d\t %d\n", *(time + i),*(risk + i),*(event + i));
//   }
// }

// void print_pair(double *x, double *y, int n, char *x_lab, char *y_lab)
// {
//   int i;
//   char cc;
//   for(i=0;i<n;i++)
//   {
//     if((n%31)) Rprintf("%s\t%s\n",x_lab,y_lab);
//     if((n%32)) scanf("continue? %c?\n",&cc);
//     Rprintf("%g\t %g\n", *(x + i),*(y + i));
//   }
// }

int cmpr_IND_VAL(const void *xx, const void *yy)
{
  IND_VAL *x, *y;
  x = (IND_VAL *)xx;
  y = (IND_VAL *)yy;
  return(1*(x->value > y->value) - 1*(x->value < y->value));
}
