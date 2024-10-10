#include<R.h>
#include<Rmath.h>
#include<R_ext/Utils.h>
#include "PwrGSD_mem.h"

void cpy(double *x, double *y, int n);
void zero(double *x, int n);

void gsd_dens(double *frac, double *xeff, double *gqxw, int *pngq, int *pnlooks, double *xx, double *dF, double *x1c, double *dF1c)
{
  int ngq, nlooks, i, k, j;
  double x, Phi_x, u, du, fkm1, fk, sqrDfk, x_r_k, dx_r_k, x_c_km1_, dx_c_km1, x_c_k_, dx_c_k, dx_a_k;
  double *x_c_km1, *dF_c_km1, *x_c_k, *dF_c_k, *x_c_km1_sv, *dF_c_km1_sv;

  ngq = *pngq;
  nlooks = *pnlooks;

  x_c_km1     = Calloc(ngq, double);
  x_c_km1_sv  = Calloc(ngq, double);
  dF_c_km1    = Calloc(ngq, double);
  dF_c_km1_sv = Calloc(ngq, double);
  x_c_k       = Calloc(ngq, double);
  dF_c_k      = Calloc(ngq, double);

  x = *xeff;
  Phi_x = pnorm5(x, 0.0, 1.0, 1, 0);
  fkm1=0.0;
  fk = *frac;
  sqrDfk = pow(fk-fkm1, 0.5);

  dx_r_k = dx_a_k = 1.0;

  for(i=0;i<ngq;i++)
  {
    u = (1.0+ (*(gqxw+i)))/2.0;
    du = (*(gqxw+ngq+i))/2.0;

    *(xx + i) = x_r_k = qnorm(Phi_x*(1-u) + u, 0.0, 1.0, 1, 0);
    
    dx_r_k = (1.0 - Phi_x)*du/dnorm4(x_r_k, 0.0, 1.0, 0.0);

    *(x_c_km1 + i) = x_c_km1_ = qnorm(Phi_x*u, 0.0, 1.0, 1, 0);
    dx_c_km1 = Phi_x*du/dnorm4(x_c_km1_, 0.0, 1.0, 0.0);

    *(dF + i) = dnorm(x_r_k/sqrDfk, 0.0, 1.0, 0)/sqrDfk*dx_r_k;
    *(dF_c_km1 + i) = dnorm4(x_c_km1_/sqrDfk, 0.0, 1.0, 0)/sqrDfk*dx_c_km1;
  }
  cpy(x_c_km1, x1c, ngq);
  cpy(dF_c_km1, dF1c, ngq);

  for(k=1;k<nlooks;k++)
  {
    x = *(xeff + k);
    Phi_x = pnorm5(x, 0.0, 1.0, 1, 0);
    fkm1 = fk;
    fk = *(frac + k);
    sqrDfk = pow(fk - fkm1, 0.5);
    zero(dF_c_k, ngq);
    for(i=0;i<ngq;i++)
    {
      u = (1.0+ (*(gqxw+i)))/2.0;
      du = (*(gqxw+ngq+i))/2.0;

      *(xx + k*ngq + i) = x_r_k = qnorm(Phi_x*(1-u) + u, 0.0, 1.0, 1, 0);
      dx_r_k = (1.0 - Phi_x)*du/dnorm(x_r_k, 0.0, 1.0, 0.0);
      *(x_c_k + i) = x_c_k_ = qnorm(Phi_x*u, 0.0, 1.0, 1, 0);
      dx_c_k = Phi_x*du/dnorm(x_c_k_, 0.0, 1.0, 0.0);
      for(j=0;j<ngq;j++)
      {
        *(dF + k*ngq + i) = *(dF + k*ngq + i) + dnorm4((x_r_k-(*(x_c_km1+j)))/sqrDfk,0.0,1.0,0)/sqrDfk*(*(dF_c_km1+j))*dx_r_k;
        *(dF_c_k + i) = *(dF_c_k + i) + dnorm4((x_c_k_-(*(x_c_km1+j)))/sqrDfk,0.0,1.0,0)/sqrDfk*(*(dF_c_km1+j))*dx_c_k;
      }
    }
    cpy(x_c_km1, x_c_km1_sv, ngq);
    cpy(dF_c_km1, dF_c_km1_sv, ngq);
    cpy(x_c_k, x_c_km1, ngq);
    cpy(dF_c_k, dF_c_km1, ngq);
  }
  
  Free(x_c_km1);
  Free(x_c_km1_sv);
  Free(dF_c_km1);
  Free(dF_c_km1_sv);
  Free(x_c_k);
  Free(dF_c_k);
}

void cpy(double *x, double *y, int n)
{
  int i;
  for(i=0;i<n;i++) *(y + i) = *(x + i);
}

void zero(double *x, int n)
{
  int i;
  for(i=0;i<n;i++) *(x + i) = 0.0;
}
