/*  Copyright (C) 2004	    Grant Izmirlian
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  Synopsis:
 */
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<R.h>
#include<Rmath.h>
#define normut 8.20953615160139

void cmptdnsty(int *dofu, int *pnlook, double *pfrac_k, double *pfrac_kp1,
               double *mu_o, double *mu_n, double *Psiab_o, double *Psiab_n,
               double *Psiminfa_o, double *Psiminfa_n, double *Psibinf_o,
               double *Psibinf_n, double *gqx, double *gqw, int *pngqnodes,
               double *bold, double *bnew, double *store, double *xgrid, 
	       double *dxgrid, double *dens);

// add header for alternate weighting function below 'flemhar':

// add header for alternate spending function below 'obrien' and 'pocock':

// begin main
void Dnsties(int *dofu,int *pnlook,int *pnstat,int *sided,
		 int *pngqnodes,double *pgqx,double *pgqw,double *pinffrac,
		 double *pbounds,double *mu,double *xgrids,double *dxgrids,
		 double *denss)
{
  int *pnthslook;
  int nlook, j, k, l, i, ngqnodes, ngq2, nstat;

  double *pInfTold,*pInfTnew,*pbold,*pbnew,*mu_o,*mu_n,*Psiab_o,*Psiab_n,*Psiminfa_o;
  double *Psiminfa_n,*Psibinf_o,*Psibinf_n,*store,*xgrid,*dxgrid,*dens;

  ngqnodes = *pngqnodes;
  ngq2 = 2*ngqnodes;
  nstat = *pnstat;
  nlook = *pnlook;

  Psiab_o = (double *)Calloc(ngqnodes, double);
  Psiab_n = (double *)Calloc(ngqnodes, double);
  Psibinf_o = (double *)Calloc(ngqnodes, double);
  Psibinf_n = (double *)Calloc(ngqnodes, double);
  Psiminfa_o = (double *)Calloc(ngqnodes, double);
  Psiminfa_n = (double *)Calloc(ngqnodes, double);
  pInfTold = (double *)Calloc(1, double);
  pInfTnew = (double *)Calloc(1, double);
  pbold = (double *)Calloc(2, double);
  pbnew = (double *)Calloc(2, double);
  pnthslook = (int *)Calloc(2, int);
  mu_o = (double *)Calloc(1, double);
  mu_n = (double *)Calloc(1, double);
  store = (double *)Calloc(1, double);
  xgrid = (double *)Calloc(ngq2, double);
  dxgrid = (double *)Calloc(ngq2, double);
  dens = (double *)Calloc(ngq2, double);

  for(j=0;j<nstat;j++){
    for(l=0;l<ngqnodes;l++) {
      *(Psiab_o+l) = 0.0;
      *(Psiab_n+l) = 0.0;
      *(Psibinf_o+l) = 0.0;
      *(Psibinf_n+l) = 0.0;
      *(Psiminfa_o+l) = 0.0;
      *(Psiminfa_n+l) = 0.0;
    }

    k=0;
    *pInfTold = *(pinffrac + nlook*j);
    *mu_o = *(mu + nlook*j);
    *pbold = *(pbounds + nlook*j);
    *(pbold + 1) = *(pbounds + nlook*nstat + nlook*j);
    while(k<nlook){
      if(k<nlook-1) {
	*pInfTnew = *(pinffrac + nlook*j + k + 1);
	*mu_n = *(mu + nlook*j + k + 1);
	*pbnew = *(pbounds + nlook*j + k + 1);
	*(pbnew + 1) = *(pbounds + nlook*nstat + nlook*j + k + 1);
      }	
      else{
	*pInfTnew = 2.0;
	*mu_n = 0.0;
	*pbnew = normut;
	*(pbnew + 1) = -normut;
      }
      *pnthslook = k+1;
      cmptdnsty(dofu,pnthslook,pInfTold,pInfTnew,mu_o,mu_n,Psiab_o,Psiab_n, 
		Psiminfa_o,Psiminfa_n,Psibinf_o,Psibinf_n,pgqx,pgqw,
		pngqnodes,pbold,pbnew,store,xgrid,dxgrid,dens);

      for(i=0;i<ngq2;i++){
        *(xgrids  + ngq2*nlook*j + ngq2*k + i) = *(xgrid   + i);
	*(dxgrids + ngq2*nlook*j + ngq2*k + i) = *(dxgrid  + i);
	*(denss   + ngq2*nlook*j + ngq2*k + i) = *(dens    + i);
      }

      *pInfTold=*pInfTnew;
      *mu_o = *mu_n;
      *pbold = *pbnew;
      *(pbold + 1) = *(pbnew + 1);
      k++;
    }//end while loop
  } //different stats loop end

  Free(Psiminfa_o);
  Free(Psiminfa_n);
  Free(pInfTold);
  Free(pInfTnew);
  Free(pbold);
  Free(pbnew);
  Free(pnthslook);
  Free(mu_o);
  Free(mu_n);
  Free(store);
  Free(xgrid);
  Free(dxgrid);
  Free(dens);
}
