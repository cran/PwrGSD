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
#include<R.h>
#include<Rmath.h>

#define normut 8.20953615160139
#define MAX(x,y) (x > y ? x : y)
#define wchidx(x,xgrid,n,l) l=0;while(*(xgrid+l)<=x && l<n) l++

void printmat(double *pA, int nr, int nc, char *name);
void printmati(int *pA, int nr, int nc, char *name);

void grpseqbnds(int *dofu, int *nbf, int *nbnd, int *nsf, double *rho, int *pnthslook,
                double *palphatot, double *palpha, double *psimin, int *dlact, 
		double *pInfTold, double *pInfTnew, double *pInfTold_ii, double *pInfTnew_ii, 
		double *px, double *py, double *ptmp, double *pintgrndx, double *pgqxw, 
		int *pngqnodes, double *mufu, double *pbold, double *pbnew, int *mybound, int *prev);

void StCu2Bnds(double *pmu,double *pfrac,double *palpha,int *psided,double *prho,int *pef,double *b);

void cmptpwrgsd(int *dofu,int *pnlook,double *pfrac_k,double *pfrac_kp1,double *mu_o,
                double *mu_n,double *Psiab_o,double *Psiab_n,double *Psiminfa_o,
                double *Psiminfa_n,double *Psibinf_o,double *Psibinf_n,double *gqxw,
                int *pngqnodes,double *bold,double *bnew,double *palpha0,double *palpha1);

// add header for alternate weighting function below 'flemhar':

// add header for alternate spending function below 'obrien' and 'pocock':

void driftfu(int *ints,double *accru,double *accrat,double *tlook,double *ppar,double *gqxw,
	     double *th0,double *h0,double *lrrf,double *thc0,double *hc0,double *thc1,
	     double *hc1,int *wttyp,double *mufu,int *puserVend,double *Vend);

void drift(int *ints,double *accru,double *accrat,double *tlook,double *ppar,double *gqxw,
	   double *th0,double *h0,double *th1,double *h1,double *thc0,double *hc0,
	   double *thc1,double *hc1,double *tlA0,double *lA0,double *tlB0,double *lB0,
	   double *tlA1,double *lA1,double *tlB1,double *lB1,double *thA0,double *hA0,
	   double *thB0,double *hB0,double *thA1,double *hA1,double *thB1,double *hB1,
	   int *wttyp,double *RR,int *pnjmp,double *InfFrac,double *InfFrac_ii,double *mu,
	   double *Var_uw,double *Var,double *Eta,int *puserVend,double *Vend);

// begin main
void AsyPwrGSD(int *ints,double *dbls,double *pttlook,double *palphatot,double *lrrf,
	       double *bHay,double *ppar,double *pgqxw,double *tcut0,double *h0,
	       double *tcut1,double *h1,double *tcutc0,double *hc0,double *tcutc1,
	       double *hc1,double *tcutd0A,double *hd0A,double *tcutd0B,double *hd0B,
	       double *tcutd1A,double *hd1A,double *tcutd1B,double *hd1B,double *tcutx0A,
	       double *hx0A,double *tcutx0B,double *hx0B,double *tcutx1A,double *hx1A,
	       double *tcutx1B,double *hx1B,int *wttyp,double *Vend,double *pinffrac,
	       double *pinffrac_ii,double *pbounds,double *mufu,double *mu,
	       double *palpha0vec,double *palpha1vec,double *RR,int *pnjmp,double *Var_uw,
	       double *Var, double *Eta, int *t_idx, double *betabdry, double *bstar)
{
  int *pn,*pntot,*pnthslook,*nbnd,*nsf,*gradual,*pnlook,*pnstat,*pngqnodes,*dofu,*dlact,*pncut0;
  int *pncut1,*pncutc0,*pncutc1,*sided,*pncutd0A,*pncutd0B,*pncutd1A,*pncutd1B,*pef,*pncutx0A;
  int *pncutx0B,*pncutx1A,*pncutx1B,*puserVend,*mybounds,*spend_info_k,*qis1orQ,*prev,*nbf;
  int nlook,ncut0,ncut1,ncut,ii,j,k,kacte,kactf,l,ngqnodes,nstat,istat,ef,ixxx,flag,idx;
  int ijmp,njmp,userhazfu,spend_info,krchd_flag,nbnd_e_sv,nbnd_f_sv;

  double *ptlook,*val,*pInfTold,*pInfTnew,*pInfTold_ii,*pInfTnew_ii,*psimin,*palpha,*pbold,*pbnew;
  double *px,*py,*ptmp,*pintgrndx,*mu_o,*mu_n,*Psiab_o,*Psiab_n,*Psiminfa_o,*Psiminfa_n,*Psibinf_o;
  double *Psibinf_n,*palpha10,*palpha11,*q,*tjmp,*rr,*CS_q_dEta,*rho,*rho_sc,*accru,*accrat,*mufuforSC;
  double atotsv_e,atotsv_f,fmin,sqrf,s,stlde,b,Eta_old,CS_,xntrial,xntrialhlf,V_,f1,dt;
  double dEta,fp,f_k,a_tmp,b_tmp,f_krchd,f_krchd_ii,f_k_ii,vend;

/* dbls <- c(rho.Efficacy,rho.Futility,accru,accrat,spend.info.p) */
/* void AsyPwrGSD(int *ints,double *rho,double *spend_info_p,double *accrual,*/

  pnlook       = ints;
  nlook        = *pnlook;
  pnstat       = ints +  1;
  pngqnodes    = ints +  2;
  pncut0       = ints +  3;
  pncut1       = ints +  4;
  pncutc0      = ints +  5;
  pncutc1      = ints +  6;
  pncutd0A     = ints +  7;
  pncutd0B     = ints +  8;
  pncutd1A     = ints +  9;
  pncutd1B     = ints + 10;
  pncutx0A     = ints + 11;
  pncutx0B     = ints + 12;
  pncutx1A     = ints + 13;
  pncutx1B     = ints + 14;
  gradual      = ints + 15;
  nbnd         = ints + 16;            /* 16           through 16+2*nlook-1 */
  nsf          = ints + 16+2*nlook;    /* 16+2*nlook   through 16+4*nlook-1 */
  dofu         = ints + 16+4*nlook;
  userhazfu    = *(ints+16+4*nlook+1);
  spend_info   = *(ints+16+4*nlook+2);
  puserVend    = ints + 16+4*nlook+3;
  mybounds     = ints + 16+4*nlook+4;  /* 16+4*nlook+4 through 16+4*nlook+5 */
  spend_info_k = ints + 16+4*nlook+6;
  qis1orQ      = ints + 16+4*nlook+7; 
  sided        = ints + 16+4*nlook+8;
  nbf          = ints + 16+4*nlook+9;

  accru  = dbls;
  accrat = dbls + 1;
  rho    = dbls + 2;          /* 2         through 2+2*nlook-1 */
  rho_sc = dbls + 2+2*nlook;  /* 2+2*nlook through 2+4*nlook-1 */

  /* for now */
  *(nbnd+1)   = *(nbnd + nlook);
  *(nsf+1)    = *(nsf + nlook);
  *(rho+1)    = *(rho + nlook);
  *(rho_sc+1) = *(rho_sc + nlook);
  /* until I implement changing boundary types */


  ngqnodes = *pngqnodes;
  nstat = *pnstat;

  xntrial = *accru * *accrat;
  xntrialhlf = pow(xntrial,0.5);

  pef        = (int *)Calloc(1,int);
  pnthslook  = (int *)Calloc(2,int);
  dlact      = (int *)Calloc(2,int);
  prev       = (int *)Calloc(1,    int);

  px         = (double *)Calloc(2*ngqnodes,double);
  py         = (double *)Calloc(2*ngqnodes,double);
  ptmp       = (double *)Calloc(2*ngqnodes,double);
  pintgrndx  = (double *)Calloc(2*ngqnodes,double);
  Psiminfa_o = (double *)Calloc(ngqnodes,double);
  Psiminfa_n = (double *)Calloc(ngqnodes,double);
  ptlook     = (double *)Calloc(1,double);
  pInfTold   = (double *)Calloc(2,double);
  pInfTnew   = (double *)Calloc(1,double);
  pInfTold_ii= (double *)Calloc(2,double);
  pInfTnew_ii= (double *)Calloc(1,double);
  palpha     = (double *)Calloc(2,double);
  palpha10   = (double *)Calloc(1,double);
  palpha11   = (double *)Calloc(1,double);
  pbold      = (double *)Calloc(2,double);
  pbnew      = (double *)Calloc(2,double);
  psimin     = (double *)Calloc(1,double);
  mu_o       = (double *)Calloc(1,double);
  mu_n       = (double *)Calloc(1,double);
  mufuforSC  = (double *)Calloc(2,double);

  *prev = 0;

  drift(ints,accru,accrat,pttlook,ppar,pgqxw,tcut0,h0,tcut1,h1,tcutc0,hc0,tcutc1,hc1,
	tcutd0A,hd0A,tcutd0B,hd0B,tcutd1A,hd1A,tcutd1B,hd1B,tcutx0A,hx0A,tcutx0B,hx0B,
	tcutx1A,hx1A,tcutx1B,hx1B,wttyp,RR,pnjmp,pinffrac,pinffrac_ii,mu,Var_uw,Var,Eta,
	puserVend,Vend);

  /* Error probability spending information on Variance Scale                                  */
  if(spend_info==0)
    for(k=0;k<nlook*nstat;k++) *(pinffrac_ii + k) = *(pinffrac + k);


  /* Error probability spending information on Events Scale                                    */
  /* IF(SPEND_INFO==1) do nothing; pinffrac_ii already equals prop'n events                    */


  /* Hybrid: linear switch from Variance Scale to Events Scale starting after analysis '*spend_info_k' */
  if(spend_info==2)
    for(j=0;j<nstat;j++){
      krchd_flag=0;
      f_krchd = 0.0;
      f_krchd_ii = 0.0;
      for(k=0;k<nlook;k++){
        f_k = *(pinffrac + nlook*j + k);
        f_k_ii = *(pinffrac_ii + nlook*j + k);
	if(k >= *spend_info_k){
	  if(krchd_flag==0){
            f_krchd = *(pinffrac + nlook*j + *spend_info_k);
            f_krchd_ii = *(pinffrac_ii + nlook*j + *spend_info_k);
	    krchd_flag=1;
	  }
          *(pinffrac_ii + nlook*j + k) = f_krchd +  (1.0 - f_krchd)/(1.0 - f_krchd_ii) * (f_k_ii - f_krchd_ii);
	}
	else  *(pinffrac_ii + nlook*j + k) = f_k;
      }
    }
  /* Error probability spending information on Calender Time Scale                              */
  if(spend_info==3)
    for(j=0;j<nstat;j++){
      for(k=0;k<nlook;k++)
        *(pinffrac_ii + nlook*j + k) = *(pttlook + k)/(*(pttlook + nlook - 1));    
    }

  for(l=0;l<nlook*nstat;l++) *(mufu + l) = 0.0;

  if(*dofu == 1){
    if(userhazfu==0)
      driftfu(ints,accru,accrat,pttlook,ppar,pgqxw,tcut0,h0,lrrf,tcutc0,hc0,tcutc1,hc1,
	      wttyp,mufu,puserVend,Vend);
     else
       for(l=0;l<nlook*nstat;l++) *(mufu+l) = *(mu+l);
  }

  njmp = *pnjmp;
  tjmp = (double *)Calloc(njmp,double);
   
  rr = (double *)Calloc(njmp,double);
  CS_q_dEta = (double *)Calloc(nstat*njmp,double);
  q = (double *)Calloc(nstat*njmp,double);

  if(*sided==-1 || *sided==-2)    
    for(l=0;l<nlook*nstat;l++) {
      *(mufu+l) = *(mufu+l) * (-1.0);
      *(mu+l) = *(mu+l) * (-1.0);
    }

  atotsv_e = *palphatot;
  atotsv_f = *(palphatot+1);

  nbnd_e_sv = *nbnd;
  nbnd_f_sv = *(nbnd + 1);

  *psimin = 1.0e-15;

  for(j=0;j<nstat;j++){
    k=0;
    kacte =0;
    kactf =0;
    *dlact=0;
    *(dlact+1)=0;

    *pInfTold = 0.0;
    *(pInfTold + 1)=0.0;

    *pInfTold_ii = 0.0;
    *(pInfTold_ii+1) = 0.0;

    if(*mybounds==0 || nbnd_e_sv==3){
      if(*nbnd==1 || nbnd_e_sv==3) *pbold = normut;
      if(*nbnd==2) *pbold = *bHay;
    }

    if(*(mybounds+1)==0 || nbnd_f_sv==3){
      if(*(nbnd+1)==1 || nbnd_f_sv==3) *(pbold+1) = -normut;
      if(*(nbnd+1)==2) *(pbold+1) = *bHay;
    }

    *palphatot = atotsv_e;
    *(palphatot + 1) = atotsv_f;
    *palpha = 0.0;
    *(palpha+1) = 0.0;
    for(l=0;l<2*ngqnodes;l++) {
      *(px+l) = 0.0;
      *(py+l) = 0.0;
      *(ptmp+l) = 0.0;
      *(pintgrndx+l) = 0.0;
    }
    while(k<nlook){
      *pInfTnew = *(pinffrac + nlook*j + k);
      *pInfTnew_ii = *(pinffrac_ii + nlook*j + k);
      /* use the Stochastic Curtailment procedure to construct the efficacy boundary */
      if(nbnd_e_sv==3){
        *pef = 0;
        *mufuforSC = *(mufu + nlook*j + k);
        *(mufuforSC + 1) = *(mufu + nlook*j + nlook - 1);
        StCu2Bnds(mufuforSC,pInfTnew,palphatot,sided,rho_sc,pef,pbounds+nlook*j+k);
        *mybounds = 1;
        *nbnd = 1;
      }
      /* use the Stochastic Curtailment procedure to constuct the futility boundary */
      if(nbnd_f_sv==3){
        *pef = 1;
        *mufuforSC = *(mufu + nlook*j + k);
        *(mufuforSC + 1) = *(mufu + nlook*j + nlook - 1);
        StCu2Bnds(mufuforSC,pInfTnew,palphatot,sided,rho_sc+1,pef,pbounds+nstat*nlook+nlook*j+k);
        /*-----------------------------------------------------------------------------------------------------------------------|
        | Rprintf("j:%d, k:%d, mu_k:%g, mu_end:%g, fnew:%g, alpha:%g, sided:%d, rho:%g, pef:%d, b:%g\n",                         |
        |          j,k,*mufuforSC,*(mufuforSC+1),*pInfTnew,* palphatot,*sided,*(rho_sc+1),*pef,*(pbounds+nstat*nlook+nlook*j+k));|
        |-----------------------------------------------------------------------------------------------------------------------*/
        *(mybounds + 1) = 1;
        *(nbnd + 1) = 1;
      }
      *pbnew = normut;
      *(pbnew+1) = -normut;
      if(*mybounds==1) *pbnew = *(pbounds + nlook*j + k);
      if(*(mybounds+1)==1) *(pbnew+1) =  *(pbounds + nstat*nlook + nlook*j + k);
      
      *ptlook = *(pttlook+k);
      *pnthslook = kacte+1;
      *(pnthslook+1) = kactf+1;

      flag = 0;

      grpseqbnds(dofu,nbf,nbnd,nsf,rho, pnthslook,palphatot,palpha,psimin,
                 dlact,pInfTold,pInfTnew,pInfTold_ii,pInfTnew_ii,px,py,ptmp,
		 pintgrndx,pgqxw,pngqnodes,mufu + nlook*j + k,pbold,pbnew,mybounds,prev);

      flag=1;
      if(*dlact==1){
        if (*nbnd == 1){
          b_tmp = (*mybounds ? *(pbounds + nlook*j + k): *pbnew);
	  *(pbounds + nlook*j + k) = b_tmp;
          *pbold = b_tmp;
	}
        if (*nbnd == 2) {
          if (1.0-*pInfTnew_ii >= 1e-6) {
            *(pbounds + nlook*j + k) = *bHay;
            *pbold = *bHay;
            *palphatot = *palphatot - *palpha;
          }
          else *(pbounds + nlook*j + k) = *pbnew;
        }
        *(palpha0vec + nlook*j + k) = *palpha;
        *pInfTold = *pInfTnew;
        *pInfTold_ii = *pInfTnew_ii;
        kacte += 1;
      }
      else{
	b_tmp = (*mybounds ? *(pbounds + nlook*j + k) : normut);
        *(pbounds + nlook*j + k) = b_tmp;
        *pbold = b_tmp;
        *(palpha0vec + nlook*j + k) = *psimin;
/*
        *pInfTold = *pInfTnew;
*/
        *pInfTold_ii = *pInfTnew_ii;

      }
      if(*dofu==1 && *(dlact+1)==1){
        if(*(nbnd+1)==1) {
	  b_tmp = (*(mybounds+1) && 1.0 - *pInfTnew_ii >= 1e-6 ? *(pbounds + nstat*nlook + nlook*j + k) : *(pbnew+1));
          *(pbounds + nstat*nlook + nlook*j + k) = b_tmp;
          *(pbold+1) = b_tmp;
        }
        if(*(nbnd+1)==2) {
          if(1.0-*pInfTnew_ii >= 1e-6){
            *(pbounds + nstat*nlook + nlook*j + k) = *bHay;
            *(pbold+1) = *bHay;
            *(palphatot+1) = *(palphatot+1) - *(palpha+1);
          }
          else *(pbounds + nstat*nlook + nlook*j + k) = *(pbnew+1);
        }
        *(palpha0vec + nstat*nlook + nlook*j + k) = *(palpha+1);
        *(pInfTold + 1) = *pInfTnew;
	*(pInfTold_ii + 1) = *pInfTnew_ii;
        kactf += 1;
      }
      if(*dofu==1 && *(dlact+1)==0){
	b_tmp = (*(mybounds+1) ? *(pbounds + nstat*nlook + nlook*j + k) : -normut);
        *(pbounds + nstat*nlook + nlook*j + k) = b_tmp;
        *(pbold+1) = b_tmp;
        *(palpha0vec + nstat*nlook + nlook*j + k) = *psimin;
/*
        *(pInfTold + 1) = *pInfTnew;
*/
	*(pInfTold_ii + 1) = *pInfTnew_ii;	

      }
      *(pinffrac + nlook*j + k) = *pInfTnew;
      *(pinffrac_ii + nlook*j + k) = *pInfTnew_ii;
      k++;
    }//end while loop
  } //different stats loop end

  for(j=0;j<nstat;j++){
    for(l=0;l<ngqnodes;l++) {
      *(px+l) = 0.0;
      *(py+l) = 0.0;
      *(ptmp+l) = 0.0;
      *(pintgrndx+l) = 0.0;
      *(Psiminfa_o+l) = 0.0;
      *(Psiminfa_n+l) = 0.0;
    }
    *palpha10 = 0.0;
    *palpha11 = 0.0;
    Psiab_o = px;
    Psiab_n = py;
    Psibinf_o = ptmp;
    Psibinf_n = pintgrndx;

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
      cmptpwrgsd(dofu,pnthslook,pInfTold,pInfTnew,mu_o,mu_n,Psiab_o,Psiab_n,
		 Psiminfa_o,Psiminfa_n,Psibinf_o,Psibinf_n,pgqxw,
		 pngqnodes,pbold,pbnew,palpha10,palpha11);

      *(palpha1vec + nlook*nstat + nlook*j + k) = *palpha10;
      *(palpha1vec + nlook*j + k) = *palpha11;

      *pInfTold=*pInfTnew;
      *mu_o = *mu_n;
      *pbold = *pbnew;
      *(pbold + 1) = *(pbnew + 1);
      k++;
    }//end while loop
  } //different stats loop end

  for(istat=0;istat<nstat;istat++)
  {
    s = 0.0;
    stlde = 0.0;
    Eta_old = 0.0;
    for(ijmp=0;ijmp<njmp;ijmp++) {
      *(tjmp + ijmp) = *(RR + ijmp);
      dEta = *(Eta + njmp*istat + ijmp) - Eta_old;
      Eta_old = *(Eta + njmp*istat + ijmp);
      *(rr + ijmp) = *(RR + 2*njmp + ijmp);
      stlde += log(*(rr + ijmp)) * dEta;
      s += log(*(RR + 4*njmp + ijmp)) * dEta;
    }
    *(bstar + istat)= stlde/Eta_old;
    *(bstar + nstat + istat) = s/Eta_old;
    for(ijmp=0;ijmp<njmp;ijmp++) 
      *(q + njmp*istat + ijmp) = log(*(rr + ijmp))/(*(bstar+istat));
  }

  for(k=0;k<nlook;k++){
    wchidx(*(pttlook+k),tjmp,njmp,l);
    *(t_idx + k) = l-1;
  }
  for(istat=0;istat<nstat;istat++){  
    Eta_old = 0.0;
    CS_ = 0.0;
    for(ijmp=0;ijmp<njmp;ijmp++){
      CS_ += *(q + njmp*istat + ijmp) * (*(Eta + njmp*istat + ijmp) - Eta_old);
      *(CS_q_dEta + njmp*istat + ijmp) = CS_; 
      Eta_old = *(Eta + njmp*istat + ijmp);
    }
  }
  
  for(ef=0;ef<2;ef++)
    for(istat=0;istat<nstat;istat++)
    {
      vend = *(Var + njmp*istat + *(t_idx + nlook - 1));
      for(k=0;k<nlook;k++)
      {
        b = *(pbounds + nstat*nlook*ef + nlook*istat + k) * (1.0 - 2 * (*sided==-1||*sided==-2));
        *(betabdry + nstat*nlook*ef + nlook*istat + k) = 
          pow(vend,0.5) * b/(xntrialhlf* *(CS_q_dEta + njmp*istat + *(t_idx + k)));
      }
    }

  *nbnd = nbnd_e_sv;
  *(nbnd + 1) = nbnd_f_sv;

  Free(px);
  Free(py);
  Free(ptmp);
  Free(pintgrndx);
  Free(Psiminfa_o);
  Free(Psiminfa_n);
  Free(ptlook);
  Free(pInfTold);
  Free(pInfTnew);
  Free(pInfTold_ii);
  Free(pInfTnew_ii);
  Free(palpha);
  Free(palpha10);
  Free(palpha11);
  Free(pbold);
  Free(pbnew);
  Free(pnthslook);
  Free(pef);
  Free(prev);
  Free(psimin);
  Free(dlact);
  Free(mu_o);
  Free(mu_n);
  Free(tjmp);
  Free(rr);
  Free(CS_q_dEta);
  Free(q);
  Free(mufuforSC);
}

void printmat(double *pA,int nr,int nc,char *name)
{
  int j,k;

  Rprintf("%s = \n",name);
  for (j=0;j<nr;j++) {
    for (k=0;k<nc;k++) Rprintf("%g ",*(pA + nr*k + j));
    Rprintf("\n");
  }
}

void printmati(int *pA,int nr,int nc,char *name)
{
  int j,k;

  Rprintf("%s = \n",name);
  for (j=0;j<nr;j++) {
    for (k=0;k<nc;k++) Rprintf("%d ",*(pA + nr*k + j));
    Rprintf("\n");
  }
}
