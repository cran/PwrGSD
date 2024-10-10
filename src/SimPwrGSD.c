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
#define EPS 1e-10
#include<R.h>
#include <R_ext/Utils.h>
#include "PwrGSD_mem.h"

/* MACROS */
#define COMPH(xh,h,H,n,l) *H=0.0;for(l=1;l<n;l++) *(H+l)=*(H+l-1)+ *(h+l-1) * (*(xh+l) - *(xh+l-1))
#define hatX(x,xh,h,h_,n,l) l=0;while((*(xh+l)<=x) && (l<n)) l++; h_ = *(h+l-1)
#define hHatX(x,xh,h,H,h_,H_,n,l) l=0;while((*(xh+l))<=x && (l<n)) l++; h_ = *(h+l-1); H_ = *(H+l-1) + h_ * (x-*(xh+l-1))
#define HIatW(W,xh,h,H,HI_,n,l) l=0;while((*(H+l)<=W) && (l<n)) l++; HI_ = *(xh+l-1) + (W - *(H+l-1))/(*(h+l-1))
#define MIN(x,y) (x<=y ? x : y)
#define MAX(x,y) (x>=y ? x : y)

#define normut 8.20953615160139

#include "PwrGSD.h"

typedef struct{
  int index;
  double time;
  int event;
  int arm;
} itea;

WtFun flemhar, sflemhar, ramp, *wtfun;

CmprFun compitea, CmprDbl, *cmprdbl, *f;

void    randfromh(int *pn, double *tcut, double *h, int *pncut, 
                  double *t);
void    randhcdtl(int *pn, double *tcut, double *h, int *pncut, double *tend,
	          double *tcutdA, double *hdA, int *pncutdA, double *tdA,
	          double *tcutdB, double *hdB, int *pncutdB, double *tdB,
                  double *tcutxA, double *hxA, int *pncutxA, double *tcutxB, 
		  double *hxB, int *pncutxB, int *gradual, int *code, double *t);
void       handle(int *pn, double *tlook,  double *u, double *t0, double *t1,
                  double *tc0, double *tc1, itea *YY, int *pntot, int *pntimes);

void      driftfu(int *ints,double *accru,double *accrat,double *tlook,double *ppar,double *gqxw,
                  double *th0,double *h0,double *lrrf,double *thc0,double *hc0,double *thc1,
                  double *hc1,int *wttyp,double *mufu,int *puserVend,double *Vend);

void    cpblocked(itea *Yord, int *pn, double *time, int *nrisk, int *nevent, int *pntimes, int *pnevtypes, int *pnblocks);

void      commonx(double *x1,double *h1,int *pn1,double *x2,double *h2,int *pn2,
	          double *x,double *hh1,double *hh2,int *pn);

void       wlrstat(double *time, int *nrisk, int *nevent, double *wt, int *pntimes, double *UQ,
                   double *varQ, double *m1, double *UQt, double *varQt, double *var1t);
void       ISDstat(double *time, int *nrisk, int *nevent, int *pntimes, double *wt, double *stat,
                   double *var);
/*------------------------------------------------------------------------------------------------------------
   If you wish to add another choice of weight function, it must be coded in C, and pasted at the bottom of this 
   file, or placed into the src subdirectory of this package. The calling sequence of your code must follow the 
   prototype convention according to the type definition 'wlrstat', shown above. Then and add the name of the 
   function to the following list of names declared as type 'WtFun'. 
-------------------------------------------------------------------------------------------------------------*/

void   grpseqbnds(int *dofu, int *nbf, int *nbnd, int *nsf, double *rho, int *pnthslook, 
                  double *palphatot, double *palpha, double *psimin, int *dlact, 
		  double *pInfTold, double *pInfTnew, double *pInfTold_ii, double *pInfTnew_ii, 
		  double *px, double *py, double *ptmp, double *pintgrndx, double *pgqxw, 
		  int *pngqnodes, double *mu, double *pbold, double *pbnew, int *mybounds);

void StCu2Bnds(double *pmu, double *pfrac, double *pzcrit, double *prho,int *pef, double *b);
void project_end(double *u, double *t0, double *t1, double *tc0, double *tc1, int *pn, double *t_proj, 
                 double *v_Tend_proj, double *m_Tend_proj);

void     printmat(double *pA, int nr, int nc, char *name);
void     printmati(double *pA, int nr, int nc, char *name);
void     printbtre(double *time, int *nrisk, int *nevent, int nt, int fromi, int toi, int nb);
int R_isnancpp(double x);

/* BEGIN MAIN */
void    SimPwrGSD(int *ints,double *dbls, double *pttlook,double *palphatot,double *lrrf,
                  double *bHay, double *bendSC, double *ppar,double *pgqxw,double *tcut0,double *h0,
		  double *tcut1, double *h1,double *tcutc0,double *hc0,double *tcutc1,double *hc1,
		  double *tcutd0A, double *hd0A,double *tcutd0B,double *hd0B,double *tcutd1A,
		  double *hd1A, double *tcutd1B,double *hd1B,double *tcutx0A,double *hx0A,
		  double *tcutx0B,double *hx0B,double *tcutx1A,double *hx1A,double *tcutx1B,
		  double *hx1B,double *t0,double *t1,double *tc0,double *tc1,double *td0A,
		  double *td0B,double *td1A, double *td1B,int *code,double *u,double *TOS,
		  int *Event,int *Arm,double *time, int *nrisk,int *nevent,int *pntimes,
		  double *avginffrac,double *avginffrac_ii, double *avgbounds,double *mufu,
		  double *palphavec,int *pRejAcc,int *kstop, double *duration,double *pStatTend,
		  double *pVarTend,double *pmTend1, double *pstatK,double *pvarK,double *v_Tend_proj,
		  double *m_Tend_proj)
{
  int *pntimesk,*pn,*pntot,*pnthslook,*nbnd,*nsf,*pnsim,*pnlook,*pnstat,*pngqnodes,*pncut0;
  int *pncut1,*pncutc0,*pncutc1,*pncutd0A,*pncutd0B,*pncutd1A,*pncutd1B,*pncutx0A,*pncutx0B;
  int *pncutx1A,*pncutx1B,*psided,*gradual,*dofu,*dlact,*pnblocks,*nrisk_,*nevent_,*puserVend;
  int *spend_info_k,*mybounds,*pef,*nbf,*stattype,*wttyp,*do_proj;
  int ngq2,nstlk,nstlk2,ntrial,n,nlook,sided_,ii,i,j,k,kacte,kactf,l,ntimes;
  int ntimesk,ngqnodes,nsim,nstat,RejNull,AccNull,totev,totev_k,idx,nppar=0,csumnppar;
  int userhazfu,spend_info,krchd_flag,evnts_krchd,nbnd_e_sv,nbnd_f_sv,pnevty,isbad=0,do_proj_;

  itea *Yord;

  double *pstatk,*pvark,*pinffrac,*pinffrac_ii,*pbounds,*ptlook,*wt,*pInfTold,*pInfTnew,*pInfTold_ii;
  double *pInfTnew_ii,*par,*palpha,*pbold,*pbnew,*px,*py,*ptmp,*pintgrndx,*statk,*vark,*m1k,*stat,*var,*m1;
  double *etaold,*etanew,*psimin,*rho,*accru,*accrat,*time_,*Vend,*rho_sc,*mufuforSC; 
  double *UQt,*varQt,*var1t,*t_proj;
  double atotsv_e,atotsv_f,wlrsgn,wlrZ,f_k,f_k_ii,tlook_=0.0,b_tmp,var_krchd,f_krchd;
  double f_krchd_ii, t_end;

  pnlook       = ints;
  nlook        = *pnlook;
  pnstat       = ints +  1;
  nstat        = *pnstat;
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

  nbnd         = ints + 16;          /* 16           through 16+2*nlook-1 */
  nsf          = ints + 16+2*nlook;  /* 16+2*nlook   through 16+4*nlook-1 */
  dofu         = ints + 16+4*nlook;
  userhazfu    = *(ints+16+4*nlook+1);
  spend_info   = *(ints+16+4*nlook+2);
  pnsim        = ints + 16+4*nlook+3;
  mybounds     = ints + 16+4*nlook+4; /* 16+4*nlook+4 through 16+4*nlook+5  */
  spend_info_k = ints + 16+4*nlook+6;
  //  qis1orQ      = ints + 16+4*nlook+7;
  psided       = ints + 16+4*nlook+8; 
  nbf          = ints + 16+4*nlook+9;
  stattype     = ints + 16+4*nlook+10;
  wttyp        = ints + 16+4*nlook+10+nstat;
  do_proj      = ints + 16+4*nlook+10+2*nstat;

  do_proj_ = *do_proj;
/* 
    ints <- c(nlook,nstat,NGaussQ,ncut0,ncut1,ncutc0,ncutc1,ncutd0A,ncutd0B,ncutd1A,
              ncutd1B,ncutx0A,ncutx0B,ncutx1A,ncutx1B,gradual,nbnd.e,nbnd.f,nsf.e,nsf.f,
              dofu,use.rhaz.fu,spend.info,Nsim,is.myE,is.myF,spend.info.k,qProp.one.or.Q,sided)

    dbls <- c(accru,accrat, rho.Efficacy,rho.Futility,rho.Eff.SC, rho.Fut.SC)
*/

  accru   = dbls;
  accrat  = dbls+1;
  rho     = dbls+2;          /* 2         through 2+2*nlook-1 */
  rho_sc  = dbls+2+2*nlook;  /* 2+2*nlook through 2+4*nlook-1 */
  t_end   = *(pttlook + nlook - 1);

  ngqnodes = *pngqnodes;
  ngq2 = 2*ngqnodes;
  //  pgqx = pgqxw;
  //  pgqw = pgqxw + ngqnodes;
  nsim = *pnsim;

  nstlk = nstat*nlook;
  nstlk2 = 2*nstat*nlook;
  ntrial = (int)((*accru)*(*accrat));
  ntrial = ntrial + (ntrial%2);
  n = ntrial/2;

  Yord        = Calloc(ntrial,  itea);
  pinffrac    = Calloc(nstlk, double);
  pinffrac_ii = Calloc(nstlk, double);
  pbounds     = Calloc(nstlk2,double);
  pstatk      = Calloc(nstat, double);
  pvark       = Calloc(nstat, double);
  par         = Calloc(3,     double);
  wt          = Calloc(ntrial, double);
  px          = Calloc(ngq2,  double);
  py          = Calloc(ngq2,  double);
  ptmp        = Calloc(ngq2,  double);
  pintgrndx   = Calloc(ngq2,  double);
  statk       = Calloc(1,     double);
  vark        = Calloc(1,     double);
  m1k         = Calloc(1,     double);
  stat        = Calloc(1,     double);
  var         = Calloc(1,     double);
  m1          = Calloc(1,     double);
  etaold      = Calloc(1,     double);
  etanew      = Calloc(1,     double);
  ptlook      = Calloc(1,     double);
  pInfTold    = Calloc(2,     double);
  pInfTnew    = Calloc(1,     double);
  pInfTold_ii = Calloc(2,     double);
  pInfTnew_ii = Calloc(1,     double);
  palpha      = Calloc(2,     double);
  pbold       = Calloc(2,     double);
  pbnew       = Calloc(2,     double);
  psimin      = Calloc(1,     double);
  Vend        = Calloc(nstat, double);
  mufuforSC   = Calloc(2,     double);
  t_proj      = Calloc(3,     double);

  pntimesk    = Calloc(1,     int);
  pnblocks    = Calloc(1,     int);
  pnthslook   = Calloc(2,     int);
  pn          = Calloc(1,     int);
  pntot       = Calloc(1,     int);
  dlact       = Calloc(2,     int);
  puserVend   = Calloc(1,     int);
  pef         = Calloc(1,     int);

  *pn = n;
  *pnblocks = 2;

  *psimin = 6.416208e-17;

  sided_ = *psided;

  for(l=0;l<2*ngqnodes;l++) {
    *(px+l) = 0.0;
    *(py+l) = 0.0;
    *(ptmp+l) = 0.0;
    *(pintgrndx+l) = 0.0;
  }
  *stat = 0.0;
  *var = 0.0;
  *m1 = 0.0;
  *statk = 0.0;
  *vark = 0.0;
  *m1k = 0.0;
  *palpha = 0.0;
  *(palpha+1) = 0.0;

  for(l=0;l<nstat;l++) *(Vend +l) = 1.0;
  *puserVend = 0;

  /* Calculate drift under design alternative based upon 'lrrf' which equals log(RR.Futility) */
  /* using the analytic variance  -- assuming  q(t) = 1 for all t                             */

  if(*dofu == 1)
    if(userhazfu==0)
      driftfu(ints,accru,accrat,pttlook,ppar,pgqxw,tcut0,h0,lrrf,tcutc0,hc0,tcutc1,hc1,
              wttyp,mufu,puserVend,Vend);

  for(j=0;j<nstat;j++) 
    if(sided_ < 0)
      for(l=0;l<nlook;l++) *(mufu + nlook*j + l) = *(mufu + nlook*j + l) * (-1.0);

  /* BEGIN simulation loop */
  GetRNGstate();
  ii=0;
  while(ii<nsim){
    /* times of death due to lung cancer in arm 0  */
    /* generated conditional upon td0Ai and td0Bi. */
    randhcdtl(pn, tcut0, h0, pncut0, pttlook+nlook-1, tcutd0A, hd0A, pncutd0A, 
	      td0A,tcutd0B, hd0B, pncutd0B, td0B, tcutx0A, hx0A, pncutx0A, 
	      tcutx0B, hx0B, pncutx0B, gradual, code, t0);

    /* times of death due to lung cancer in arm 1  */
    /* generated conditional upon td1Ai and td1Bi. */
    randhcdtl(pn, tcut1, h1, pncut1, pttlook+nlook-1, tcutd1A, hd1A, pncutd1A, 
	      td1A, tcutd1B, hd1B, pncutd1B, td1B, tcutx1A, hx1A, pncutx1A, 
              tcutx1B, hx1B, pncutx1B, gradual, code+n, t1);

    for(l=0;l<n;l++) *(code+n+l) = -1 * *(code+n+l);

    /* censoring times due to death of other cause, arms 0 and 1.    */
    /* censoring of this type as well as censoring due to alive      */
    /* at end of trial are handled outside in the function 'handle'. */
    randfromh(pn, tcutc0, hc0, pncutc0, tc0);
    randfromh(pn, tcutc1, hc1, pncutc1, tc1);

    for(l=0;l<ntrial;l++) *(u+l) = *accru * unif_rand();

    atotsv_e = *palphatot;
    atotsv_f = *(palphatot+1);

    nbnd_e_sv = *nbnd;
    nbnd_f_sv = *(nbnd + 1);

    csumnppar = 0;

    for(j=0;j<nstat;j++){
      /* ----------------------------------------------------------
         Add new weight functions to this hard-coded menu.
         note: integer vector 'wttyp' of length 'nstat' contains 
               weight function choices, currently 0, 1, or 2 
      -----------------------------------------------------------*/
      if(*(wttyp+j) == 0) {
	  wtfun = &flemhar;
	  nppar = 2;
      }
      if(*(wttyp+j) == 1) {
          wtfun = &sflemhar;
	  nppar = 3;
      }
      if(*(wttyp+j) == 2) {
          wtfun = &ramp;
	  nppar = 1;
      }

      for(k=0;k<nppar;k++)
        *(par +k) = *(ppar + csumnppar + k);

      *ptlook = *(pttlook+nlook-1);
      handle(pn, ptlook, u, t0, t1, tc0, tc1, Yord, pntot, pntimes);

      ntimes = *pntimes;

      time_   = Calloc(ntimes, double);
      nrisk_  = Calloc(2*ntimes, int);
      nevent_ = Calloc(2*ntimes, int);
      UQt     = Calloc(ntimes, double);
      varQt   = Calloc(ntimes, double);
      var1t   = Calloc(ntimes, double);

      pnevty = 1;
      cpblocked(Yord, pntot, time_, nrisk_, nevent_, pntimes, &pnevty, pnblocks);
 
      totev = 0;
      for(i=0;i<2*ntimes;i++) totev += *(nevent_ + i);

      (*wtfun)(time_, nrisk_, nevent_, pntimes, par, wt);
      if(*(stattype + j)==0) wlrstat(time_,nrisk_,nevent_,wt,pntimes,stat,var,m1,UQt,varQt,var1t);
      if(*(stattype + j)==1) ISDstat(time_, nrisk_, nevent_, pntimes, wt, stat, var);

      *(pStatTend + nsim*j + ii) = *stat/pow(ntrial, 0.5);
      *(pVarTend + nsim*j + ii) = (*var)/ntrial;
      *(pmTend1 + nsim*j + ii) = (*m1)/ntrial;

      if(ii==nsim-1){
	for(l=0;l<ntimes;l++){
	  *(time+l) = *(time_ + l);
	  *(nrisk+l) = *(nrisk_+l);

	  *(nevent+l) = *(nevent_+l);
	}

        for(l=0;l<(*pntot);l++){
	  idx = (Yord+l)->index;
	  *(TOS+idx) = (Yord+l)->time;
	  *(Event+idx) = (Yord+l)->event;
	  *(Arm+idx) = (Yord+l)->arm;
        }
      }
 
      Free(time_);
      Free(nrisk_);
      Free(nevent_);
      Free(UQt);
      Free(varQt);
      Free(var1t);

      k=0;                         
      kacte =0;
      kactf =0;
      *dlact=0;
      *(dlact+1)=0;
      RejNull=0;
      AccNull=0;

      *pInfTold          = 0.0;
      *(pInfTold + 1)    = 0.0;
      *pInfTold_ii       = 0.0;
      *(pInfTold_ii + 1) = 0.0;

      if((*mybounds==0) || (nbnd_e_sv==3)){
        if((*nbnd==1) || (nbnd_e_sv==3)) *pbold = normut;   
        if(*nbnd==2) *pbold = *bHay;
      }
      if((*(mybounds+1)==0) || (nbnd_f_sv==3)){
        if((*(nbnd+1)==1) || (nbnd_f_sv==3)) *(pbold+1) = -normut;
        if(*(nbnd+1)==2) *(pbold+1) = *bHay;
      }

      *palphatot = atotsv_e;
      *(palphatot+1) = atotsv_f;

      *palpha = 0.0;
      *(palpha+1) = 0.0;
      for(l=0;l<2*ngqnodes;l++) {
        *(px+l) = 0.0;
        *(py+l) = 0.0;
        *(ptmp+l) = 0.0;
        *(pintgrndx+l) = 0.0;
      }

      *etaold = 0.0;
      krchd_flag=0;
      evnts_krchd=0.0;
      var_krchd=0.0;
      f_krchd = 0.0;
      f_krchd_ii = 0.0;

      isbad = 0;
      while((k<nlook) && (1-RejNull) && (1-AccNull) && (!isbad)){
        *ptlook = *(pttlook+k);
        tlook_ = *ptlook;
        *pnthslook = kacte+1;
	*(pnthslook+1) = kactf+1;

        handle(pn, ptlook, u, t0, t1, tc0, tc1, Yord, pntot, pntimesk);

	ntimesk = *pntimesk;
	time_   = Calloc(ntimesk, double);
	nrisk_  = Calloc(2*ntimesk, int);
	nevent_ = Calloc(2*ntimesk, int);
        UQt     = Calloc(ntimes, double);
        varQt   = Calloc(ntimes, double);
        var1t   = Calloc(ntimes, double);

        cpblocked(Yord, pntot, time_, nrisk_, nevent_, pntimesk, &pnevty, pnblocks);
 
        (*wtfun)(time_, nrisk_, nevent_, pntimesk, par, wt);
        if(*(stattype +j)==0) wlrstat(time_,nrisk_,nevent_,wt,pntimesk,statk,vark,m1k,UQt,varQt,var1t);
        if(*(stattype +j)==1) ISDstat(time_, nrisk_, nevent_, pntimesk, wt, statk, vark);

        wlrZ = *statk/pow(*vark, 0.5);
        *pInfTnew = *vark/(*var);
        totev_k = 0;
        for(i=0;i<2*ntimesk;i++) totev_k += *(nevent_ + i);

        if((k>= *spend_info_k) && (krchd_flag==0)) {
          var_krchd = *vark;
	  evnts_krchd =totev_k;
          krchd_flag=1;
	}

        /* use the Stochastic Curtailment procedure to construct the efficacy boundary */
        if(nbnd_e_sv==3){
          *pef = 0;
          *mufuforSC = *(mufu + nlook*j + k);
          *(mufuforSC + 1) = *(mufu + nlook*j + nlook - 1);
	  StCu2Bnds(mufuforSC, pInfTnew, bendSC, rho_sc, pef, pbounds+nlook*j+k);
          *mybounds = 1;
          *nbnd = 1;
        }
        /* use the Stochastic Curtailment procedure to constuct the futility boundary */
        if(nbnd_f_sv==3){
          *pef = 1;
          *mufuforSC = *(mufu + nlook*j + k);
          *(mufuforSC + 1) = *(mufu + nlook*j + nlook - 1);
	  StCu2Bnds(mufuforSC, pInfTnew, bendSC, rho_sc+1, pef, pbounds+nstat*nlook+nlook*j+k);
          *(mybounds + 1) = 1;
          *(nbnd + 1) = 1;
        }
        *pbnew = normut;
        *(pbnew+1) = -normut;
        if(*mybounds==1) *pbnew = *(pbounds + nlook*j + k);
        if(*(mybounds+1)==1) *(pbnew+1) =  *(pbounds + nstat*nlook + nlook*j + k);

        f_k = *pInfTnew;
	f_k_ii = ((double)totev_k)/((double)totev);
	f_krchd = var_krchd/(*var);
        f_krchd_ii = ((double)evnts_krchd)/((double)totev);

        /* Error probability spending information on Variance Scale */
        if(spend_info==0) 
	  *pInfTnew_ii = f_k;

        /* Error probability spending information on Events Scale */
	if(spend_info==1)
	  *pInfTnew_ii = f_k_ii; 

        /* Hybrid: linear switch from Variance Scale to Events Scale starting after analysis '*spend_info_k' */
        if(spend_info==2)
	{
          if(k>= *spend_info_k) 
	  {
	    *pInfTnew_ii = f_krchd +  (1.0 - f_krchd)/(1.0 - f_krchd_ii) * (f_k_ii - f_krchd_ii);
	  }
          else *pInfTnew_ii = f_k;
	}
	
        /* Error probability spending information on Calender Time Scale */
	if(spend_info==3)
	  *pInfTnew_ii = *ptlook/(*(pttlook + nlook -1));

	grpseqbnds(dofu,nbf,nbnd,nsf,rho,pnthslook,palphatot,palpha,psimin,dlact,
		   pInfTold,pInfTnew,pInfTold_ii,pInfTnew_ii,px,py,ptmp,pintgrndx,
		   pgqxw,pngqnodes,mufu + nlook*j + k,pbold,pbnew, mybounds);

	/* test using *pbnew */
	wlrsgn = (wlrZ >= 0.0 ? 1.0 : -1.0);
	if(abs(sided_) !=2)
	  RejNull = (wlrZ >= *pbnew)*(sided_ ==1) + (wlrZ <= ((*pbnew)*(-1.0)))*(sided_ ==-1);
	if(abs(sided_) ==2)
	  RejNull = (wlrZ * wlrsgn >= *pbnew);
	if(*dofu==1){
	  if(abs(sided_) !=2)
	    AccNull = (wlrZ <= *(pbnew+1))*(sided_ ==1) + (wlrZ >= ((*(pbnew+1))*(-1.0)))*(sided_ ==-1);
	  if(abs(sided_) ==2)
	    AccNull = (wlrZ * wlrsgn <= *(pbnew+1));
	}
/*	
        Rprintf("ii,k,nthslook,kacte,kactf,last,aold,anew,bold,bnew,tlook_,max(TOS),wlrZ,RejNull,AccNull\n");   
        Rprintf("%d,%d,%d,%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%d,%d\n",ii,k,*pnthslook,kacte,kactf,*(pbold+1),*(pbnew+1),
		*pbold,*pbnew,tlook_,*(time_+ *pntimesk - 1),wlrZ,RejNull,AccNull);   
*/
	if(*dlact==1){
	  if(*nbnd==1) {
            b_tmp = (*mybounds ? *(pbounds + nlook*j + k): *pbnew);
	    *(pbounds + nlook*j + k) = b_tmp;
	    *pbold = b_tmp;
	  }
          if(*nbnd==2){
	    if(1.0-*pInfTnew_ii >= 1e-6) {
	      *(pbounds + nlook*j + k) = *bHay;
	      *pbold = *bHay;
	      *palphatot = *palphatot - *palpha;
	    }
	    else *(pbounds + nlook*j + k) = *pbnew;
	  }
	  *(palphavec + nlook*j + k) = *palpha;
	  *pInfTold = *pInfTnew;
	  *pInfTold_ii = *pInfTnew_ii;
          kacte += 1;
	}
	else {
	  b_tmp = (*mybounds ? *(pbounds + nlook*j + k) : normut);
	  *(pbounds + nlook*j + k) = b_tmp; 
	  *pbold = b_tmp;
	  *(palphavec + nlook*j + k) = *psimin;
/*
	  *pInfTold = *pInfTnew;
*/
	  *pInfTold_ii = *pInfTnew_ii;

	}
	if((*dofu==1) && (*(dlact+1)==1)){
	  if(*(nbnd+1)==1) {
	    b_tmp = ((*(mybounds+1)) && (1.0 - *pInfTnew_ii >= 1e-6) ? *(pbounds + nstat*nlook + nlook*j + k) : *(pbnew+1));
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
	  *(palphavec + nstat*nlook + nlook*j + k) = *(palpha+1);
	  *(pInfTold + 1) = *pInfTnew;
	  *(pInfTold_ii + 1) = *pInfTnew_ii;
	  kactf += 1;
	}
	if((*dofu==1) && (*(dlact+1)==0)){
	  b_tmp = (*(mybounds+1) ? *(pbounds + nstat*nlook + nlook*j + k) : -normut);
	  *(pbounds + nstat*nlook + nlook*j + k) = b_tmp;
	  *(pbold+1) = b_tmp;
	  *(palphavec + nstat*nlook + nlook*j + k) = *psimin;

/*
	  *(pInfTold + 1) = *pInfTnew;
*/
	  *(pInfTold_ii +1) = *pInfTnew_ii;

	}
	*(pinffrac + nlook*j + k) = *pInfTnew;
	*(pinffrac_ii + nlook*j + k) = *pInfTnew_ii;
	*etaold = *etanew;

        isbad = isbad || R_isnancpp(*(pbounds + nlook*j + k));
        isbad = isbad || R_isnancpp(*(pbounds + nstat*nlook + nlook*j + k));
        isbad = isbad || R_isnancpp(*(pinffrac + nlook*j + k));
        isbad = isbad || R_isnancpp(*(pinffrac_ii + nlook*j + k));

        if(!isbad)
        {
            *(avgbounds + nlook*j + k) = *(avgbounds + nlook*j + k) + *(pbounds + nlook*j + k);
            *(avgbounds + nstat*nlook + nlook*j + k) = *(avgbounds + nstat*nlook + nlook*j + k) + 
                                                           *(pbounds + nstat*nlook + nlook*j + k);

            *(avginffrac + nlook*j + k) = *(avginffrac + nlook*j + k) + *(pinffrac + nlook*j + k);
            *(avginffrac_ii + nlook*j + k) = *(avginffrac_ii + nlook*j + k) + *(pinffrac_ii + nlook*j + k);
	}
	k++;
	Free(time_);
	Free(nrisk_);
	Free(nevent_);
	Free(UQt);
	Free(varQt);
	Free(var1t);
      } /* End sequential looks loop */
      if(RejNull && do_proj_)
      {
        *t_proj=*par;
        *(t_proj+1)=t_end-*accru;
        *(t_proj+2)=t_end;
        project_end(u, t0, t1, tc0, tc1, pn, t_proj, v_Tend_proj+nsim*j+ii, m_Tend_proj+nsim*j+ii);
      }
      *(pstatk+j) = *statk;
      *(pvark+j) = *vark;
      *(pstatK+nsim*j+ii) = *statk;
      *(pvarK+nsim*j+ii) = *vark;
      *(pRejAcc + nsim*j + ii) = RejNull;
      *(pRejAcc + nsim*nstat + nsim*j + ii) = AccNull;
      *(kstop + nsim*j + ii) = k;
      *(duration + nsim*j + ii) = tlook_;
      csumnppar += nppar;
    } /* END different stats loop */
    Rprintf("%g\r",round(1000*((double)ii)/((double)nsim))/100);
    if(!isbad) ii++;
  }
  /* END simulation loop */

  PutRNGstate();


  Free(Yord);

  Free(pinffrac);
  Free(pinffrac_ii);
  Free(pbounds);
  Free(pstatk);
  Free(pvark);
  Free(par);
  Free(wt);
  Free(px);
  Free(py);
  Free(ptmp);
  Free(pintgrndx);
  Free(statk);
  Free(vark);
  Free(m1k);
  Free(stat);
  Free(var);
  Free(m1);
  Free(etaold);
  Free(etanew);
  Free(ptlook);
  Free(pInfTold);
  Free(pInfTnew);
  Free(pInfTold_ii);
  Free(pInfTnew_ii);
  Free(palpha);
  Free(pbold);
  Free(pbnew);
  Free(psimin);
  Free(Vend);
  Free(mufuforSC);
  Free(UQt);
  Free(varQt);
  Free(var1t);
  Free(pntimesk);
  Free(pnblocks);
  Free(pnthslook);
  Free(pn);
  Free(pntot);
  Free(dlact);
  Free(puserVend);
  Free(pef);
}

void randfromh(int *pn, double *tcut, double *h, int *pncut, double *t)
{
  int n,ncut,i,l;
  double u,X,HI_;
  double *H;
/*-------------------------------------------------------------------------------------
     simulates variates from a peicewise exponential distribution.
    
     arguments
    
        pn: *int;   number of simulations desired
     pncut: *int;   length of vector specifications for changepoints and hazards
      tcut: *double; vector of left-hand endpoints of constant hazards interval; 
                               last right-hand endpoint (not specified) is infinity.
         h: *double; vector of hazard function values.
         t: *double; vector of length *pn containing the answer.  Must be allocated 
                               in calling routine.
---------------------------------------------------------------------------------------*/
  n = *pn;
  ncut = *pncut;
  H = Calloc(ncut,double);

  COMPH(tcut,h,H,ncut,l);

  for(i=0;i<n;i++){

    u = unif_rand();
    X = -log(u);

    HIatW(X,tcut,h,H,HI_,ncut,l);

    *(t+i) = HI_;
  }
  Free(H);
}

void randhcdtl(int *pn, double *tcut, double *h, int *pncut, double *tend,
	       double *tcutdA, double *hdA, int *pncutdA, double *tdA,
	       double *tcutdB, double *hdB, int *pncutdB, double *tdB,
               double *tcutxA, double *hxA, int *pncutxA, double *tcutxB, 
	       double *hxB, int *pncutxB, int *gradual, int *code, double *t)
{
  int n,ncut,ncutdA,ncutdB,ncutxA,ncutxB,nncutxA,nncutxB,ncutx,i,l,cd=0;
  int *pnncutxA,*pnncutxB;
  double tdAi, tdBi, td=0.0, u,X,htd,Htd,hx_td,Hx_td,HI_,htd_,Htd_;
  double hdAtend,HdAtend=0.0,hdBtend,HdBtend=0.0,pi;
  double *H,*HdA,*HdB,*HxA,*HxB,*tcutx,*hx,*Hx;
  double *ttcutxA,*hh_A,*HH_A,*hhxA,*HHxA,*ttcutxB,*hh_B,*HH_B,*hhxB,*HHxB;

/*-----------------------------------------------------------------------------------
     Specialized code.  The conditional hazard function (dropping dependance on 
     trial arm) is:
    
     gamma(x, tdA, tdB) = I(x <= tdA /\ tdB) h(x) + I(tdA < x /\ tdB) h_xA(x) 
                                                  + I(tdB < x /\ tdA) h_xB(x)
     The cummulative hazard becomes:
    
     Gamma(x, tdA, tdB) =   I(x <= tdA/\ tdB) H(x) 
                          + I(tdA < x /\ tdB)(H_xA(x) - H_xA(tdA) + H(tdA))
                          + I(tdB < x /\ tdA)(H_xB(x) - H_xB(tdB) + H(tdB))
    
     To simulate a variate, conditional upon 'tdA' and 'tdB', that has this 
     cumulative hazard function, this program starts with a unit exponential, 'X', 
     and finds the inverse of the above 'Gamma' at 'X'.
     This is done by (i)   find 'td' = 'tdA' /\ 'tdB' the minimum of 'tdA' and 'tdB'.
                     (ii)  find 'H(td)' and compare it to 'X'.
                           If 'H(td)' is larger than 'X' then find 'T' such that 
                           'H(T)' = 'X' and output 'T'
                           If 'H(td)' is less than 'X' then let 'DX' = 'X' - 'H(td)'.
                     (iii) let 'DX' = 'DX' + 'H_xL(td)'
                     (iv)  Find 'T' such that 'H_xL(T)' = 'DX' . 
                           In the above, 'L' = 'A' if td = tdA, 'L' = 'B' if td = tdB.  
                     (v)   Output T
-----------------------------------------------------------------------------------*/

  /* times to leave for reason A        */
  randfromh(pn, tcutdA, hdA, pncutdA, tdA);

  /* times to leave for reason B        */
  randfromh(pn, tcutdB, hdB, pncutdB, tdB);

  n = *pn;
  ncut = *pncut;
  ncutdA = *pncutdA;
  ncutdB = *pncutdB;
  ncutxA = *pncutxA;
  ncutxB = *pncutxB;
  nncutxA = MAX(ncut,ncutxA);
  nncutxB = MAX(ncut,ncutxB);
  ncutx = MAX(nncutxA,nncutxB);

  H =       Calloc(ncut,double);
  HxA =     Calloc(ncutxA,double);
  HxB =     Calloc(ncutxB,double);
  //  if(*gradual==1){
  HdA =     Calloc(ncutdA,double);
  HdB =     Calloc(ncutdB,double);
  ttcutxA = Calloc(nncutxA,double);
  hh_A =    Calloc(nncutxA,double);
  HH_A =    Calloc(nncutxA,double);
  hhxA =    Calloc(nncutxA,double);
  HHxA =    Calloc(nncutxA,double);
  ttcutxB = Calloc(nncutxB,double);
  hh_B =    Calloc(nncutxB,double);
  HH_B =    Calloc(nncutxB,double);
  hhxB =    Calloc(nncutxB,double);
  HHxB =    Calloc(nncutxB,double);
  tcutx =   Calloc(ncutx,double);
  hx =      Calloc(ncutx,double);
  Hx =      Calloc(ncutx,double);
  pnncutxA =Calloc(1,int);
  pnncutxB =Calloc(1,int);
  //  }

  COMPH(tcut,h,H,ncut,l);
  COMPH(tcutxA,hxA,HxA,ncutxA,l);
  COMPH(tcutxB,hxB,HxB,ncutxB,l);

  if(*gradual == 1){
    commonx(tcut,h,pncut,tcutxA,hxA,pncutxA,ttcutxA,hh_A,hhxA,pnncutxA);
    nncutxA = *pnncutxA;
    commonx(tcut,h,pncut,tcutxB,hxB,pncutxB,ttcutxB,hh_B,hhxB,pnncutxB);
    nncutxB = *pnncutxB;

    COMPH(tcutdA,hdA,HdA,ncutdA,l);
    COMPH(tcutdB,hdB,HdB,ncutdB,l);
    COMPH(ttcutxA,hh_A,HH_A,nncutxA,l);
    COMPH(ttcutxA,hhxA,HHxA,nncutxA,l);
    COMPH(ttcutxB,hh_B,HH_B,nncutxB,l);
    COMPH(ttcutxB,hhxB,HHxB,nncutxB,l);

    hHatX(*tend,tcutdA,hdA,HdA,hdAtend,HdAtend,ncutdA,l);
    hHatX(*tend,tcutdB,hdB,HdB,hdBtend,HdBtend,ncutdB,l);
  }

  for(i=0;i<n;i++){
    *(code + i) = 0;
    tdAi = *(tdA + i);
    tdBi = *(tdB + i);
    if(tdAi<=tdBi){
      td = tdAi;
      if(*gradual == 1){
        hHatX(td,tcutdA,hdA,HdA,htd_,Htd_,ncutdA,l);
	pi = 1.0;
	if(td <= *tend) pi = Htd_/HdAtend;
        for(l=0;l<nncutxA;l++) {
	  *(tcutx+l) = *(ttcutxA+l);
	  *(hx+l) = pi* (*(hh_A+l)) + (1.0-pi)* (*(hhxA+l));
	  *(Hx+l) = pi* (*(HH_A+l)) + (1.0-pi)* (*(HHxA+l));
	  ncutx = nncutxA;
        }
      }
      else{
	tcutx = tcutxA;
	hx = hxA;
	Hx = HxA;
	ncutx = *pncutxA;
      }
      cd = 1;
    }
    if(tdAi>tdBi){
      td = tdBi;
      if(*gradual == 1){
        hHatX(td,tcutdB,hdB,HdB,htd_,Htd_,ncutdB,l);
	pi = 1.0;
        if(td <= *tend) pi = Htd_/HdBtend;
        for(l=0;l<nncutxB;l++) {
	  *(tcutx+l) = *(ttcutxB+l);
	  *(hx+l) = pi* (*(hh_B+l)) + (1.0-pi)* (*(hhxB+l));
	  *(Hx+l) = pi* (*(HH_B+l)) + (1.0-pi)* (*(HHxB+l));
	  ncutx = nncutxB;
	}
      }
      else{
	tcutx = tcutxB;
	hx = hxB;
	Hx = HxB;
	ncutx = *pncutxB;
      }
      cd = 2;
    }

    hHatX(td,tcut,h,H,htd,Htd,ncut,l);
    hHatX(td,tcutx,hx,Hx,hx_td,Hx_td,ncutx,l);
    u = unif_rand();
    X = -log(u);
    if(Htd >= X){
      HIatW(X,tcut,h,H,HI_,ncut,l);
      *(t+i) = HI_;
    }
    if(Htd < X){
      *(code + i) = cd;
      X = X - Htd + Hx_td;
      HIatW(X,tcutx,hx,Hx,HI_,ncutx,l);
      *(t+i) = HI_;
    }
  }
  Free(H);
  Free(HxA);
  Free(HxB);
  if(*gradual == 1){
    Free(HdA);
    Free(HdB);
    Free(ttcutxA);
    Free(hh_A);
    Free(HH_A);
    Free(hhxA);
    Free(HHxA);
    Free(ttcutxB);
    Free(hh_B);
    Free(HH_B);
    Free(hhxB);
    Free(HHxB);
    Free(tcutx);
    Free(hx);
    Free(Hx);
    Free(pnncutxA);
    Free(pnncutxB);
  }
}

void handle(int *pn, double *tlook, double *u, double *t0, double *t1, double *tc0,
            double *tc1, itea *YY, int *pntot, int *pntimes)
{
  int n, n0, n1, ntot, i, j, ndths;
  double tout, tmpp;

  n = *pn;

  ndths = 0;
  j=0;
  for(i=0;i<n;i++){
    tout = (*tlook - *(u+i) > 0.0 ? *tlook - *(u+i) : 0.0);
    if(*(u+i)<*tlook){
      tmpp = (*(tc0+i) < tout ? *(tc0+i) : tout);
      (YY+j)->index = j;
      (YY+j)->time = (*(t0+i) <= tmpp ? *(t0+i) : tmpp);
      (YY+j)->event = (*(t0+i) <= tmpp ? 1 : 0);
      ndths += (YY+j)->event;
      (YY+j)->arm = 0;
      j++;
    }
  }
  n0 = j;
  j=0;
  for(i=0;i<n;i++){
    tout = (*tlook - *(u+n+i) > 0.0 ? *tlook - *(u+n+i) : 0.0);
    if(*(u+n+i)<*tlook){
      tmpp = (*(tc1+i) < tout ? *(tc1+i) : tout);
      (YY+n0+j)->index = n0 + j;
      (YY+n0+j)->time = (*(t1+i) <= tmpp ? *(t1+i) : tmpp);
      (YY+n0+j)->event = (*(t1+i) <= tmpp ? 1 : 0);
      ndths += (YY+n0+j)->event;
      (YY+n0+j)->arm = 1;
      j++;
    }
  }
  n1 = j;
  ntot = n0+n1;
  *pntot = ntot;
  *pntimes = ndths;
}

void commonx(double *x1,double *h1,int *pn1,double *x2,double *h2,int *pn2,
	     double *x,double *hh1,double *hh2,int *pn)
{
  int l,j,n1,n2,n,nx;
  double yold,ynew,h1_,h2_;
  double *y;
  n1 = *pn1;
  n2 = *pn2;
  n=n1+n2;
  y = Calloc(n,double);
  for(l=0;l<n1;l++) *(y+l) = *(x1+l);
  for(l=0;l<n2;l++) *(y+n1+l) = *(x2+l);
  cmprdbl = &CmprDbl;
  qsort(y, n, sizeof(double), cmprdbl);
  j=1;
  yold = *y;
  *x = yold;
  for(l=1;l<n;l++){
    ynew = *(y+l);
    if(ynew!=yold) {
      *(x+j) = ynew;
      j++;
    }
    yold = ynew;
  }
  nx = j;
  for(j=0;j<nx;j++){
    hatX(*(x+j),x1,h1,h1_,n1,l);
    hatX(*(x+j),x2,h2,h2_,n2,l);
    *(hh1+j) = h1_;
    *(hh2+j) = h2_;
  }
  *pn = nx;  
  Free(y);
}

int CmprDbl(const void *x, const void *y)
{
  double *xx, *yy;
  xx = (double *) x;
  yy = (double *) y;
  return(1*(*xx > *yy) - 1*(*xx < *yy));
}

void printbtre(double *time, int *nrisk, int *nevent, int nt, int fromi, int toi, int nb)
{
    int i,j,n;
    n = toi - fromi + 1;
    for(i=0;i<n;i++){
	Rprintf("[%d] %g    ",fromi+i, *(time +fromi-1 +i));
	for(j=0;j<nb;j++) Rprintf("%d   ", *(nrisk + nb*(fromi-1 +i) + j));
	for(j=0;j<nb;j++) Rprintf("%d   ", *(nevent +nb*(fromi-1 +i) + j));
	Rprintf("\n");
    }
}

/*--------------------------------------------------------------------------------------------
 function definitions for

   void wlrstat(double *time, int *nrisk, int *nevent, int *pntimes, double *wt, double *stat,
                double *var, double *eta)
   void flemhar(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt)
   void sflemhar(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt)
   void ramp(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt)

   are contained in the file WtdLogRank.c 
----------------------------------------------------------------------------------------------*/



