#include<stdlib.h>
#include<R.h>
#include "PwrGSD_mem.h"
//
//MACROS 
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define COMPH(xh,h,H,n,l) *H=0.0;for(l=1;l<n;l++) *(H+l)=*(H+l-1)+ *(h+l-1) * (*(xh+l) - *(xh+l-1))
#define hHatX(x,xh,h,H,h_,H_,n,l,flg) l=0; flg=1; while(flg){                                                   \
                                                    flg=1*(l<n); if(flg) flg=1*(*(xh+l)<=x); l+=flg;            \
                                                  }                                                             \
                                                  h_ = *(h+l-1); H_ = *(H+l-1) + h_ * (x-*(xh+l-1))
#define HIatW(W,xh,h,H,HI_,n,l,flg) l=0; flg=1; while(flg){                                                     \
                                                  flg=1*(l<n); if(flg) flg=1*(*(H+l)<=W); l+=flg;               \
                                                }                                                               \
                                                HI_ = *(xh+l-1) + (W - *(H+l-1))/(*(h+l-1))
#define wchidx(x,xgrid,n,l,flg) l=0; flg=1; while(flg){                                                         \
                                                flg=1*(l<n); if(flg) flg=1*((*(xgrid+l)<=x) && (l<n)); l+=flg;	\
                                            } /*   don't remove this line or add anything onto the end please  */
#define unique(x,n,y,m,l,old) m = 0;                           \
                              l = 0;                           \
                              old = -10 + *x;                  \
                              while(l<n){		       \
                                if(fabs(*(x+l)-old)<=1e-16) {  \
                                  old = *(x+l);                \
                                  l++;                         \
                                }                              \
                                else {                         \
                                  *(y+m) = *(x+l);             \
                                  old = *(x+l);                \
                                  l++;                         \
                                  m++;                         \
                                }                              \
			      }   /*   don't remove this line or add anything onto the end please */

//
//

void printmat(double *,int,int,char *);
void printmati(int *, int,int,char *);

void htilde(double *x,int *nx,double *gqx,double *gqw,int *ngq,
            double *xh,double *h,int *nh,double *xhA,double *hA,int *nhA,
            double *xhB,double *hB,int *nhB,double *xlA,double *lA,int *nlA,
            double *xlB,double *lB,int *nlB,int *gradual,double *tend,
	    double *ftlde,double *Stlde,double *htlde);

typedef int cmprfun(const void *, const void *);
cmprfun cmpr;

void drift(int *ints,double *accru,double *accrat,double *tlook,double *ppar,double *gqxw,
	   double *th0,double *h0,double *th1,double *h1,double *thc0,double *hc0,
	   double *thc1,double *hc1,double *tlA0,double *lA0,double *tlB0,double *lB0,
	   double *tlA1,double *lA1,double *tlB1,double *lB1,double *thA0,double *hA0,
	   double *thB0,double *hB0,double *thA1,double *hA1,double *thB1,double *hB1,
	   int *wttyp,double *RR,int *pnnjmp,double *InfFrac,double *InfFrac_ii,
	   double *mu,double *betastar, double *Var_uw, double *Var, double *Eta,
           int *puserVend,double *Vend)
{
  int nnstat,nngq,l,flg,j,j0,j1,istat,nnh0,nnh1,nnsm,nnjmp,nnhc0,nnhc1,nmx,flaguserVE;
  int nnlook,ilook,ijmp,ntrial,isumppar,nppar=0;
  int *ntlook,*nstat,*ngq,*nh0,*nh1,*nhc0,*nhc1,*nlA0,*nlB0,*nlA1,*nlB1,*nhA0,*nhB0,*nhA1;
  int *nhB1,*gradual,*one;

  double t_ENR,t_END,vend_uw,vend,xi_,w,Q,dV_uw,dV,dMU,Beta,ans_MU=0.0,hc0_,Hc0_,hc1_;
  double Hc1_,Sc,S_LR,Stlde,ftlde,V_uw,V,et=0.0, old, t_jmp,ans_IF_uw, ans_IF;
  double *xi,*ftlde0,*Stlde0,*htlde0,*ftlde1,*Stlde1,*htlde1,*Hc0,*Hc1,*V_END_uw,*V_END;
  double *tjump,*tjmp,*tend,*gqx,*gqw,*Qstop,*fftlde0,*SStlde0,*hhtlde0,*fftlde1,*SStlde1,*hhtlde1;

  ntlook  = ints;
  nstat   = ints +  1;
  ngq     = ints +  2;
  nh0     = ints +  3;
  nh1     = ints +  4; 
  nhc0    = ints +  5;
  nhc1    = ints +  6;
  nlA0    = ints +  7;
  nlB0    = ints +  8; 
  nlA1    = ints +  9;
  nlB1    = ints + 10; 
  nhA0    = ints + 11;
  nhB0    = ints + 12; 
  nhA1    = ints + 13; 
  nhB1    = ints + 14;
  gradual = ints + 15;

  ntrial = (int)(*accru * *accrat);
  ntrial += (ntrial % 2);

  nnlook = *ntlook;
  nnstat = *nstat;
  t_ENR = *accru;
  t_END = *(tlook+nnlook-1);
  nngq = *ngq;
  gqx = gqxw;
  gqw = gqxw+nngq;
  nnh0 = *nh0;
  nnh1 = *nh1;
  nnsm = nnh0 + nnh1 + nnlook - 2;
  nnhc0 = *nhc0;
  nnhc1 = *nhc1;

  Hc0      = Calloc(nnhc0, double);
  Hc1      = Calloc(nnhc1, double);
  V_END_uw = Calloc(nnstat, double);
  V_END    = Calloc(nnstat, double);

  tjump    = Calloc(nnsm, double);
  tjmp     = Calloc(nnsm, double);
  tend     = Calloc(1, double);
  xi       = Calloc(nngq, double);
  Qstop    = Calloc(nnstat, double);
  one      = Calloc(1, int);
  *one = 1;
  Q=0.0;

  if(nnh0>1) for(l=0;l<nnh0-1;l++) *(tjump + l) = *(th0 + l + 1);
  if(nnh1>1) for(l=0;l<nnh1-1;l++) *(tjump + nnh0 - 1 + l) = *(th1 + l + 1);
  for(l=0;l<nnlook;l++) *(tjump + nnh0 + nnh1 - 2 + l) = *(tlook + l);

  qsort(tjump, nnsm, sizeof(double), &cmpr);
  unique(tjump, nnsm, tjmp, nnjmp, l, old);
  *pnnjmp = nnjmp;
  *tend = t_END;

  nmx = MAX(nngq, nnjmp);

  ftlde0 =  Calloc(nmx, double);
  Stlde0 =  Calloc(nmx, double);
  htlde0 =  Calloc(nmx, double);
  ftlde1 =  Calloc(nmx, double);
  Stlde1 =  Calloc(nmx, double);
  htlde1 =  Calloc(nmx, double);

  fftlde0 = Calloc(1, double);
  SStlde0 = Calloc(1, double);
  hhtlde0 = Calloc(1, double);
  fftlde1 = Calloc(1, double);
  SStlde1 = Calloc(1, double);
  hhtlde1 = Calloc(1, double);

  *ftlde0=0.0;
  *Stlde0=1.0;
  *htlde0=0.0;
  *ftlde1=0.0;
  *Stlde1=1.0;
  *htlde1=0.0;

  htilde(tjmp,pnnjmp,gqx,gqw,ngq,th0,h0,nh0,thA0,hA0,nhA0,thB0,hB0,nhB0,tlA0,lA0,
	 nlA0,tlB0,lB0,nlB0,gradual,tend,ftlde0,Stlde0,htlde0);

  htilde(tjmp,pnnjmp,gqx,gqw,ngq,th1,h1,nh1,thA1,hA1,nhA1,thB1,hB1,nhB1,tlA1,lA1,
	 nlA1,tlB1,lB1,nlB1,gradual,tend,ftlde1,Stlde1,htlde1);

  for(l=0;l<nnjmp;l++){
    *(RR+l) = *(tjmp+l);
    *(RR+nnjmp+l) = *(htlde0+l);
    *(RR+2*nnjmp+l) = *(htlde1+l)/(*(htlde0+l));
    wchidx(*(tjmp+l),th0,nnh0,j0,flg);
    j0 = MAX(j0-1,0);
    *(RR+3*nnjmp+l) = *(h0+j0);
    wchidx(*(tjmp+l),th1,nnh1,j1,flg);
    j1 = MAX(j1-1,0);
    *(RR+4*nnjmp+l) = *(h1+j1)/(*(h0+j0));
  }

  for(l=0;l<nngq;l++) *(xi+l) = (*(gqx+l)+1.0)*t_END/2.0;

  htilde(xi,ngq,gqx,gqw,ngq,th0,h0,nh0,thA0,hA0,nhA0,thB0,hB0,nhB0,tlA0,lA0,
	 nlA0,tlB0,lB0,nlB0,gradual,tend,ftlde0,Stlde0,htlde0);

  htilde(xi,ngq,gqx,gqw,ngq,th1,h1,nh1,thA1,hA1,nhA1,thB1,hB1,nhB1,tlA1,lA1,
	 nlA1,tlB1,lB1,nlB1,gradual,tend,ftlde1,Stlde1,htlde1);

  COMPH(thc0,hc0,Hc0,nnhc0,l);
  COMPH(thc1,hc1,Hc1,nnhc1,l);

  isumppar = 0;
  flaguserVE = *puserVend;
  if(flaguserVE==0){
    for(istat=0;istat<nnstat;istat++){
      vend = 0.0;
      vend_uw = 0.0;
      if(*(wttyp+istat) == 1){
	htilde(ppar+isumppar+2,one,gqx,gqw,ngq,th0,h0,nh0,thA0,hA0,nhA0,thB0,hB0,nhB0,tlA0,lA0,
               nlA0,tlB0,lB0,nlB0,gradual,tend,fftlde0,SStlde0,hhtlde0);

        htilde(ppar+isumppar+2,one,gqx,gqw,ngq,th1,h1,nh1,thA1,hA1,nhA1,thB1,hB1,nhB1,tlA1,lA1,
               nlA1,tlB1,lB1,nlB1,gradual,tend,fftlde1,SStlde1,hhtlde1);

        Stlde = (*SStlde0 + *SStlde1)/2.0;
        *(Qstop + istat) = pow(Stlde,*(ppar+isumppar))*pow(1.0-Stlde,*(ppar+isumppar+1));
      }

      for(j=0;j<nngq;j++){
        xi_ = *(xi+j);
        w = *(gqw+j)*t_END/2.0;
        hHatX(xi_,thc0,hc0,Hc0,hc0_,Hc0_,nnhc0,l,flg);
        hHatX(xi_,thc1,hc1,Hc1,hc1_,Hc1_,nnhc1,l,flg);
        Sc = 2.0*exp(-(Hc0_ + Hc1_))/(exp(-Hc0_) + exp(-Hc1_));
        S_LR = MIN((t_END - xi_)/t_ENR,1.0);
        Stlde = (*(Stlde0+j) + *(Stlde1+j))/2.0;
        ftlde = (*(htlde0+j) * *(Stlde0+j) + *(htlde1+j) * *(Stlde1+j))/2.0;

        if(*(wttyp+istat) == 0){
          nppar = 2;
          Q = pow(Stlde,*(ppar+isumppar))*pow(1.0-Stlde,*(ppar+isumppar+1));
        }
        if(*(wttyp+istat) == 1){
          nppar = 3;
          Q = pow(Stlde,*(ppar+isumppar))*pow(1.0-Stlde,*(ppar+isumppar+1));
	  if(xi_ > *(ppar + isumppar+2)) {
            Q= *(Qstop+istat);
          }
        }
        if(*(wttyp+istat) == 2){
          nppar = 1;
          Q =MIN(xi_/(*(ppar+isumppar)),1.0);
        }
	vend_uw += 1.0/4.0 * Sc * S_LR * ftlde * w;
        vend += 1.0/4.0 * Q * Q * Sc * S_LR * ftlde * w;
      }
      *(V_END_uw + istat) = vend_uw;
      *(V_END + istat) = vend;
      isumppar += nppar;
    }
  }
  else{
    for(istat=0;istat<nnstat;istat++){
      *(V_END_uw + istat) = *(Vend + istat);
      *(V_END + istat) = *(Vend + istat);
    }
  }

  isumppar = 0;
  for(istat=0;istat<nnstat;istat++){
    ilook = 0;
    for(ijmp=0;ijmp<nnjmp;ijmp++){
      t_jmp = *(tjmp+ijmp);
      for(l=0;l<nngq;l++) *(xi+l) = (*(gqx+l)+1.0)*t_jmp/2.0;

      htilde(xi,ngq,gqx,gqw,ngq,th0,h0,nh0,thA0,hA0,nhA0,thB0,hB0,nhB0,tlA0,lA0,
	     nlA0,tlB0,lB0,nlB0,gradual,tend,ftlde0,Stlde0,htlde0);

      htilde(xi,ngq,gqx,gqw,ngq,th1,h1,nh1,thA1,hA1,nhA1,thB1,hB1,nhB1,tlA1,lA1,
	     nlA1,tlB1,lB1,nlB1,gradual,tend,ftlde1,Stlde1,htlde1);

      ans_IF_uw = 0.0;
      ans_IF = 0.0;
      ans_MU = 0.0;
      et = 0.0;
      V_uw = 0.0;
      V = 0.0;
      for(j=0;j<nngq;j++){
	xi_ = *(xi+j);
	w = *(gqw+j)*t_jmp/2.0;
	hHatX(xi_,thc0,hc0,Hc0,hc0_,Hc0_,nnhc0,l,flg);
	hHatX(xi_,thc1,hc1,Hc1,hc1_,Hc1_,nnhc1,l,flg);
	Sc = 2.0*exp(-(Hc0_ + Hc1_))/(exp(-Hc0_) + exp(-Hc1_));
	S_LR = MIN((t_jmp - xi_)/t_ENR,1.0);
	Stlde = (*(Stlde0+j) + *(Stlde1+j))/2.0;
	ftlde = (*(ftlde0+j) + *(ftlde1+j))/2.0;

	if(*(wttyp+istat) == 0){
	    nppar = 2;
	    Q = pow(Stlde,*(ppar+isumppar))*pow(1.0-Stlde,*(ppar+isumppar+1));
	}
	if(*(wttyp+istat) == 1){
	    nppar = 3;
	    Q = pow(Stlde,*(ppar+isumppar))*pow(1.0-Stlde,*(ppar+isumppar+1));
	    if(xi_ > *(ppar + isumppar+2)){
		Q= *(Qstop+istat);
	    }
	}
	if(*(wttyp+istat) == 2){
	    nppar = 1;
	    Q =MIN(xi_/(*(ppar+isumppar)),1.0);
	}

	dV_uw = Sc * S_LR * ftlde * w/4.0;
	dMU = Q * dV_uw;
	dV = Q * dMU;
	et += dMU;
	V_uw += dV_uw;
	V += dV;
	ans_IF_uw += dV_uw/(*(V_END_uw + istat));
	ans_IF += dV/(*(V_END + istat));
	Beta = log(*(htlde1+j)) - log(*(htlde0+j));
	ans_MU += Beta * dMU;
      }
      if((fabs(t_jmp - *(tlook+ilook))<=1e-16) && (ilook < nnlook)){
        *(InfFrac + nnlook*istat + ilook) = ans_IF;
	*(InfFrac_ii + nnlook*istat + ilook) = ans_IF_uw;
        *(mu + nnlook*istat + ilook) = pow(ntrial,0.5) * ans_MU/pow(*(V_END+istat),0.5);
        ilook++;
      }
      *(Var_uw + nnjmp*istat + ijmp) = V_uw;
      *(Var + nnjmp*istat + ijmp) = V;
      *(Eta + nnjmp*istat + ijmp) = et;
    }
    *(betastar + istat) = ans_MU/et;
    isumppar += nppar;
  }

  Free(Hc0);
  Free(Hc1);
  Free(V_END_uw);
  Free(V_END);
  Free(tjump);
  Free(tjmp);
  Free(tend);
  Free(xi);
  Free(Qstop);
  Free(one);
  Free(ftlde0);
  Free(Stlde0);
  Free(htlde0);
  Free(ftlde1);
  Free(Stlde1);
  Free(htlde1);
  Free(fftlde0);
  Free(SStlde0);
  Free(hhtlde0);
  Free(fftlde1);
  Free(SStlde1);
  Free(hhtlde1);
}

int cmpr(const void *x, const void *y)
{
  double *xx, *yy;
  xx = (double *) x;
  yy = (double *) y;
  return(1*((*xx > *yy) - (*xx < *yy)));
}
