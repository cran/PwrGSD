#include<R.h>

/* MACROS */
#define COMPH(xh,h,H,n,l) *H=0.0;for(l=1;l<n;l++) *(H+l)=*(H+l-1)+ *(h+l-1) * (*(xh+l) - *(xh+l-1))
#define hatX(x,xh,h,h_,n,l) l=0;while(*(xh+l)<=x && l<n) l++; h_ = *(h+l-1)
#define hHatX(x,xh,h,H,h_,H_,n,l) l=0;while(*(xh+l)<=x && l<n) l++; h_ = *(h+l-1); H_ = *(H+l-1) + h_ * (x-*(xh+l-1))
#define HIatW(W,xh,h,H,HI_,n,l) l=0;while(*(H+l)<=W && l<n) l++; HI_ = *(xh+l-1) + (W - *(H+l-1))/(*(h+l-1))
#define MIN(x,y) (x<=y ? x : y)
#define MAX(x,y) (x>=y ? x : y)

void printmat(double *,int,int,char *);
void Qmoments(double *pK, double *ph, double *ptc, double *ptr, double *ptau, double *ans);

void driftfu(int *ints,double *accru,double *accrat,double *tlook,double *ppar,double *gqxw,
	     double *th0,double *h0,double *lrrf,double *thc0,double *hc0,double *thc1,
	     double *hc1,int *wttyp,double *mufu,int *puserVend,double *Vend)
{
  int nnstat,nngq,l,j,istat,nnh0,nnh1,nnhc0,nnhc1,isumppar,istop,flaguserVE,qis1orQ;
  int nnlook,ilook,ntrial,nppar, *lbuff,*ntlook,*nstat,*ngq,*nh0,*nhc0,*nhc1,*pqis1orQ;

  double t_ENR,t_END,vend,xi_,w,t_AN,Q,dMU,Beta,ans_MU,lrrf_,h0_,H0_,S0_,h1_,H1_,S1_,hc0_,Hc0_,hc1_;
  double Hc1_,Sc,S_LR,S,f,h_avg;
  double *xi,*H0,*H1,*th1,*h1,*Hc0,*Hc1,*V_END,*gqx,*gqw,*Qstop,*atten,*Qmom_args;

  ntlook   = ints;
  nnlook = *ntlook;
  nstat    = ints + 1;
  ngq      = ints + 2;
  nh0      = ints + 3;
  nhc0     = ints + 5;
  nhc1     = ints + 6;
  pqis1orQ  = ints + 16+4*nnlook+7;

  ntrial = *accru * *accrat;
  ntrial += (ntrial % 2);

  nnstat = *nstat;
  t_ENR = *accru;
  t_END = *(tlook+nnlook-1);
  nngq = *ngq;
  gqx = gqxw;
  gqw = gqxw + nngq;
  nnh0 = *nh0;
  nnh1 = *nh0;
  nnhc0 = *nhc0;
  nnhc1 = *nhc1;
  qis1orQ = *pqis1orQ;

  xi = (double *)Calloc(nngq, double);
  H0 = (double *)Calloc(nnh0, double);
  th1 = (double *)Calloc(nnh1, double);
  h1 = (double *)Calloc(nnh1, double);
  H1 = (double *)Calloc(nnh1, double);
  Hc0 = (double *)Calloc(nnhc0, double);
  Hc1 = (double *)Calloc(nnhc1, double);
  lbuff = (int *)Calloc(1, int);
  V_END = (double *)Calloc(nnstat, double);
  Qstop = (double *)Calloc(nnstat, double);
  atten = (double *)Calloc(nnstat, double);
  Qmom_args = (double *)Calloc(7, double);

  *lbuff = nnh0;

  for(l=0;l<nnh1;l++) {
    *(th1+l) = *(th0+l);
    *(h1+l) = exp(*(lrrf + l)) * *(h0 + l);  
  }

  for(l=0;l<nngq;l++) *(xi+l) = (*(gqx+l)+1.0)*t_END/2.0;
  COMPH(th0,h0,H0,nnh0,l);
  COMPH(th1,h1,H1,nnh1,l);
  COMPH(thc0,hc0,Hc0,nnhc0,l);
  COMPH(thc1,hc1,Hc1,nnhc1,l);


  hHatX(t_END,th0,h0,H0,h0_,H0_,nnh0,l);
  hHatX(t_END,th1,h1,H1,h1_,H1_,nnh1,l);
  S0_ = exp(-H0_);
  S1_ = exp(-H1_);
  S = (S0_ + S1_)/2.0;
  h_avg = -log(S)/t_END;

  isumppar = 0;
  flaguserVE = *puserVend;
  if(flaguserVE==0){
    for(istat=0;istat<nnstat;istat++){
      *(atten + istat)=1.0;
      vend = 0.0;

      istop = 0.0;
      if(*(wttyp+istat) == 1){
        hHatX(*(ppar+isumppar+2),th0,h0,H0,h0_,H0_,nnh0,l);
        hHatX(*(ppar+isumppar+2),th1,h1,H1,h1_,H1_,nnh1,l);
        S0_ = exp(-H0_);
        S1_ = exp(-H1_);
        S = (S0_ + S1_)/2.0;
        *(Qstop + istat) = pow(S,*(ppar+isumppar))*pow(1.0-S,*(ppar+isumppar+1));
        if(qis1orQ==1){
	  /*  pK   <-> Qmom_arg    :  proportionality between other cause mortality and cancer mortality--use 30
	      ph   <-> Qmom_arg + 1:  constant cancer mortality rate to use
              ptc  <-> Qmom_arg + 2:  time at which weight function is capped
              ptr  <-> Qmom_arg + 3:  time at which randomization ends
              ptau <-> Qmom_arg + 4:  time at which trial ends
              ans  <-> Qmom_arg + 5:  (length 2) first and 2nd moments of Q.
          Qmoments(double *pK, double *ph, double *ptc, double *ptr, double *ptau, double *ans); */
          *Qmom_args       = 30.0;
          *(Qmom_args + 1) = h_avg;
          *(Qmom_args + 2) = *(ppar + isumppar+2);
          *(Qmom_args + 3) = t_ENR;
          *(Qmom_args + 4) = t_END;
	  Qmoments(Qmom_args, Qmom_args+1, Qmom_args+2, Qmom_args+3, Qmom_args+4, Qmom_args+5);
          *(atten + istat) = *(Qmom_args + 5)/(*(Qmom_args + 6));
        }
      }

      for(j=0;j<nngq;j++){
        xi_ = *(xi+j);
        w = *(gqw+j)*t_END/2.0;
        hHatX(xi_,th0,h0,H0,h0_,H0_,nnh0,l);
        hHatX(xi_,th1,h1,H1,h1_,H1_,nnh1,l);
        S0_ = exp(-H0_);
        S1_ = exp(-H1_);
        hHatX(xi_,thc0,hc0,Hc0,hc0_,Hc0_,nnhc0,l);
        hHatX(xi_,thc1,hc1,Hc1,hc1_,Hc1_,nnhc1,l);
        Sc = 2.0*exp(-(Hc0_ + Hc1_))/(exp(-Hc0_) + exp(-Hc1_));
        S_LR = MIN((t_END - xi_)/t_ENR,1.0);
        S = (S0_ + S1_)/2.0;
        f = (h0_ * S0_ + h1_ * S1_)/2.0;
        if(*(wttyp+istat) == 0){
          nppar = 2;
          Q = pow(S,*(ppar+isumppar))*pow(1.0-S,*(ppar+isumppar+1));
        }
        if(*(wttyp+istat) == 1){
          nppar = 3;
          Q = pow(S,*(ppar+isumppar))*pow(1.0-S,*(ppar+isumppar+1));
          if(xi_ > *(ppar + isumppar+2)) {
	    Q= *(Qstop + istat);
            istop=1;
	  }
        }
        if(*(wttyp+istat) == 2){
          nppar = 1;
          Q =MIN(xi_/(*(ppar+isumppar)),1.0);
        }
        vend += 1.0/4.0 * Q * Q * Sc * S_LR * f * w;
      }
      *(V_END + istat) = vend;
      isumppar += nppar;
    }
  }
  else{
    for(istat=0;istat<nnstat;istat++) *(V_END + istat) = *(Vend + istat);
  }

  isumppar = 0;
  for(istat=0;istat<nnstat;istat++){
    for(ilook=0;ilook<nnlook;ilook++){
      t_AN = *(tlook+ilook);
      for(l=0;l<nngq;l++) *(xi+l) = (*(gqx+l)+1.0)*t_AN/2.0;

      ans_MU = 0.0;
      istop = 0;
      for(j=0;j<nngq;j++){
	xi_ = *(xi+j);
	w = *(gqw+j)*t_AN/2.0;
        hHatX(xi_,th0,h0,H0,h0_,H0_,nnh0,l);
        hHatX(xi_,th1,h1,H1,h1_,H1_,nnh1,l);
	S0_ = exp(-H0_);
        S1_ = exp(-H1_);
	hHatX(xi_,thc0,hc0,Hc0,hc0_,Hc0_,nnhc0,l);
	hHatX(xi_,thc1,hc1,Hc1,hc1_,Hc1_,nnhc1,l);
	Sc = 2.0*exp(-(Hc0_ + Hc1_))/(exp(-Hc0_) + exp(-Hc1_));
	S_LR = MIN((t_AN - xi_)/t_ENR,1.0);
	S = (S0_ + S1_)/2.0;
	f = (h0_ * S0_ + h1_ * S1_)/2.0;

	if(*(wttyp+istat) == 0){
	    nppar = 2;
	    Q = pow(S,*(ppar+isumppar))*pow(1.0-S,*(ppar+isumppar+1));
	}
	if(*(wttyp+istat) == 1){
	    nppar = 3;
	    Q = pow(S,*(ppar+isumppar))*pow(1.0-S,*(ppar+isumppar+1));
	    if(xi_ > *(ppar + isumppar+2)) {
	      Q= *(Qstop + istat);
              istop = 1;
	    }
	}
	if(*(wttyp+istat) == 2){
	    nppar = 1;
	    Q =MIN(xi_/(*(ppar+isumppar)),1.0);
	}

	dMU = Q/4.0 * Sc * S_LR * f * w;
        if(qis1orQ==1){
	  dMU *= Q;
	}
	hatX(xi_,th1,lrrf,lrrf_,nnh1,l);
	Beta = lrrf_;
	ans_MU += Beta * dMU;
      }
      *(mufu + nnlook*istat + ilook) = pow(ntrial,0.5) * *(atten + istat) * ans_MU/pow(*(V_END+istat),0.5);
    }
    isumppar += nppar;
  }

  Free(xi);
  Free(H0);
  Free(th1);
  Free(h1);
  Free(H1);
  Free(Hc0);
  Free(Hc1);
  Free(lbuff);
  Free(V_END);
  Free(Qstop);
  Free(atten);
  Free(Qmom_args);
}

void Qmoments(double *pK, double *ph, double *ptc, double *ptr, double *ptau, double *ans)
{
  int p;
  double K, h, tc, tr, tau;
  double H_tc, H_tau_m_tr, H_tau, two_p, xp, I_p, II_p, III_p;

  K = *pK;
  h = *ph;
  tc = *ptc;
  tr = *ptr;
  tau = *ptau;

  H_tc = h * tc;
  H_tau_m_tr = h * (tau - tr);
  H_tau = h * tau;

  for(p=0;p<2;p++)
  {
    two_p = (p==1 ? 2.0 : 1.0);
    xp = (double)p;

    I_p = (1.0 - exp(-(K+1.0)*H_tc))/(K+1.0) - two_p*(1.0-exp(-(K+2.0)*H_tc))/(K+2.0) + xp * (1.0-exp(-(K+3.0)*H_tc))/(K+3.0);
 
    II_p = pow(1-exp(-H_tc), xp+1.0) * (exp(-(K+1.0)*H_tc)-exp(-(K+1.0)*H_tau_m_tr))/(K+1.0);

    III_p = pow(1.0-exp(-H_tc), xp+1.0)/(H_tau - H_tau_m_tr) *(exp(-(K+1.0)*H_tau_m_tr) * (H_tau - H_tau_m_tr)/(K+1.0) -
	   (exp(-(K+1.0) * H_tau_m_tr) - exp(-(K+1.0) * H_tau))/((K+1.0)*(K+1.0)));

    *(ans + p) = (I_p + II_p + III_p)/4.0;
  }
  Rprintf("K:%g, h:%g, tc:%g, tr:%g, tau:%g, m(tau,1)=%g, m(tau,Q)=%g\n",K,h,tc,tr,tau,*ans,*(ans+1));
}

/*
  ans <- double(2)
  for(p in 1:2)
  {
    I.p <- (1-exp(-(K+1)*H.tc))/(K+1) - 2^(p-1)*(1-exp(-(K+2)*H.tc))/(K+2) + (p-1)*(1-exp(-(K+3)*H.tc))/(K+3)

    II.p <- (1-exp(-H.tc))^p * (exp(-(K+1)*H.tc)-exp(-(K+1)*H.tau.m.tr))/(K+1)

    III.p <- (1-exp(-H.tc))^p/(H.tau - H.tau.m.tr) *(exp(-(K+1)*H.tau.m.tr ) * (H.tau - H.tau.m.tr)/(K+1) -
             (exp(-(K+1) * H.tau.m.tr) - exp(-(K+1) * H.tau))/(K+1)^2)

    ans[p] <- (I.p + II.p + III.p)/4
  }
  names(ans) <- c("M_n(\tau, 1)", "V_n(\tau)")
  ans
}
*/
