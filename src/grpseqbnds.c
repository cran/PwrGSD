 #include<R.h>
#include<Rmath.h>

#define MIN(x,y) (x <= y ? x : y);
#define ratiodnorm(x, y) exp(-0.5*(x*x - y*y))

void  grpseqbndsL(int *pef, double (*spfu)(double frac, double alphatot, double rho), 
                  double *rho, int *islast, int *pnlook, double *palphtot, double *palpha, 
                  double *pfmin, int *dlact, double *pfracold, double *pfracnew, double *pfracold_ii, 
		  double *pfracnew_ii, double *x, double *y, double *tmp, double *intgrndx, 
		  double *gqxw, int *pngqnodes, double *mufu, double *bold, 
		  double *bnew, int *mybound);
void      updateL(int *nbf, int *dofu, int *pef, int *pnlook, double *pfracold, double *pfracnew, 
		  double *x, double *y, double *tmp, double *intgrndx, double *gqxw, 
		  int *pngqnodes, double *mufu, double *bnew);
void  grpseqbndsH(int *islast, int *pnlook, double *palphtot, double *palpha,
		  double *pfracold, double *pfracnew, double *x, double *y, double *tmp,
		  double *intgrndx, double *gqxw, int *pngqnodes, double *mufu, double *bin, 
		  double *blast);
void      updateH(int *dofu, int *islast, int *pnlook, double *pfracold, double *pfracnew, 
		  double *x, double *y, double *tmp, double *intgrndx, double *gqxw, 
		  int *pngqnodes, double *mufu, double *bin, double *blast);
void      printmat(double *pA, int nr, int nc, char *name);
double     obrien(double frac, double alphtot, double rho);
double     pocock(double frac, double alphtot, double rho);
double    powersp(double frac, double alphtot, double rho);
double    (*spfu)(double frac, double alphtot, double rho);

/* 
double pointers, pfracold, pfracnew, pfracold_ii, pfracnew_ii, point to the current
and prior values pointed to by pinffrac and pinffrac_ii.  This allows for the two information
scales: one to spend alpha on and the second to use internally to rescale the Z scores to
brownian motion scale achieving a process with independent increments.  NOTE:  pinffrac_ii is
the alpha spending scale and pinffrac is the original variance scale.
*/

void grpseqbnds(int *dofu, int *nbf, int *nbnd, int *nsf, double *rho, int *pnlook, double *palphtot, 
		double *palpha, double *psimin, int *dlact, double *pfracold, double *pfracnew, 
		double *pfracold_ii, double *pfracnew_ii, double *x, double *y, double *tmp, 
		double *intgrndx, double *gqxw, int *pngqnodes, double *mufu, double *bold, 
		double *bnew, int *mybound)
{
  double denom;
  double *pfmin;
  int ef, ngq;
  int *pef, *islast;

  islast = (int *)Calloc(1, int);
  pef = (int *)Calloc(1, int);
  pfmin = (double *)Calloc(1, double);

  *islast = (1.0 - *pfracnew_ii < 1.0e-6);

  ngq = *pngqnodes;

  for(ef=0;ef<1+*dofu;ef++){
    *pef = ef;
    if(*(nsf+ef)==1) {
      spfu = &obrien;
      denom = qnorm5(1.0-*psimin,0.0,1.0,1,0);
      *pfmin = qnorm5(1.0-*(palphtot+ef)/2.0,0.0,1.0,1,0)/denom;
      *pfmin = *pfmin* *pfmin;
    }
    if(*(nsf+ef)==2) {
      spfu = &pocock;
      *pfmin = (exp(*psimin/(*(palphtot+ef)))-1.0)/(exp(1.0)-1.0);
    }
    if(*(nsf+ef)==3) {
      spfu = &powersp;
      *pfmin = pow(*psimin/(*(palphtot+ef)),1.0/(*(rho+ef)));
    }
    if(*(nbnd+ef)==1 || *(nbnd+ef)==3 || *(nbnd+ef)==4)
      grpseqbndsL(pef, spfu, rho+ef, islast, pnlook+ef, palphtot+ef, palpha+ef, pfmin, 
		  dlact+ef, pfracold + ef, pfracnew, pfracold_ii + ef, pfracnew_ii, x+ef*ngq, y+ef*ngq, 
		  tmp+ef*ngq, intgrndx+ef*ngq, gqxw, pngqnodes, mufu, bold, bnew, mybound);
    if(*(nbnd+ef)==2){
      *(dlact+ef) = 1;
      grpseqbndsH(islast, pnlook+ef, palphtot+ef, palpha+ef, pfracold + ef, pfracnew, 
		  x+ef*ngq, y+ef*ngq, tmp+ef*ngq, intgrndx+ef*ngq, 
		  gqxw, pngqnodes, mufu, bold+ef, bnew+ef);
    }
  }
  if(*islast==0){
    for(ef=0;ef<1+*dofu;ef++){
      *pef = ef;
      if((*(nbnd+ef)==1 || *(nbnd+ef)==3 || *(nbnd+ef)==4) && *(dlact+ef)==1)
	updateL(nbf, dofu, pef, pnlook+ef, pfracold + ef, pfracnew, x+ef*ngq, y+ef*ngq, 
		tmp+ef*ngq, intgrndx+ef*ngq, gqxw, pngqnodes, mufu, bnew);
      if(*(nbnd+ef)==2)
	updateH(dofu, islast, pnlook+ef, pfracold + ef, pfracnew, x+ef*ngq, y+ef*ngq, 
		tmp+ef*ngq, intgrndx+ef*ngq, gqxw, pngqnodes, mufu, bold+ef, bnew);
    }
  }
//printf("nbnd:%d, islast:%d, nlook:%d, alphatot:%g, alpha:%g, fold:%g, fnew:%g, bold:%g, bnew:%g\n",
//	 *nbnd, *islast, *pnlook, *palphtot, *palpha, *pfracold, *pfracnew, *bold, *bnew);

  Free(islast);
  Free(pef);
  Free(pfmin);
}

void grpseqbndsL(int *pef, double (*spfu)(double frac, double alphatot, double rho), 
                 double *rho, int *islast, int *pnlook, double *palphtot, double *palpha, 
                 double *pfmin, int *dlact, double *pfracold, double *pfracnew, double *pfracold_ii, 
		 double *pfracnew_ii, double *x, double *y, double *tmp, double *intgrndx, 
		 double *gqxw, int *pngqnodes, double *mufu, double *bold, double *bnew,
		 int *mybound)
{
  double vsmall=1.0e-6, vvsmall=1.0e-15, ltone=7.0, utzero=18.66, sw, x_, dx_;
  double psimin, aold, anew, sqrf, sqrdf, b, Phib,bl,bu,berr,aerr,aerrsgn,intgrl,yy;
  double *gqx, *gqw;
  int nlook,nlkm1,ifault,i,j,ef,ngqnodes,hangs,zero=0, one=1;

  ngqnodes=*pngqnodes;
  gqx = gqxw;
  gqw = gqxw + ngqnodes;

  nlook = *pnlook;
  nlkm1 = nlook - 1;
  ef = *pef;
  sw = (double) ef;
  psimin = (*spfu)(*pfmin, *palphtot, *rho);

  *dlact = 0;
  aold = 0.0;
  anew = psimin;
  if(*pfracold_ii > *pfmin) aold = (*spfu)(*pfracold_ii, *palphtot, *rho);
  if(*pfracnew_ii > *pfmin || *(mybound + ef)==1) {
    anew = (*spfu)(*pfracnew_ii, *palphtot, *rho);
    *dlact = 1;
  }

  *palpha = anew - aold;

  sqrf = pow(*pfracnew,0.5);
  sqrdf = pow(*pfracnew - *pfracold,0.5);

  if(*dlact==1 && (*islast==0 || ef==0) && *(mybound + ef)==0){
    if(nlook==1)
      b = qnorm5(*palpha,0.0,1.0,ef,zero) + sw * *mufu/sqrf;
    else{
      bl=(1.0-sw) * vsmall + sw * *bold;
      bu=(1.0-sw) *  *bold + sw * *(bold+1);
      b = (bl+bu)/2.0;
      berr = (bu-bl)/2.0;
      aerr = 1.0;
      hangs=0;
      while((berr>vsmall||aerr>vvsmall)&&hangs<300){
	Phib = pnorm5(sqrf * b - sw * *mufu,0.0,1.0,one,zero);
	intgrl = 0.0;
	for (j=0;j<ngqnodes;j++){
	  x_  = (1.0-sw) * ((1.0-*(gqx+j))/2.0 * Phib + (1.0+*(gqx+j))/2.0) +
	             sw  *  (1.0+*(gqx+j))/2.0 * Phib;
	  dx_ = (1.0-sw) * (1.0-Phib)/2.0*(*(gqw+j)) + sw *  Phib/2.0*(*(gqw+j));
	  yy = qnorm5(x_,0.0,1.0,one,zero);
	  for(i=0;i<ngqnodes;i++)
	    intgrl += ratiodnorm((yy-*(x+i))/sqrdf, yy) * dx_/sqrdf * (*(intgrndx+i));
	}
	aerr = *palpha - intgrl;
	aerrsgn = (aerr >= 0.0 ? 1.0 : -1.0);
	aerr = ((double) aerrsgn) * aerr;
	if(aerrsgn<0) bl = b;
	else bu = b;
	b = (bl+bu)/2.0;
	berr=fabs(bu-bl)/2.0;
	hangs++;
      }
    }
    *(bnew+ef)=b;
  }
  if(*(mybound + ef)==1){
    if(*islast==1 && ef==1) *(bnew+1) = *bnew;
    b = *(bnew+ef);
    Phib = pnorm5(sqrf * b - sw * *mufu,0.0,1.0,one,zero);
    intgrl = 0.0;
    for (j=0;j<ngqnodes;j++){
      x_  = (1.0-sw) * ((1.0-*(gqx+j))/2.0 * Phib + (1.0+*(gqx+j))/2.0) +
	         sw  *  (1.0+*(gqx+j))/2.0 * Phib;
      dx_ = (1.0-sw) * (1.0-Phib)/2.0*(*(gqw+j)) + sw *  Phib/2.0*(*(gqw+j));
      yy = qnorm5(x_,0.0,1.0,one,zero);
      for(i=0;i<ngqnodes;i++)
	intgrl += ratiodnorm((yy-*(x+i))/sqrdf, yy) * dx_/sqrdf * (*(intgrndx+i));
    }
    *palpha = intgrl;
  }
  if(*(mybound + 1)==0 && *islast==1 && ef==1) {
    *(bnew+1) = *bnew;
    b = *bnew;
    Phib = pnorm5(sqrf * b - *mufu,0.0,1.0,one,zero);
    intgrl = 0.0;
    for (j=0;j<ngqnodes;j++){
      x_  = (1.0+*(gqx+j))/2.0 * Phib;
      dx_ = Phib/2.0*(*(gqw+j));
      yy = qnorm5(x_,0.0,1.0,one,zero);
      for(i=0;i<ngqnodes;i++)
	intgrl += ratiodnorm((yy-*(x+i))/sqrdf, yy) * dx_/sqrdf * (*(intgrndx+i));
    }
    *palpha = intgrl;
  }
}

void grpseqbndsH(int *islast, int *pnlook, double *palphtot, double *palpha,
		 double *pfracold, double *pfracnew, double *x, double *y, double *tmp,
		 double *intgrndx, double *gqxw, int *pngqnodes, double *mufu,
		 double *bin, double *blast)
{
  double vsmall=1.0e-6, vvsmall=1.0e-15, ltone=7.0, utzero=18.66;
  double sqrf, sqrdf, b, Phib,bl,bu,berr,aerr,aerrsgn,intgrl,yy;
  double *gqx, *gqw;
  int nlook,nlkm1,ifault,i,j,ngqnodes,hangs,zero=0, one=1;

  ngqnodes=*pngqnodes;
  gqx = gqxw;
  gqw = gqxw + ngqnodes;

  nlook = *pnlook;
  nlkm1 = nlook - 1;
  sqrf = pow(*pfracnew,0.5);
  sqrdf = pow(*pfracnew - *pfracold,0.5);
  if(nlook==1) {
    if(*islast==0) {
      *palpha = 1.0 - pnorm5(*bin, 0.0, 1.0, one, zero);
      b = *bin;
    }
    if(*islast==1) {
      *palpha = *palphtot;
      *blast = qnorm5(1.0-*palpha,0.0,1.0,one,zero);
      b = *bin;
    }
  }
  
  if(nlook>1&&(*islast==0)){
    b = *bin;
    Phib = pnorm5(sqrf * b,0.0,1.0,one,zero);
    intgrl = 0.0;
    for (j=0;j<ngqnodes;j++){
      yy = qnorm5((1.0-*(gqx+j))/2.0 * Phib + (1.0+*(gqx+j))/2.0,0.0,1.0,one,zero);
      for(i=0;i<ngqnodes;i++)
	intgrl += ratiodnorm((yy-*(x+i))/sqrdf, yy) *(1-Phib)/2.0*(*(gqw+j))/sqrdf * (*(intgrndx+i));
    }
    *palpha = intgrl;
  }
  
  if(nlook>1&&(*islast==1)){
    bl=vsmall;
    bu=*bin;
    b = (bl+bu)/2.0;
    berr = (bu-bl)/2.0;
    aerr = 1.0;
    hangs=0;
    while((berr>vsmall||aerr>vvsmall)&&hangs<300){
      Phib = pnorm5(sqrf * b,0.0,1.0,one,zero);
      intgrl = 0.0;
      for (j=0;j<ngqnodes;j++){
	yy = qnorm5((1.0-*(gqx+j))/2.0 * Phib + (1.0+*(gqx+j))/2.0,0.0,1.0,one,zero);
	for(i=0;i<ngqnodes;i++)
	  intgrl += ratiodnorm((yy-*(x+i))/sqrdf, yy) *(1-Phib)/2.0*(*(gqw+j))/sqrdf * (*(intgrndx+i));
      }
      aerr = *palphtot - intgrl;
      aerrsgn = (aerr >= 0.0 ? 1.0 : -1.0);
      aerr = aerrsgn * aerr;
      if(aerrsgn<0) bl = b;
      else bu = b;
      b = (bl+bu)/2.0;
      berr=(bu-bl)/2.0;
      hangs++;
    }
    *blast = b;
    *palpha = intgrl;
  }
}

// Alpha Spending functions for stopping boundary constuction:
//1.  Obrien-Flemming Spending: 
double obrien(double frac, double alphtot, double rho)
{
  int one=1,zero=0;
  return(2.0*pnorm5(qnorm5(1.0-alphtot/2.0,0.0,1.0,one,zero)/pow(frac,0.5),0.0,1.0,zero,zero));
}

//2.  Pocock Spending: 
double pocock(double frac, double alphtot, double rho)
{
  double tmp,e,ans;
  e = exp(1);
  tmp = alphtot * log(1.0+(e-1.0)*frac);
  ans = (tmp <= alphtot ? tmp : alphtot);
  return(ans);
}

//3.  Power Spending:
double powersp(double frac, double alphtot, double rho)
{
  double ans;
  ans = MIN(alphtot, alphtot*pow(frac, rho));
  return(ans);
}

// add definintion for alternate spending function here...
// don't forget to make corresponding changes in this program and the
// R function interface to expand the menu of choices, and don't forget
// to add include the corresponding function prototype declaration above.


void updateL(int *nbf, int *dofu, int *pef, int *pnlook, double *pfracold, double *pfracnew, 
	     double *x, double *y, double *tmp, double *intgrndx, double *gqxw,
	     int *pngqnodes, double *mufu, double *bnew)
{
  int ngq, nlook, i, j, one=1, zero=0, nlkm1, ef;
  double a, b, Phia, Phib, sqrf, sqrdf, sw; 
  double *gqx, *gqw;

  ngq = *pngqnodes;
  gqx = gqxw;
  gqw = gqxw + ngq;

  nlook = *pnlook;
  nlkm1 = nlook -1;
  ef = *pef;
  sw = (double) ef;
  sqrf = pow(*pfracnew,0.5);
  sqrdf = pow(*pfracnew - *pfracold,0.5);
  a=-9.0;
  if(*dofu==1 && (*nbf==0 || ef==1)){
    a = *(bnew+1);
    Phia = pnorm5(sqrf * a - sw * *mufu,0.0,1.0,one,zero);
  }
  else Phia = 0.0;
  b = *bnew;
  Phib = pnorm5(sqrf * b - sw * *mufu,0.0,1.0,one,zero);
  if(nlook==1){
    for(i=0;i<ngq;i++){
      *(y+i) = qnorm5((1.0-*(gqx+i))/2.0 * Phia + (1.0+*(gqx+i))/2.0 * Phib,0.0,1.0,one,zero);
      *(tmp+i)= ratiodnorm(*(y+i)/sqrdf, *(y+i))*(Phib-Phia)/2.0*(*(gqw+i))/sqrdf;
    }
  }
  else
    for(j=0;j<ngq;j++) {
      *(tmp+j) = 0.0;
      *(y+j) = qnorm5((1.0-*(gqx+j))/2.0*Phia + (1.0+*(gqx+j))/2.0*Phib,0.0,1.0,one,zero);
      for(i=0;i<ngq;i++)
	*(tmp+j) += ratiodnorm((*(y+j) - *(x+i))/sqrdf, *(y+j)) * (Phib-Phia)/2.0 * (*(gqw+j))/sqrdf * 
	            (*(intgrndx+i));
    }

  for(i=0;i<ngq;i++) {
    *(intgrndx+i) = *(tmp+i);
    *(x+i) = *(y+i);
  }
}


void updateH(int *dofu, int *islast, int *pnlook, double *pfracold, double *pfracnew, 
	     double *x, double *y, double *tmp, double *intgrndx, double *gqxw,
	     int *pngqnodes, double *mufu, double *bin, double *blast)
{
  int nlook, ngq, i, j, one=1, zero=0, ef;
  double a, b, Phia, Phib, sqrf, sqrdf;
  double *gqx, *gqw;

  nlook = *pnlook;
  ngq = *pngqnodes;
  gqx = gqxw;
  gqw = gqxw + ngq;

  sqrf = pow(*pfracnew,0.5);
  sqrdf = pow(*pfracnew - *pfracold,0.5);
  if(*dofu==1){
    a = *(blast+1);
    Phia = pnorm5(sqrf * a,0.0,1.0,one,zero);
  }
  else Phia = 0.0;
  b = *bin;
  if(*islast==1&nlook>1) b = *blast;
  Phib=pnorm5(sqrf * b,0.0,1.0,one,zero);

  if(nlook==1)
    for(i=0;i<ngq;i++){
      *(y+i) = qnorm5((1.0-*(gqx+i))/2.0*Phia + (1.0+*(gqx+i))/2.0*Phib,0.0,1.0,one,zero);
      *(tmp+i)= ratiodnorm(*(y+i)/sqrdf, *(y+i))*(Phib-Phia)/2.0*(*(gqw+i))/sqrdf;
    }
  else
    for(j=0;j<ngq;j++) {
      *(tmp+j) = 0.0;
      *(y+j) = qnorm5((1.0-*(gqx+j))/2.0*Phia + (1.0+*(gqx+j))/2.0*Phib,0.0,1.0,one,zero);
      for(i=0;i<ngq;i++)
	*(tmp+j) += ratiodnorm((*(y+j)-*(x+i))/sqrdf,*(y+j))*(Phib-Phia)/2.0*(*(gqw+j))/sqrdf 
	            * (*(intgrndx+i));
    }

  for(i=0;i<ngq;i++) {
    *(intgrndx+i) = *(tmp+i);
    *(x+i) = *(y+i);
  }
}
