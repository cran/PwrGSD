/* the Dnsty of the statistic (k, x) for k=1, 2,..,K and x in the rejection region should be 
   easily obtained from the pointer 'intgrndx' below here inside grpseqbnds.c
   step 1: pair grpseqbnds down to the bare minimum-- i.e. I don't need the Haybittle
           stuff-- I don't need futility stuff, I don't even need a drift.  Pull these
           all out and verify that grpseqbnds (this program is called Dnsty) still
           gives boundary values which agree with grpseqbnds
   step 2: lose the 'compute the boundary' part--boundary will be provided as an argument,
   step 3: collect the vectors intgrndx at each analysis k=1, 2,...,K into a matrix and
           return.  Make sure that I didn't lose a constant somewhere. 
*/
#include<R.h>
#include<Rmath.h>

#define MIN(x,y) (x <= y ? x : y);
#define ratiodnorm(x, y) exp(-0.5*(x*x - y*y))

void   dnstyKern(double (*spfu)(double frac, double alphatot, double rho), 
                 double *rho, int *islast, int *pnlook, double *palphtot, double *palpha, 
                 double *pfmin, int *dlact, double *pfracold, double *pfracnew, double *pfracold_ii, 
	         double *pfracnew_ii, double *xc, double *y, double *tmp, double *dF, double *gqxw, 
                 int *pngqnodes, double *pxc, double *pxr, double *bold, double *bnew, int *mybound);
void updateIntgr(int *nbf, int *pnlook, double *pfracold, double *pfracnew, double *xc, 
                  double *y, double *tmp, double *dF, double *gqxw, int *pngqnodes, 
                  double *bnew);
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

void dnsty(int *nbf, int *nbnd, int *nsf, double *rho, int *pnlook, double *palphtot, 
	   double *palpha, double *psimin, int *dlact, double *pfracold, double *pfracnew, 
	   double *pfracold_ii, double *pfracnew_ii, double *xc, double *y, double *tmp, 
           double *dF, double *gqxw, int *pngqnodes, double *pxc, double *pxr, double *bold, 
           double *bnew, int *mybound)
{
  double denom;
  double *pfmin;
  int i,ngq;
  int *islast;

  islast = (int *)Calloc(1, int);
  pfmin = (double *)Calloc(1, double);

  *islast = (1.0 - *pfracnew_ii < 1.0e-6);

  ngq = *pngqnodes;
  /*  ef=0 dofu=0 */
  if(*nsf==1) {
    spfu = &obrien;
    denom = qnorm5(1.0-*psimin,0.0,1.0,1,0);
    *pfmin = qnorm5(1.0-*palphtot/2.0,0.0,1.0,1,0)/denom;
    *pfmin = *pfmin* *pfmin;
  }
  if(*nsf==2) {
    spfu = &pocock;
    *pfmin = (exp(*psimin/(*palphtot))-1.0)/(exp(1.0)-1.0);
  }
  if(*nsf==3) {
    spfu = &powersp;
    *pfmin = pow(*psimin/(*palphtot),1.0/(*rho));
  }
  if(*nbnd==1 || *nbnd==3 || *nbnd==4)
    dnstyKern(spfu, rho, islast, pnlook, palphtot, palpha, pfmin, dlact, pfracold, pfracnew, 
              pfracold_ii, pfracnew_ii, xc, y, tmp, dF, gqxw, pngqnodes, pxc, pxr, bold, bnew, mybound);

  /*take *nbnd==2 out as a possibility (Haybittle) */

  if(*islast==0){
    /*  ef=0 dofu=0 */
    if((*nbnd==1 || *nbnd==3 || *nbnd==4) && *dlact==1)
      updateIntgr(nbf, pnlook, pfracold, pfracnew, xc, y, tmp, dF, gqxw, pngqnodes, bnew);
  }

//printf("nbnd:%d, islast:%d, nlook:%d, alphatot:%g, alpha:%g, fold:%g, fnew:%g, bold:%g, bnew:%g\n",
//	 *nbnd, *islast, *pnlook, *palphtot, *palpha, *pfracold, *pfracnew, *bold, *bnew);

  Free(islast);
  Free(pfmin);
}

void dnstyKern(double (*spfu)(double frac, double alphatot, double rho), 
               double *rho, int *islast, int *pnlook, double *palphtot, double *palpha, 
               double *pfmin, int *dlact, double *pfracold, double *pfracnew, double *pfracold_ii, 
	       double *pfracnew_ii, double *xc, double *y, double *tmp, double *dF, double *gqxw, 
               int *pngqnodes, double *pxc, double *pxr, double *bold, double *bnew, int *mybound)
{
  double vsmall=1.0e-6, vvsmall=1.0e-15, ltone=7.0, utzero=18.66, sw, xc_, dxc_, xr_, dxr_;
  double psimin, aold, anew, sqrf, sqrdf, b, Phib,bl,bu,berr,aerr,aerrsgn,intgrl,yyc_,yyr_;
  double *gqx, *gqw;
  int nlook,nlkm1,ifault,i,j,ngqnodes,hangs,zero=0, one=1;

  ngqnodes=*pngqnodes;
  gqx = gqxw;
  gqw = gqxw + ngqnodes;

  nlook = *pnlook;
  nlkm1 = nlook - 1;

  /* dofu=0, *pef=0, ef=0, sw=0.0 */
  psimin = (*spfu)(*pfmin, *palphtot, *rho);

  *dlact = 0;
  aold = 0.0;
  anew = psimin;
  if(*pfracold_ii > *pfmin) aold = (*spfu)(*pfracold_ii, *palphtot, *rho);
  if(*pfracnew_ii > *pfmin || *mybound==1) {
    anew = (*spfu)(*pfracnew_ii, *palphtot, *rho);
    *dlact = 1;
  }

  *palpha = anew - aold;

  sqrf = pow(*pfracnew,0.5);
  sqrdf = pow(*pfracnew - *pfracold,0.5);

  if(*dlact==1 && *mybound==0){
    if(nlook==1)
      b = qnorm5(*palpha,0.0,1.0,zero,zero);
    else{
      bl=vsmall;
      bu=*bold;
      b = (bl+bu)/2.0;
      berr = (bu-bl)/2.0;
      aerr = 1.0;
      hangs=0;
      while((berr>vsmall||aerr>vvsmall)&&hangs<300){
	Phib = pnorm5(sqrf * b,0.0,1.0,one,zero);
	intgrl = 0.0;
	for (j=0;j<ngqnodes;j++){
	  xr_  = ((1.0-*(gqx+j))/2.0 * Phib + (1.0+*(gqx+j))/2.0);
	  dxr_ = (1.0-Phib)/2.0*(*(gqw+j));
	  yyr_ = qnorm5(xr_,0.0,1.0,one,zero);
         
	  for(i=0;i<ngqnodes;i++)
	    intgrl += ratiodnorm((yyr_ - *(xc+i))/sqrdf, yyr_) * dxr_/sqrdf * (*(dF + i));
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
    *bnew=b;

    for(i=0;i<ngqnodes;i++)
    {
      Phib = pnorm5(sqrf * b,0.0,1.0,one,zero);
      xr_  = (1.0-*(gqx+i))/2.0 * Phib + (1.0+*(gqx+i))/2.0;
      dxr_ = (1.0-Phib)/2.0*(*(gqw+i));
      yyr_ = qnorm5(xr_,0.0,1.0,one,zero);
      *(pxr + i) = yyr_;

      xc_  = (1.0+*(gqx+i)) * Phib/2.0;
      dxc_ = Phib/2.0*(*(gqw+i));
      yyc_ = qnorm5(xc_,0.0,1.0,one,zero);
      *(pxc + i) = yyc_;
    }
  }
  if(*mybound==1){
    b = *bnew;
    Phib = pnorm5(sqrf * b,0.0,1.0,one,zero);
    intgrl = 0.0;
    for (j=0;j<ngqnodes;j++){
      xr_  = (1.0-*(gqx+j))/2.0 * Phib + (1.0+*(gqx+j))/2.0;
      dxr_ = (1.0-Phib)/2.0*(*(gqw+j));
      yyr_ = qnorm5(xr_,0.0,1.0,one,zero);

      for(i=0;i<ngqnodes;i++)
	intgrl += ratiodnorm((yyr_ - *(xc + i))/sqrdf, yyr_) * dxr_/sqrdf * (*(dF + i));
    }
    *palpha = intgrl;
  }
}

void updateIntgr(int *nbf, int *pnlook, double *pfracold, double *pfracnew, double *xc, double *y, 
                 double *tmp, double *dF, double *gqxw, int *pngqnodes, double *bnew)
{
  int ngq, nlook, i, j, one=1, zero=0, nlkm1;
  double b, Phib, sqrf, sqrdf; 
  double *gqx, *gqw;

  ngq = *pngqnodes;
  gqx = gqxw;
  gqw = gqxw + ngq;

  nlook = *pnlook;
  nlkm1 = nlook -1;
  /* ef=0, sw=0.0, *dofu=0 */
  sqrf = pow(*pfracnew,0.5);
  sqrdf = pow(*pfracnew - *pfracold,0.5);

  b = *bnew;
  Phib = pnorm5(sqrf * b,0.0,1.0,one,zero);
  if(nlook==1){
    for(i=0;i<ngq;i++){
      *(y+i) = qnorm5((1.0+*(gqx+i))/2.0 * Phib,0.0,1.0,one,zero);
      *(tmp+i)= ratiodnorm(*(y+i)/sqrdf, *(y+i))*Phib/2.0*(*(gqw+i))/sqrdf;
    }
  }
  else
    for(j=0;j<ngq;j++) {
      *(tmp+j) = 0.0;
      *(y+j) = qnorm5((1.0+*(gqx+j))/2.0*Phib,0.0,1.0,one,zero);
      for(i=0;i<ngq;i++)
	*(tmp+j) += ratiodnorm((*(y+j) - *(xc + i))/sqrdf, *(y+j)) * Phib/2.0 * (*(gqw+j))/sqrdf * (*(dF + i));
    }

  for(i=0;i<ngq;i++) {
    *(dF + i) = *(tmp+i);
    *(xc + i) = *(y+i);
  }
}
