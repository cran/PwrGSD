#include<R.h>
//
//MACROS 
#define COMPH(xh,h,H,n,l) *H=0.0;for(l=1;l<n;l++) *(H+l)=*(H+l-1)+ *(h+l-1) * (*(xh+l) - *(xh+l-1))
#define hHatX(x,xh,h,H,h_,H_,n,l) l=0;while(*(xh+l)<=x && l<n) l++; h_ = *(h+l-1); H_ = *(H+l-1) + h_ * (x-*(xh+l-1))
#define HIatW(W,xh,h,H,HI_,n,l) l=0;while(*(H+l)<=W && l<n) l++; HI_ = *(xh+l-1) + (W - *(H+l-1))/(*(h+l-1))
//
//
void htildeConst(double *x,int *nx,double *h,double *hA,double *hB,
                 double *lA,double *lB,double *Stilde,double *htlde);

void htilde(double *x,int *nx,double *gqx,double *gqw,int *ngq,
            double *xh,double *h,int *nh,double *xhA,double *hA,int *nhA,
            double *xhB,double *hB,int *nhB,double *xlA,double *lA,int *nlA,
            double *xlB,double *lB,int *nlB,int *gradual,double *tend,
	    double *ftlde,double *Stlde,double *htlde)
{
  double xx,xiA,xiB,h_xiA,H_xiA,h_xiB,H_xiB,hA_xiA,HA_xiA,hB_xiB,HB_xiB;
  double h_x,H_x,hA_x,HA_x,hB_x,HB_x,DS,SlAx,SlBx,LAtend,LBtend;
  double lA_,LA_,lB_,LB_,lA_xiB,LA_xiB,lB_xiA,LB_xiA,SlA_xiB,SlB_xiA;
  double LA_B,tauA_B,h_AB,H_AB,hA_AB,HA_AB,LB_A,tauB_A,h_BA,H_BA,hB_BA,HB_BA;
  double hA_x_,hB_x_,D_HA,D_HB,piA,piB,thA,thB,wA,wB,phiA_B,wA_B,phiB_A,wB_A,h_,H_;
  double SI_1, SI_2, SI_3, SI_4, SI_5, fI_1, fI_2, fI_3, fI_4, fI_5, tmp;
  double *H, *HA, *HB, *LA, *LB;
  int nnx, nngq, nnh, nnhA, nnhB, nnlA, nnlB, i, l, j, k, idx;

  nnx = *nx;
  nngq = *ngq;
  nnh = *nh;
  nnhA = *nhA;
  nnhB = *nhB;
  nnlA = *nlA;
  nnlB = *nlB;

  H = (double *)Calloc(nnh, double);
  HA = (double *)Calloc(nnhA, double);
  HB = (double *)Calloc(nnhB, double);
  LA = (double *)Calloc(nnlA, double);
  LB = (double *)Calloc(nnlB, double);

  COMPH(xh,h,H,nnh,l);
  COMPH(xhA,hA,HA,nnhA,l);
  COMPH(xhB,hB,HB,nnhB,l);
  COMPH(xlA,lA,LA,nnlA,l);
  COMPH(xlB,lB,LB,nnlB,l);
  
  hHatX(*tend,xlA,lA,LA,lA_,LAtend,nnlA,l);
  hHatX(*tend,xlB,lB,LB,lB_,LBtend,nnlB,l);
  
  for(i=0;i<nnx;i++){
    xx = *(x+i);
    hHatX(xx, xh, h, H, h_, H_, nnh,l);
    hHatX(xx,xlA,lA,LA,lA_,LA_,nnlA,l);
    hHatX(xx,xlB,lB,LB,lB_,LB_,nnlB,l);
    hHatX(xx,xh,h,H,h_x,H_x,nnhB,l);
    hHatX(xx,xhA,hA,HA,hA_x,HA_x,nnhA,l);
    hHatX(xx,xhB,hB,HB,hB_x,HB_x,nnhB,l);

    SlAx = exp(-LA_);
    SlBx = exp(-LB_);

    SI_1 = exp(-H_) * SlAx * SlBx;
    fI_1 = h_ * SI_1;

    SI_2 = 0.0;
    fI_2 = 0.0;
    SI_3 = 0.0;
    fI_3 = 0.0;
    SI_4 = 0.0;
    fI_4 = 0.0;
    SI_5 = 0.0;
    fI_5 = 0.0;
    for(j=0;j<nngq;j++){
      thA = (1.0 - SlAx) * (1.0 + *(gqx + j))/2.0;  
      wA = *(gqw + j)*(1.0 - SlAx)/2.0;
      LA_ = -log(1.0-thA);

      thB = (1.0 - SlBx) * (1.0 + *(gqx + j))/2.0;  
      wB = *(gqw + j)*(1.0 - SlBx)/2.0;
      LB_ = -log(1.0-thB);

      HIatW(LA_, xlA, lA, LA, xiA, nnlA, l);
      HIatW(LB_, xlB, lB, LB, xiB, nnlB, l);

      hHatX(xiA, xh, h, H, h_xiA, H_xiA, nnh,l);
      hHatX(xiB, xh, h, H, h_xiB, H_xiB, nnh,l);

      hHatX(xiA,xhA,hA,HA,hA_xiA,HA_xiA,nnhA,l);
      hHatX(xiB,xhB,hB,HB,hB_xiB,HB_xiB,nnhB,l);

      piA = piB = 0.0;
      if(*gradual==1){
	piA = LA_/LAtend;
	piB = LB_/LBtend;
      }
 
      hA_x_ = piA * h_x + (1.0 - piA) * hA_x;
      D_HA = piA * (H_x - H_xiA) + (1.0 - piA) * (HA_x - HA_xiA);
      DS = exp(-(D_HA + H_xiA)) * wA;
      SI_2 += DS;
      fI_2 += hA_x_ * DS;

      hB_x_ = piB * h_x + (1.0 - piB) * hB_x;
      D_HB = piB * (H_x - H_xiB) + (1.0 - piB) * (HB_x - HB_xiB);
      DS = exp(-(D_HB + H_xiB)) * wB;
      SI_3 += DS;
      fI_3 += hB_x_ * DS;
      
      hHatX(xiB, xlA, lA, LA, lA_xiB, LA_xiB, nnlA, l);
      hHatX(xiA, xlB, lB, LB, lB_xiA, LB_xiA, nnlB, l);
      
      SlA_xiB = exp(-LA_xiB);
      SlB_xiA = exp(-LB_xiA);

      for(k=0;k<nngq;k++){
	phiA_B = (1.0 - SlA_xiB) * (1.0 + *(gqx + k))/2.0;
	wA_B = *(gqw + k) * (1.0 - SlA_xiB)/2.0;
	LA_B = -log(1.0 - phiA_B);
	HIatW(LA_B, xlA, lA, LA, tauA_B, nnlA, l);
	hHatX(tauA_B, xh, h, H, h_AB, H_AB, nnh, l);
	hHatX(tauA_B, xhA, hA, HA, hA_AB, HA_AB, nnhA, l);

	phiB_A = (1.0 - SlB_xiA) * (1.0 + *(gqx + k))/2.0;
	wB_A = *(gqw + k) * (1.0 - SlB_xiA)/2.0;
	LB_A = -log(1.0 - phiB_A);
	HIatW(LB_A, xlB, lB, LB, tauB_A, nnlB, l);
	hHatX(tauB_A, xh, h, H, h_BA, H_BA, nnh, l);
	hHatX(tauB_A, xhB, hB, HB, hB_BA, HB_BA, nnhB, l);

	piA=piB=0.0;
	if(*gradual==1){
	  piA = LA_B/LAtend;
	  piB = LB_A/LBtend;
	}

        hA_x_ = piA * h_x + (1.0 - piA) * hA_x;
	D_HA = piA * (H_x - H_AB) + (1.0 - piA) * (HA_x - HA_AB);
	DS = exp(-(D_HA + H_AB)) * wA_B * wB;
	SI_4 += DS;
	fI_4 += hA_x_ * DS;

        hB_x_ = piB * h_x + (1.0 - piB) * hB_x;
	D_HB = piB * (H_x - H_BA) + (1.0 - piB) * (HB_x - HB_BA);
	DS = exp(-(D_HB + H_BA)) * wB_A * wA;
	SI_5 += DS;
	fI_5 += hB_x_ * DS;
      }
    }
    SI_2 = SI_2 * SlBx;
    fI_2 = fI_2 * SlBx;
    SI_3 = SI_3 * SlAx;
    fI_3 = fI_3 * SlAx;

    *(Stlde + i) = (SI_1 + SI_2 + SI_3 + SI_4 + SI_5);
    *(ftlde+i) = (fI_1 + fI_2 + fI_3 + fI_4 + fI_5);
    *(htlde + i) = *(ftlde+i)/(*(Stlde+i));
  }
  Free(H); 
  Free(HA);
  Free(HB);
  Free(LA);
  Free(LB);
}

void htildeConst(double *x,int *nx,double *h,double *hA,double *hB,double *lA,double *lB,
		 double *Stlde,double *htlde)
{
  double xx;
  double SI_1, SI_2, SI_3, SI_4, SI_5, fI_1, fI_2, fI_3, fI_4, fI_5;
  int nnx, i;
  nnx = *nx;

  for(i=0;i<nnx;i++){
    xx = *(x+i);
    SI_1 = exp(-(*h + *lA + *lB)* xx);
    SI_2 = *lA * exp(-(*hA + *lB) * xx)/(*h - *hA + *lA) * 
           (1.0 - exp(-(*h - *hA + *lA) * xx));
    SI_3 = *lB * exp(-(*hB + *lA) * xx)/(*h - *hB + *lB) *
           (1.0 - exp(-(*h - *hB + *lB) * xx));
    SI_4 = *lA * exp(- *hA * xx)/(*h - *hA + *lA) *
           (1.0 - exp(- *lB * xx) - *lB/(*h - *hA + *lA + *lB) *
	    (1.0 - exp(-(*h - *hA + *lA + *lB) * xx)));
    SI_5 = *lB * exp(- *hB * xx)/(*h - *hB + *lB) *
           (1.0 - exp(- *lA * xx) - *lA/(*h - *hB + *lB + *lA) *
            (1.0 - exp(-(*h - *hB + *lB + *lA) * xx)));
    
    fI_1 = *h * SI_1;
    fI_2 = *hA * SI_2;
    fI_3 = *hB * SI_3;
    fI_4 = *hA * SI_4;
    fI_5 = *hB * SI_5;

    *(Stlde + i) = SI_1 + SI_2 + SI_3 + SI_4 + SI_5;
    *(htlde + i) = (fI_1 + fI_2 + fI_3 + fI_4 + fI_5)/(*(Stlde + i));
  }
}
