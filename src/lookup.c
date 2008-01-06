//MACROS
#define fatx(x,xgrid,f,f0,f_,n,l) l=n;while((*(xgrid+l-1)>x) && l>0) l--; f_=f0; if(l!=0) f_= *(f+l-1)
#define MIN(x,y) (x<=y ? x : y)
#define MAX(x,y) (x>=y ? x : y)

void lookup(double *xgrid, double *ygrid, int *pngrid, double *x, int *pnx, 
	    double *py0, double *yatx, int *index)
{
  int ngrid, nx, i, l;
  double y0,y_;
  ngrid = *pngrid;
  nx = *pnx;
  y0 = *py0;
  
  for(i=0;i<nx;i++){
    fatx(*(x+i), xgrid, ygrid, y0, y_, ngrid, l);
    *(yatx+i) = y_;
    *(index+i) = l;
  }
}
