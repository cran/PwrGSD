#include<math.h>
#include<R.h>

void Finv(double *tgrid, double *hgrid, long *ntgrid, double *ugrid, long *nugrid, double *xgrid)
{
  long nt, ntm1, nu, i, j;
  double ressv, res;
  nt = *ntgrid;
  ntm1 = nt-1;
  nu = *nugrid;
  for(i=0;i<nu;i++){
    res = -log(1 - *(ugrid+i));
    j=0;
    while(res>0&&j<ntm1){
      ressv = res;
      res -= *(hgrid+j)*(*(tgrid+j+1)- *(tgrid+j));
      j++;
    }
    if(res<0) 
      j--;
    else
      ressv = res;
    *(xgrid + i) = ressv/(*(hgrid + j)) + *(tgrid + j);
  }
}
