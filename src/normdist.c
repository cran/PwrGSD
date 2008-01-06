#include<R.h>
#include<Rmath.h>

double ratiodnorm(double x,double y)
{
  return exp(-0.5*(x*x - y*y));
}
