#include<R.h>

void    randfromh(long *pn, double *tcut, double *h, long *pncut,
                  double *t);
void    randhcdtl(long *pn, double *tcut, double *h, long *pncut,
                  double *tcutxA, double *hxA, long *pncutxA, double *tdA,
                  double *tcutxB, double *hxB, long *pncutxB, double *tdB,
                  long *code, double *t);

void trandfromh(long *pn, double *tcut, double *h, long *pncut,
                  double *t)
{
  GetRNGstate();
  randfromh(pn, tcut,  h, pncut, t);
  PutRNGstate();
}

void trandhcdtl(long *pn, double *tcut, double *h, long *pncut,
                  double *tcutxA, double *hxA, long *pncutxA, double *tdA,
                  double *tcutxB, double *hxB, long *pncutxB, double *tdB,
                  long *code, double *t)
{
  GetRNGstate();
  randhcdtl(pn,tcut,h,pncut,tcutxA,hxA,pncutxA,tdA,tcutxB,hxB,pncutxB,tdB,
            code,t);
  PutRNGstate();
}
