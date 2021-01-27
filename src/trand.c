#include<R.h>

void randfromh(int *pn, double *tcut, double *h, int *pncut, double *t);
void    randhcdtl(int *pn, double *tcut, double *h, int *pncut, double *tend,
	          double *tcutdA, double *hdA, int *pncutdA, double *tdA,
	          double *tcutdB, double *hdB, int *pncutdB, double *tdB,
                  double *tcutxA, double *hxA, int *pncutxA, double *tcutxB, 
		  double *hxB, int *pncutxB, int *gradual, int *code, double *t);

void trandfromh(int *pn, double *tcut, double *h, int *pncut, double *t)
{
  GetRNGstate();
  randfromh(pn, tcut,  h, pncut, t);
  PutRNGstate();
}

void trandhcdtl(int *pn, double *tcut, double *h, int *pncut, double *tend,
	        double *tcutdA, double *hdA, int *pncutdA, double *tdA,
	        double *tcutdB, double *hdB, int *pncutdB, double *tdB,
                double *tcutxA, double *hxA, int *pncutxA, double *tcutxB, 
		double *hxB, int *pncutxB, int *gradual, int *code, double *t)
{
  GetRNGstate();

  randhcdtl(pn, tcut, h, pncut, tend, tcutdA, hdA, pncutdA, tdA, 
            tcutdB, hdB, pncutdB, tdB, tcutxA, hxA, pncutxA, tcutxB, hxB, pncutxB,
	    gradual, code, t);

  PutRNGstate();
}
