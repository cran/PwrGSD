#include <R.h>

void agghaz(double *tagg, double *time, int *nrisk, int *nevent, int *pndth, int *pnb,
	    double *timea, int *nriska, int *neventa, int *pnagg)
{
/*
  Blocked survival data of form 'time', 'nrisk', 'nevent'
  with 'nrisk' and 'nevent' blocked into 'nb' columns 
  with a a row corresponding to an event from any of the blocks
  (zeros elsewhere) so that one common 'time' axis is possible. 

  This routine aggregates the events, event times and numbers
  at risk into windows of length 'tagg'.  In order that this
  has the desired impact, to stabilize answers at the end 
  of the follow-up period, blocking should start at the end, 
  with the left over odd piece at the beginning.
  
 */
    int i,jj,l,ndths,nagg,nb;
    int *DN;
    double t_max,t_old, DT;
    ndths = *pndth;
    nb = *pnb;
    DN = (int *)Calloc(nb, int);

    t_max = *(time + ndths - 1);
    nagg = (int) (floor(t_max/(*tagg)) + 1.0);
    *pnagg = nagg;
    l=0;
    t_old = t_max;
    for(jj=0;jj<nb;jj++) *(DN+jj) = 0;
    for (i=0;i<ndths;i++){
      DT = t_old - *(time + ndths-1-i);
      for(jj=0;jj<nb;jj++) *(DN+jj) = *(DN+jj) + *(nevent + ndths*jj + ndths-1-i);
      if(DT >= *tagg || i==(ndths-1)){
        *(timea + nagg-1-l) = *(time + ndths-1-i);
	t_old = *(time + ndths-1-i);
	for(jj=0;jj<nb;jj++){
	  *(nriska + nagg*jj + nagg-1-l) = *(nrisk + ndths*jj + ndths-1-i);
          *(neventa + nagg*jj + nagg-1-l) = *(DN + jj);
	  *(DN+jj) = 0;
	}
        l++;
      }
    }
    Free(DN);
}

