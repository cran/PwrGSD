#include<R.h>
#include<Rmath.h>
#define MIN(x,y) (x<y ? x : y)
#define MAX(x,y) (x>y ? x : y)

/* arguments are expected to be of the following sizes:                */
/*---------------------------------------------------------------------*/
/* pmu    = double(2)      (mufu_k, mufu_end)                          */
/* pfrac  = double(1)      current pinffrac (variance info) vector     */
/* palpha = double(2)      alphatot vector (efficacy & futility)       */
/* psided = double(1)      -1,1, or 2 : sign and sided-ness of test    */ 
/* prho   = double(1)      the SC stopping criterion                   */
/* pef    = integer(1)     is this an efficacy (0) or futility (1) bnd */
/* b      = double(1)      output: current boundary point              */

/* Pr{ X(1) > b_m | X(f_k) } = 1-pnorm((b_m - X(fk))/(1-f_k)^0.5) */
/* Pr{ X(1) < a_m | X(f_k) } = pnorm((a_m - X(fk) - (mu(t_m)-mu(t_k)))/(1-f_k)^0.5) */

/* Pr{  Pr{ X(1) > b_m | X(f_l) } < rho, l=1,..,k-1, Pr{ X(1) > b_m | X(f_k) } > rho } 
   = Pr{ (b_m - X(f_l))/(1-f_l)^0.5 > qnorm(1-rho) , l=1,...,k-1,   (b_m - X(f_k))/(1-f_k)^0.5 < qnorm(1-rho) } 
   = Pr{ X(f_l) < b_m - (1-f_l)^0.5 *qnorm(1-rho) , l=1,...,k-1,  X(f_k) > b_m - (1-f_k)^0.5 *qnorm(1-rho) }

so b_k = (b_m - (1-f_k)^0.5 *qnorm(1-rho))/f_k^0.5
*/

/* Pr{  Pr{ X(1) < a_m | X(f_l) } < rho, l=1,..,k-1, Pr{ X(1) < a_m | X(f_k) } > rho } 
   = Pr{ (a_m - X(f_l) - (mu(t_m)-mu(t_l)))/(1-f_l)^0.5 < qnorm(rho),l=1,...,k-1,
          (a_m - X(f_k) - (mu(t_m)-mu(t_k)))/(1-f_k)^0.5 > qnorm(rho)             }

   = Pr{ X(f_l) > a_m - (mu(t_m)-mu(t_l)) - (1-f_l)^0.5 *qnorm(rho),l=1,...,k-1,
         X(f_k) < a_m - (mu(t_m)-mu(t_k)) - (1-f_k)^0.5 *qnorm(rho)              }

so a_k = (a_m - (mu(t_m)-mu(t_k)) - (1-f_k)^0.5 *qnorm(rho))/f_k^0.5

*/

void StCu2Bnds(double *pmu, double *pfrac, double *pzcrit, double *prho, int *pef, double *b)
{
  int ef, k;
  double rho, be_end, mu_end, f_k, mu_k;

  ef = *pef;
  rho = *prho;
  be_end = *pzcrit;
  mu_end = *(pmu + 1);

  f_k = *pfrac;
  f_k = MIN(MAX(f_k, 1.0e-16), 1.0-1.0e-16);
  mu_k = *pmu;
  if(ef==0)
    *b = (be_end-pow(1.0-f_k,0.5)*qnorm5(1.0-rho,0.0,1.0,1,0))/pow(f_k,0.5);
  if(ef==1)
    *b = (be_end-pow(1.0-f_k,0.5)*qnorm5(rho,0.0,1.0,1,0)-(mu_end-mu_k))/pow(f_k,0.5);
}
