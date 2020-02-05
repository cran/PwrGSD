typedef void WtFun(double *time, int *nrisk, int *nevent, int *pntimes, double *par, double *wt);

typedef int CmprFun(const void *x, const void *y);

extern WtFun flemhar, sflemhar, ramp, *wtfun;

extern CmprFun compitea, CmprDbl, *cmprdbl, *f;

