#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(kdwd)(double *qval, double *Kmat, int *nobs, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlam, int *npass, int *jerr, double *alpmat);
extern void F77_NAME(kqdwd)(double *qvals, int *qlen, double *Kmat, int *nobs, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlams, int *npasses, int *jerrs, double *alpmats);
extern void F77_NAME(ldwd)(double *qval, double *Xmat, int *nobs, int *np, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlam, int *npass, int *jerr, double *btmat);
extern void F77_NAME(lqdwd)(double *qvals, int *qlen, double *Xmat, int *nobs, int *np, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlams, int *npasses, int *jerrs, double *btmats);
extern void F77_NAME(wkdwd)(double *qval, double *Kmat, double *wts, int *nobs, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlam, int *npass, int *jerr, double *alpmat);
extern void F77_NAME(wldwd)(double *qval, double *Xmat, double *wts, int *nobs, int *np, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlam, int *npass, int *jerr, double *btmat);

static const R_FortranMethodDef FortranEntries[] = {
    {"kdwd",  (DL_FUNC) &F77_NAME(kdwd),  13},
    {"kqdwd", (DL_FUNC) &F77_NAME(kqdwd), 14},
    {"ldwd",  (DL_FUNC) &F77_NAME(ldwd),  14},
    {"lqdwd", (DL_FUNC) &F77_NAME(lqdwd), 15},
    {"wkdwd", (DL_FUNC) &F77_NAME(wkdwd), 14},
    {"wldwd", (DL_FUNC) &F77_NAME(wldwd), 15},
    {NULL, NULL, 0}
};

void R_init_kerndwd(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
