#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(dd_fill1d)(double *vec, int *DIMP, double *parms, int *II);
extern void F77_NAME(dd_initmod)(void (*steadyparms)(int *, double *));
extern void F77_NAME(dd_runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(dd_runmodbw)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(dd_runmodtd)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);

static const R_FortranMethodDef FortranEntries[] = {
  {"dd_fill1d", (DL_FUNC) &F77_NAME(dd_fill1d),  4},
  {"dd_initmod", (DL_FUNC) &F77_NAME(dd_initmod),  1},
  {"dd_runmod", (DL_FUNC) &F77_NAME(dd_runmod),  6},
  {"dd_runmodbw", (DL_FUNC) &F77_NAME(dd_runmodbw),  6},
  {"dd_runmodtd", (DL_FUNC) &F77_NAME(dd_runmodtd),  6},
  {NULL, NULL, 0}
};


/* C bindings */
extern SEXP dd_integrate_odeint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dd_integrate_bw_odeint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dd_logliknorm1_odeint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dd_logliknorm2_odeint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"dd_integrate_odeint", (DL_FUNC) &dd_integrate_odeint, 6},
  {"dd_integrate_bw_odeint", (DL_FUNC) &dd_integrate_bw_odeint, 6},
  {"dd_logliknorm1_odeint", (DL_FUNC) &dd_logliknorm1_odeint, 6},
  {"dd_logliknorm2_odeint", (DL_FUNC) &dd_logliknorm2_odeint, 6},
  {NULL, NULL, 0}
};


void R_init_DDD(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
