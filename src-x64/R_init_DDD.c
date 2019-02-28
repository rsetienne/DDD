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

void R_init_DDD(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}