#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP cooksdObs(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cooksdSubset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP covratioCalc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP covtraceCalc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cxxmatsub(SEXP, SEXP);
extern SEXP mdffitsSubset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"cooksdObs",     (DL_FUNC) &cooksdObs,     5},
  {"cooksdSubset",  (DL_FUNC) &cooksdSubset,  6},
  {"covratioCalc",  (DL_FUNC) &covratioCalc,  5},
  {"covtraceCalc",  (DL_FUNC) &covtraceCalc,  5},
  {"cxxmatsub",     (DL_FUNC) &cxxmatsub,     2},
  {"mdffitsSubset", (DL_FUNC) &mdffitsSubset, 6},
  {NULL, NULL, 0}
};

void R_init_HLMdiag(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
