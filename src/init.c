#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GPareto_EHI_2d_wrap_Rcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP GPareto_exipsi_Rcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP GPareto_hvolume2d_Rcpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GPareto_EHI_2d_wrap_Rcpp", (DL_FUNC) &GPareto_EHI_2d_wrap_Rcpp, 4},
    {"GPareto_exipsi_Rcpp",      (DL_FUNC) &GPareto_exipsi_Rcpp,      4},
    {"GPareto_hvolume2d_Rcpp",   (DL_FUNC) &GPareto_hvolume2d_Rcpp,   3},
    {NULL, NULL, 0}
};

void R_init_GPareto(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
