#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

void F77_NAME(xcorr2d_f)(int m, int n, int p, int q, int k, int l, double *a, double *b, double *cc);
void F77_NAME(conv2d_f)(int m, int n, int p, int q, int k, int l, double *a, double *b, double *conv);
void F77_NAME(cape_f)(double *t_parcel, double *dwpt_parcel, double *mr_parcel, int nlevels, double *p_profile,
                      double *t_profile, double *mr_profile, int vtc, int nresult, double *result);

extern SEXP c_xcorr2d_f(SEXP dim_a, SEXP a, SEXP dim_b, SEXP b)
{
    int m = INTEGER(dim_a)[0];
    int n = INTEGER(dim_a)[1];
    //
    int p = INTEGER(dim_b)[0];
    int q = INTEGER(dim_b)[1];
    //
    int k = m + p - 1;
    int l = n + q - 1;

    SEXP ret;
    PROTECT(ret = allocMatrix(REALSXP, k, l));
    F77_CALL(xcorr2d_f)
    (m, n, p, q, k, l, REAL(a), REAL(b), REAL(ret));
    UNPROTECT(1);
    return (ret);
}

extern SEXP c_conv2d_f(SEXP dim_a, SEXP a, SEXP dim_b, SEXP b)
{
    int m = INTEGER(dim_a)[0];
    int n = INTEGER(dim_a)[1];
    //
    int p = INTEGER(dim_b)[0];
    int q = INTEGER(dim_b)[1];
    //
    int k = m + p - 1;
    int l = n + q - 1;

    SEXP ret;
    PROTECT(ret = allocMatrix(REALSXP, k, l));
    F77_CALL(conv2d_f)
    (m, n, p, q, k, l, REAL(a), REAL(b), REAL(ret));
    UNPROTECT(1);
    return (ret);
}

extern SEXP c_cape_f(SEXP t_parcel, SEXP dwpt_parcel, SEXP mr_parcel,
                     SEXP p_profile, SEXP t_profile, SEXP mr_profile,
                     SEXP DIM)
{
    int nlevels = INTEGER(DIM)[0];
    int nret = INTEGER(DIM)[1];
    int vtc = INTEGER(DIM)[2];

    SEXP ret;
    /* Order: CAPE, CIN, p_LCL, p_LFC */
    PROTECT(ret = allocVector(REALSXP, 4));
    F77_CALL(cape_f)
    (REAL(t_parcel), REAL(dwpt_parcel), REAL(mr_parcel), nlevels,
     REAL(p_profile), REAL(t_profile), REAL(mr_profile), vtc,
     nret, REAL(ret));
    UNPROTECT(1);
    return (ret);
}

static const R_CallMethodDef CallEntries[] = {
    {"c_xcorr2d_f", (DL_FUNC)&c_xcorr2d_f, 4},
    {"c_conv2d_f", (DL_FUNC)&c_conv2d_f, 4},
    {"c_cape_f", (DL_FUNC)&c_cape_f, 7},
    {NULL, NULL, 0}};

static const R_FortranMethodDef FEntries[] = {
    {NULL, NULL, 0}};

void R_init_AquaFortR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
