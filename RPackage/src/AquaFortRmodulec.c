#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL

void F77_NAME(xcorr2d_f)(int m, int n, int p, int q, int k, int l, double *a, double *b, double *cc);
void F77_NAME(conv2d_f)(int m, int n, int p, int q, int k, int l, double *a, double *b, double *conv);
void F77_NAME(cape_f)(double *t_parcel, double *dwpt_parcel, double *mr_parcel, int nlevels, double *p_profile,
                      double *t_profile, double *mr_profile, int vtc, int nresult, double *result);

extern SEXP c_xcorr2d_f(SEXP a, SEXP b)
{
    int m = Rf_nrows(a);
    int n = Rf_ncols(a);
    //
    int p = Rf_nrows(b);
    int q = Rf_ncols(b);
    //
    int k = m + p - 1;
    int l = n + q - 1;

    SEXP ret;
    PROTECT(ret = Rf_allocMatrix(REALSXP, k, l));
    F77_CALL(xcorr2d_f)
    (m, n, p, q, k, l, REAL(a), REAL(b), REAL(ret));
    UNPROTECT(1);
    return (ret);
}

extern SEXP c_conv2d_f(SEXP a, SEXP b)
{
    int m = Rf_nrows(a);
    int n = Rf_ncols(a);
    //
    int p = Rf_nrows(b);
    int q = Rf_ncols(b);
    //
    int k = m + p - 1;
    int l = n + q - 1;

    SEXP ret;
    PROTECT(ret = Rf_allocMatrix(REALSXP, k, l));
    F77_CALL(conv2d_f)
    (m, n, p, q, k, l, REAL(a), REAL(b), REAL(ret));
    UNPROTECT(1);
    return (ret);
}

extern SEXP c_cape_f(SEXP t_parcel, SEXP dwpt_parcel, SEXP mr_parcel,
                     SEXP p_profile, SEXP t_profile, SEXP mr_profile,
                     SEXP vtc)
{
    int nlevels = Rf_length(t_profile);
    int vtc2 = INTEGER(vtc)[0];
    int nret = 4;

    SEXP ret;
    /* Order: CAPE, CIN, p_LCL, p_LFC */
    PROTECT(ret = Rf_allocVector(REALSXP, 4));
    F77_CALL(cape_f)
    (REAL(t_parcel), REAL(dwpt_parcel), REAL(mr_parcel), nlevels,
     REAL(p_profile), REAL(t_profile), REAL(mr_profile), vtc2,
     nret, REAL(ret));
    UNPROTECT(1);
    return (ret);
}

static const R_CallMethodDef CallEntries[] = {
    {"c_xcorr2d_f", (DL_FUNC)&c_xcorr2d_f, 2},
    {"c_conv2d_f", (DL_FUNC)&c_conv2d_f, 2},
    {"c_cape_f", (DL_FUNC)&c_cape_f, 7},
    {NULL, NULL, 0}};


void R_init_AquaFortR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
