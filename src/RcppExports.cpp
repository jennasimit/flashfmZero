// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// odds_no_sharing
double odds_no_sharing(const double kappa, const NumericVector p, const int ndis);
RcppExport SEXP _flashfmZero_odds_no_sharing(SEXP kappaSEXP, SEXP pSEXP, SEXP ndisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type ndis(ndisSEXP);
    rcpp_result_gen = Rcpp::wrap(odds_no_sharing(kappa, p, ndis));
    return rcpp_result_gen;
END_RCPP
}
// concat
String concat(CharacterVector x);
RcppExport SEXP _flashfmZero_concat(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(concat(x));
    return rcpp_result_gen;
END_RCPP
}
// calcDij
NumericVector calcDij(const int i, const int j, const NumericMatrix& Cr12, const NumericVector& Vr1, const NumericVector& Vr2);
RcppExport SEXP _flashfmZero_calcDij(SEXP iSEXP, SEXP jSEXP, SEXP Cr12SEXP, SEXP Vr1SEXP, SEXP Vr2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int >::type j(jSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Cr12(Cr12SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Vr1(Vr1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Vr2(Vr2SEXP);
    rcpp_result_gen = Rcpp::wrap(calcDij(i, j, Cr12, Vr1, Vr2));
    return rcpp_result_gen;
END_RCPP
}
// pairs
NumericMatrix pairs(const int M);
RcppExport SEXP _flashfmZero_pairs(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(pairs(M));
    return rcpp_result_gen;
END_RCPP
}
// indices
NumericVector indices(const int M);
RcppExport SEXP _flashfmZero_indices(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(indices(M));
    return rcpp_result_gen;
END_RCPP
}
// Mlogsum1
NumericMatrix Mlogsum1(const NumericMatrix& x);
RcppExport SEXP _flashfmZero_Mlogsum1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Mlogsum1(x));
    return rcpp_result_gen;
END_RCPP
}
// Mlogsum
double Mlogsum(const NumericMatrix& x);
RcppExport SEXP _flashfmZero_Mlogsum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Mlogsum(x));
    return rcpp_result_gen;
END_RCPP
}
// calcDelta
double calcDelta(const IntegerVector& imods, const int N, const int M, const double dcon, const List Vr, const List Cr, const NumericMatrix& c2);
RcppExport SEXP _flashfmZero_calcDelta(SEXP imodsSEXP, SEXP NSEXP, SEXP MSEXP, SEXP dconSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP c2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type imods(imodsSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type c2(c2SEXP);
    rcpp_result_gen = Rcpp::wrap(calcDelta(imods, N, M, dcon, Vr, Cr, c2));
    return rcpp_result_gen;
END_RCPP
}
// calcQD3
double calcQD3(const int mod1, const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const NumericMatrix& Nqq, const NumericVector& Ldcon12, const NumericMatrix& c2, const List lPP, const int Nsame);
RcppExport SEXP _flashfmZero_calcQD3(SEXP mod1SEXP, SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP NqqSEXP, SEXP Ldcon12SEXP, SEXP c2SEXP, SEXP lPPSEXP, SEXP NsameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type mod1(mod1SEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Nqq(NqqSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon12(Ldcon12SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    Rcpp::traits::input_parameter< const int >::type Nsame(NsameSEXP);
    rcpp_result_gen = Rcpp::wrap(calcQD3(mod1, N, nummods, Vr, Cr, dcon, keep, Nqq, Ldcon12, c2, lPP, Nsame));
    return rcpp_result_gen;
END_RCPP
}
// ppadjT3
NumericVector ppadjT3(const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const NumericMatrix& Nqq, const NumericVector& Ldcon12, const NumericMatrix& c2, const List lPP, const int Nsame);
RcppExport SEXP _flashfmZero_ppadjT3(SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP NqqSEXP, SEXP Ldcon12SEXP, SEXP c2SEXP, SEXP lPPSEXP, SEXP NsameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Nqq(NqqSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon12(Ldcon12SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    Rcpp::traits::input_parameter< const int >::type Nsame(NsameSEXP);
    rcpp_result_gen = Rcpp::wrap(ppadjT3(N, nummods, Vr, Cr, dcon, keep, Nqq, Ldcon12, c2, lPP, Nsame));
    return rcpp_result_gen;
END_RCPP
}
// calcQD4
double calcQD4(const int mod1, const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const NumericMatrix& Nqq, const NumericVector& Ldcon12, const NumericVector& Nq3, const NumericVector& Ldcon123, const List Cr3, const List lPP, const int Nsame);
RcppExport SEXP _flashfmZero_calcQD4(SEXP mod1SEXP, SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP NqqSEXP, SEXP Ldcon12SEXP, SEXP Nq3SEXP, SEXP Ldcon123SEXP, SEXP Cr3SEXP, SEXP lPPSEXP, SEXP NsameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type mod1(mod1SEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Nqq(NqqSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon12(Ldcon12SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Nq3(Nq3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon123(Ldcon123SEXP);
    Rcpp::traits::input_parameter< const List >::type Cr3(Cr3SEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    Rcpp::traits::input_parameter< const int >::type Nsame(NsameSEXP);
    rcpp_result_gen = Rcpp::wrap(calcQD4(mod1, N, nummods, Vr, Cr, dcon, keep, Nqq, Ldcon12, Nq3, Ldcon123, Cr3, lPP, Nsame));
    return rcpp_result_gen;
END_RCPP
}
// ppadjT4
NumericVector ppadjT4(const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const NumericMatrix& Nqq, const NumericVector& Ldcon12, const NumericVector& Nq3, const NumericVector& Ldcon123, const List Cr3, const List lPP, const int Nsame);
RcppExport SEXP _flashfmZero_ppadjT4(SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP NqqSEXP, SEXP Ldcon12SEXP, SEXP Nq3SEXP, SEXP Ldcon123SEXP, SEXP Cr3SEXP, SEXP lPPSEXP, SEXP NsameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Nqq(NqqSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon12(Ldcon12SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Nq3(Nq3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon123(Ldcon123SEXP);
    Rcpp::traits::input_parameter< const List >::type Cr3(Cr3SEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    Rcpp::traits::input_parameter< const int >::type Nsame(NsameSEXP);
    rcpp_result_gen = Rcpp::wrap(ppadjT4(N, nummods, Vr, Cr, dcon, keep, Nqq, Ldcon12, Nq3, Ldcon123, Cr3, lPP, Nsame));
    return rcpp_result_gen;
END_RCPP
}
// calcQD5
double calcQD5(const int mod1, const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const NumericMatrix& Nqq, const NumericVector& Ldcon12, const NumericVector& Nq3, const NumericVector& Ldcon123, const List Cr3, const NumericVector& Nq4, const NumericVector& Ldcon1234, const List Cr4, const List lPP, const int Nsame);
RcppExport SEXP _flashfmZero_calcQD5(SEXP mod1SEXP, SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP NqqSEXP, SEXP Ldcon12SEXP, SEXP Nq3SEXP, SEXP Ldcon123SEXP, SEXP Cr3SEXP, SEXP Nq4SEXP, SEXP Ldcon1234SEXP, SEXP Cr4SEXP, SEXP lPPSEXP, SEXP NsameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type mod1(mod1SEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Nqq(NqqSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon12(Ldcon12SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Nq3(Nq3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon123(Ldcon123SEXP);
    Rcpp::traits::input_parameter< const List >::type Cr3(Cr3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Nq4(Nq4SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon1234(Ldcon1234SEXP);
    Rcpp::traits::input_parameter< const List >::type Cr4(Cr4SEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    Rcpp::traits::input_parameter< const int >::type Nsame(NsameSEXP);
    rcpp_result_gen = Rcpp::wrap(calcQD5(mod1, N, nummods, Vr, Cr, dcon, keep, Nqq, Ldcon12, Nq3, Ldcon123, Cr3, Nq4, Ldcon1234, Cr4, lPP, Nsame));
    return rcpp_result_gen;
END_RCPP
}
// ppadjT5
NumericVector ppadjT5(const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const NumericMatrix& Nqq, const NumericVector& Ldcon12, const NumericVector& Nq3, const NumericVector& Ldcon123, const List Cr3, const NumericVector& Nq4, const NumericVector& Ldcon1234, const List Cr4, const List lPP, const int Nsame);
RcppExport SEXP _flashfmZero_ppadjT5(SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP NqqSEXP, SEXP Ldcon12SEXP, SEXP Nq3SEXP, SEXP Ldcon123SEXP, SEXP Cr3SEXP, SEXP Nq4SEXP, SEXP Ldcon1234SEXP, SEXP Cr4SEXP, SEXP lPPSEXP, SEXP NsameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Nqq(NqqSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon12(Ldcon12SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Nq3(Nq3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon123(Ldcon123SEXP);
    Rcpp::traits::input_parameter< const List >::type Cr3(Cr3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Nq4(Nq4SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Ldcon1234(Ldcon1234SEXP);
    Rcpp::traits::input_parameter< const List >::type Cr4(Cr4SEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    Rcpp::traits::input_parameter< const int >::type Nsame(NsameSEXP);
    rcpp_result_gen = Rcpp::wrap(ppadjT5(N, nummods, Vr, Cr, dcon, keep, Nqq, Ldcon12, Nq3, Ldcon123, Cr3, Nq4, Ldcon1234, Cr4, lPP, Nsame));
    return rcpp_result_gen;
END_RCPP
}
// calcQD6fast
double calcQD6fast(const int mod1, const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const List lPP);
RcppExport SEXP _flashfmZero_calcQD6fast(SEXP mod1SEXP, SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP lPPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type mod1(mod1SEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    rcpp_result_gen = Rcpp::wrap(calcQD6fast(mod1, N, nummods, Vr, Cr, dcon, keep, lPP));
    return rcpp_result_gen;
END_RCPP
}
// ppadjT6fast
NumericVector ppadjT6fast(const int N, const NumericVector& nummods, const List Vr, const List Cr, const double dcon, const List keep, const List lPP);
RcppExport SEXP _flashfmZero_ppadjT6fast(SEXP NSEXP, SEXP nummodsSEXP, SEXP VrSEXP, SEXP CrSEXP, SEXP dconSEXP, SEXP keepSEXP, SEXP lPPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nummods(nummodsSEXP);
    Rcpp::traits::input_parameter< const List >::type Vr(VrSEXP);
    Rcpp::traits::input_parameter< const List >::type Cr(CrSEXP);
    Rcpp::traits::input_parameter< const double >::type dcon(dconSEXP);
    Rcpp::traits::input_parameter< const List >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< const List >::type lPP(lPPSEXP);
    rcpp_result_gen = Rcpp::wrap(ppadjT6fast(N, nummods, Vr, Cr, dcon, keep, lPP));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flashfmZero_odds_no_sharing", (DL_FUNC) &_flashfmZero_odds_no_sharing, 3},
    {"_flashfmZero_concat", (DL_FUNC) &_flashfmZero_concat, 1},
    {"_flashfmZero_calcDij", (DL_FUNC) &_flashfmZero_calcDij, 5},
    {"_flashfmZero_pairs", (DL_FUNC) &_flashfmZero_pairs, 1},
    {"_flashfmZero_indices", (DL_FUNC) &_flashfmZero_indices, 1},
    {"_flashfmZero_Mlogsum1", (DL_FUNC) &_flashfmZero_Mlogsum1, 1},
    {"_flashfmZero_Mlogsum", (DL_FUNC) &_flashfmZero_Mlogsum, 1},
    {"_flashfmZero_calcDelta", (DL_FUNC) &_flashfmZero_calcDelta, 7},
    {"_flashfmZero_calcQD3", (DL_FUNC) &_flashfmZero_calcQD3, 12},
    {"_flashfmZero_ppadjT3", (DL_FUNC) &_flashfmZero_ppadjT3, 11},
    {"_flashfmZero_calcQD4", (DL_FUNC) &_flashfmZero_calcQD4, 14},
    {"_flashfmZero_ppadjT4", (DL_FUNC) &_flashfmZero_ppadjT4, 13},
    {"_flashfmZero_calcQD5", (DL_FUNC) &_flashfmZero_calcQD5, 17},
    {"_flashfmZero_ppadjT5", (DL_FUNC) &_flashfmZero_ppadjT5, 16},
    {"_flashfmZero_calcQD6fast", (DL_FUNC) &_flashfmZero_calcQD6fast, 8},
    {"_flashfmZero_ppadjT6fast", (DL_FUNC) &_flashfmZero_ppadjT6fast, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_flashfmZero(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
