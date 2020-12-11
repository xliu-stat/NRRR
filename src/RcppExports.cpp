// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// pinv_my
Rcpp::List pinv_my(arma::mat Bg);
RcppExport SEXP _NRRR_pinv_my(SEXP BgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Bg(BgSEXP);
    rcpp_result_gen = Rcpp::wrap(pinv_my(Bg));
    return rcpp_result_gen;
END_RCPP
}
// rrr_my
Rcpp::List rrr_my(arma::mat X, arma::mat Y, int r);
RcppExport SEXP _NRRR_rrr_my(SEXP XSEXP, SEXP YSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(rrr_my(X, Y, r));
    return rcpp_result_gen;
END_RCPP
}
// BLCD
Rcpp::List BLCD(arma::mat X, arma::mat Y, arma::mat Xl0, arma::mat Yl0, arma::mat Ag0, arma::mat Bg0, arma::mat Al0, arma::mat Bl0, int n, int r, int d, int p, int rx, int ry, int jx, int jy, int maxiter, double conv);
RcppExport SEXP _NRRR_BLCD(SEXP XSEXP, SEXP YSEXP, SEXP Xl0SEXP, SEXP Yl0SEXP, SEXP Ag0SEXP, SEXP Bg0SEXP, SEXP Al0SEXP, SEXP Bl0SEXP, SEXP nSEXP, SEXP rSEXP, SEXP dSEXP, SEXP pSEXP, SEXP rxSEXP, SEXP rySEXP, SEXP jxSEXP, SEXP jySEXP, SEXP maxiterSEXP, SEXP convSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xl0(Xl0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Yl0(Yl0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ag0(Ag0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Bg0(Bg0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Al0(Al0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Bl0(Bl0SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type rx(rxSEXP);
    Rcpp::traits::input_parameter< int >::type ry(rySEXP);
    Rcpp::traits::input_parameter< int >::type jx(jxSEXP);
    Rcpp::traits::input_parameter< int >::type jy(jySEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    rcpp_result_gen = Rcpp::wrap(BLCD(X, Y, Xl0, Yl0, Ag0, Bg0, Al0, Bl0, n, r, d, p, rx, ry, jx, jy, maxiter, conv));
    return rcpp_result_gen;
END_RCPP
}
// nrrr_init_my
Rcpp::List nrrr_init_my(arma::mat Y, arma::mat X, int r, int rx, int ry, int jx, int jy, int p, int d, int n);
RcppExport SEXP _NRRR_nrrr_init_my(SEXP YSEXP, SEXP XSEXP, SEXP rSEXP, SEXP rxSEXP, SEXP rySEXP, SEXP jxSEXP, SEXP jySEXP, SEXP pSEXP, SEXP dSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type rx(rxSEXP);
    Rcpp::traits::input_parameter< int >::type ry(rySEXP);
    Rcpp::traits::input_parameter< int >::type jx(jxSEXP);
    Rcpp::traits::input_parameter< int >::type jy(jySEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(nrrr_init_my(Y, X, r, rx, ry, jx, jy, p, d, n));
    return rcpp_result_gen;
END_RCPP
}
// nrrr_est_my
Rcpp::List nrrr_est_my(arma::mat Y, arma::mat X, int rini, int r, int rx, int ry, int jx, int jy, int p, int d, int n, int maxiter, double conv, int method, double lambda);
RcppExport SEXP _NRRR_nrrr_est_my(SEXP YSEXP, SEXP XSEXP, SEXP riniSEXP, SEXP rSEXP, SEXP rxSEXP, SEXP rySEXP, SEXP jxSEXP, SEXP jySEXP, SEXP pSEXP, SEXP dSEXP, SEXP nSEXP, SEXP maxiterSEXP, SEXP convSEXP, SEXP methodSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type rini(riniSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type rx(rxSEXP);
    Rcpp::traits::input_parameter< int >::type ry(rySEXP);
    Rcpp::traits::input_parameter< int >::type jx(jxSEXP);
    Rcpp::traits::input_parameter< int >::type jy(jySEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(nrrr_est_my(Y, X, rini, r, rx, ry, jx, jy, p, d, n, maxiter, conv, method, lambda));
    return rcpp_result_gen;
END_RCPP
}
// nrrr_select_my
Rcpp::List nrrr_select_my(arma::mat Y, arma::mat X, int xr, int rfit, int jx, int jy, int p, int d, int n, int ic, int maxiter, double conv, int method, double lambda, int dimred1, int dimred2, int dimred3);
RcppExport SEXP _NRRR_nrrr_select_my(SEXP YSEXP, SEXP XSEXP, SEXP xrSEXP, SEXP rfitSEXP, SEXP jxSEXP, SEXP jySEXP, SEXP pSEXP, SEXP dSEXP, SEXP nSEXP, SEXP icSEXP, SEXP maxiterSEXP, SEXP convSEXP, SEXP methodSEXP, SEXP lambdaSEXP, SEXP dimred1SEXP, SEXP dimred2SEXP, SEXP dimred3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< int >::type rfit(rfitSEXP);
    Rcpp::traits::input_parameter< int >::type jx(jxSEXP);
    Rcpp::traits::input_parameter< int >::type jy(jySEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type ic(icSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type dimred1(dimred1SEXP);
    Rcpp::traits::input_parameter< int >::type dimred2(dimred2SEXP);
    Rcpp::traits::input_parameter< int >::type dimred3(dimred3SEXP);
    rcpp_result_gen = Rcpp::wrap(nrrr_select_my(Y, X, xr, rfit, jx, jy, p, d, n, ic, maxiter, conv, method, lambda, dimred1, dimred2, dimred3));
    return rcpp_result_gen;
END_RCPP
}
// del_rows
arma::mat del_rows(arma::mat X, arma::uvec e);
RcppExport SEXP _NRRR_del_rows(SEXP XSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(del_rows(X, e));
    return rcpp_result_gen;
END_RCPP
}
// nrrr_cv_my
Rcpp::List nrrr_cv_my(arma::mat Y, arma::mat X, arma::uvec norder, int nfold, int xr, int rfit, int xrankfix, int yrankfix, int jx, int jy, int p, int d, int n, int ic, int maxiter, int method, int dimred1, int dimred2, int dimred3, double conv, double lambda);
RcppExport SEXP _NRRR_nrrr_cv_my(SEXP YSEXP, SEXP XSEXP, SEXP norderSEXP, SEXP nfoldSEXP, SEXP xrSEXP, SEXP rfitSEXP, SEXP xrankfixSEXP, SEXP yrankfixSEXP, SEXP jxSEXP, SEXP jySEXP, SEXP pSEXP, SEXP dSEXP, SEXP nSEXP, SEXP icSEXP, SEXP maxiterSEXP, SEXP methodSEXP, SEXP dimred1SEXP, SEXP dimred2SEXP, SEXP dimred3SEXP, SEXP convSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type norder(norderSEXP);
    Rcpp::traits::input_parameter< int >::type nfold(nfoldSEXP);
    Rcpp::traits::input_parameter< int >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< int >::type rfit(rfitSEXP);
    Rcpp::traits::input_parameter< int >::type xrankfix(xrankfixSEXP);
    Rcpp::traits::input_parameter< int >::type yrankfix(yrankfixSEXP);
    Rcpp::traits::input_parameter< int >::type jx(jxSEXP);
    Rcpp::traits::input_parameter< int >::type jy(jySEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type ic(icSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type dimred1(dimred1SEXP);
    Rcpp::traits::input_parameter< int >::type dimred2(dimred2SEXP);
    Rcpp::traits::input_parameter< int >::type dimred3(dimred3SEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(nrrr_cv_my(Y, X, norder, nfold, xr, rfit, xrankfix, yrankfix, jx, jy, p, d, n, ic, maxiter, method, dimred1, dimred2, dimred3, conv, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NRRR_pinv_my", (DL_FUNC) &_NRRR_pinv_my, 1},
    {"_NRRR_rrr_my", (DL_FUNC) &_NRRR_rrr_my, 3},
    {"_NRRR_BLCD", (DL_FUNC) &_NRRR_BLCD, 18},
    {"_NRRR_nrrr_init_my", (DL_FUNC) &_NRRR_nrrr_init_my, 10},
    {"_NRRR_nrrr_est_my", (DL_FUNC) &_NRRR_nrrr_est_my, 15},
    {"_NRRR_nrrr_select_my", (DL_FUNC) &_NRRR_nrrr_select_my, 17},
    {"_NRRR_del_rows", (DL_FUNC) &_NRRR_del_rows, 2},
    {"_NRRR_nrrr_cv_my", (DL_FUNC) &_NRRR_nrrr_cv_my, 21},
    {NULL, NULL, 0}
};

RcppExport void R_init_NRRR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
