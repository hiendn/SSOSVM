// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// SquareHingeC
Rcpp::List SquareHingeC(arma::mat& YMAT, int DIM, double EPSILON, bool returnAll, bool SLOW);
RcppExport SEXP _SSOSVM_SquareHingeC(SEXP YMATSEXP, SEXP DIMSEXP, SEXP EPSILONSEXP, SEXP returnAllSEXP, SEXP SLOWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type YMAT(YMATSEXP);
    Rcpp::traits::input_parameter< int >::type DIM(DIMSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    Rcpp::traits::input_parameter< bool >::type returnAll(returnAllSEXP);
    Rcpp::traits::input_parameter< bool >::type SLOW(SLOWSEXP);
    rcpp_result_gen = Rcpp::wrap(SquareHingeC(YMAT, DIM, EPSILON, returnAll, SLOW));
    return rcpp_result_gen;
END_RCPP
}
// HingeC
Rcpp::List HingeC(arma::mat& YMAT, int DIM, double EPSILON, bool returnAll);
RcppExport SEXP _SSOSVM_HingeC(SEXP YMATSEXP, SEXP DIMSEXP, SEXP EPSILONSEXP, SEXP returnAllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type YMAT(YMATSEXP);
    Rcpp::traits::input_parameter< int >::type DIM(DIMSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    Rcpp::traits::input_parameter< bool >::type returnAll(returnAllSEXP);
    rcpp_result_gen = Rcpp::wrap(HingeC(YMAT, DIM, EPSILON, returnAll));
    return rcpp_result_gen;
END_RCPP
}
// LogisticC
Rcpp::List LogisticC(arma::mat& YMAT, int DIM, double EPSILON, bool returnAll, bool SLOW);
RcppExport SEXP _SSOSVM_LogisticC(SEXP YMATSEXP, SEXP DIMSEXP, SEXP EPSILONSEXP, SEXP returnAllSEXP, SEXP SLOWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type YMAT(YMATSEXP);
    Rcpp::traits::input_parameter< int >::type DIM(DIMSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    Rcpp::traits::input_parameter< bool >::type returnAll(returnAllSEXP);
    Rcpp::traits::input_parameter< bool >::type SLOW(SLOWSEXP);
    rcpp_result_gen = Rcpp::wrap(LogisticC(YMAT, DIM, EPSILON, returnAll, SLOW));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SSOSVM_SquareHingeC", (DL_FUNC) &_SSOSVM_SquareHingeC, 5},
    {"_SSOSVM_HingeC", (DL_FUNC) &_SSOSVM_HingeC, 4},
    {"_SSOSVM_LogisticC", (DL_FUNC) &_SSOSVM_LogisticC, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SSOSVM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
