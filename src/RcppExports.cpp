// do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fisher_test
double fisher_test(double a, double b, double c, double d);
RcppExport SEXP GOpro_fisher_test(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(fisher_test(a, b, c, d));
    return rcpp_result_gen;
END_RCPP
}
