// [[Rcpp::depends(BH)]]
#include <boost/math/distributions/hypergeometric.hpp>
#include <algorithm>

using namespace boost::math;
using namespace std;

// [[register]]
// [[Rcpp::export]]
double fisher_test(double a, double b, double c, double d) {
  double N = a + b + c + d;
  double r = a + c;
  double n = c + d;
  double max_k = min(r, n);
  double min_k = max(0, int(r + n - N));
  //dist<> name(args)
  hypergeometric_distribution<> hgd(r, n, N);
  //double no_s = pdf(hgd, c);
  double tmp_p = 0.0;
  for(int k = min_k; k < max_k + 1; k++) {
    double p = pdf(hgd, k);
    if(k <= c) tmp_p = tmp_p + p;
  }
  return tmp_p;
}
