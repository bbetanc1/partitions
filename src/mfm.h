#ifndef MFM_H
#define MFM_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "logsumexp.h"

// Note: Do not export to Rcpp as is. To do so, it should be modified to take a 0-indexed vector log_pK.
NumericVector log_Vn_MFM(NumericVector log_pK1, double gamma, int n, int upto) {
    // Compute the logs of the coefficients V_n(t) for the MFM partition distribution.
    // 
    // INPUTS:
    //   log_pK1 = log probability mass function of K --- i.e., log_pK1[k] is log(P(K=k))
    //   gamma = finite-dimensional Dirichlet parameter
    //   n = sample size
    //   upto = largest value of t for which to compute V_n(t)
    //
    // OUTPUTS:
    //   log_Vn = array of log-coefficients --- i.e., log_Vn[t] is log(V_n(t))
    // 
    double tolerance = 1e-12;
    bool warn = false;
    NumericVector log_Vn(upto+1);
    for (int t=1; t<=upto; t++) {
        if (t > n) {
            log_Vn[t] = -INFINITY;
            continue;
        }
        double a = 0, b = 0, c = -INFINITY, p = 0; 
        int k = 1;
        while ((fabs(a-c) > tolerance) || (p < 1.0-tolerance)) {
            // Note: The first condition above is false when a and c are both -INFINITY.
            if (k >= log_pK1.size()) {
                if (p < 1.0-tolerance) warn = true;
                break;
            }
            if (k >= t) {
                a = c;
                b = lgamma(k+1) - lgamma(k-t+1) - lgamma(k*gamma+n) + lgamma(k*gamma) + log_pK1[k];
                c = logsumexp(a,b);
            }
            p += exp(log_pK1[k]);
            //Rprintf("t=%d k=%d p=%.4f a=%.3f b=%.3f c=%.3f\n",t,k,p,a,b,c);
            k += 1;
        }
        //Rprintf("\n");
        log_Vn[t] = c;
    }
    if (warn) Rf_warning("(MFM) Potential loss of accuracy due to insufficiently many values of log_pK provided.");
    return log_Vn;
}

// [[Rcpp::export]]
NumericVector a_MFM(double gamma, int n) {
    NumericVector a(n+1); // 1-based indexing
    for (int nc=1; nc<=n; nc++) a[nc] = nc+gamma;
    a.erase(0); // remove dummy value at 0th position
    return a;
}


// [[Rcpp::export]]
NumericVector b_MFM(NumericVector log_pK, double gamma, int n, int upto) {
    NumericVector log_pK1 = clone(log_pK);
    log_pK1.push_front(0); // insert dummy value at 0th position, to enable 1-based indexing
    log_pK1.push_back(-INFINITY); // insert log(0) at end, to give log_Vn_MFM enough padding
    NumericVector log_Vn = log_Vn_MFM(log_pK1,gamma,n,upto+1); // 1-based
    NumericVector b(upto+1); // 1-based
    for (int t=1; t<=upto; t++) {
        b[t] = exp(log_Vn[t+1] - log_Vn[t] + log(gamma));
    }
    b.erase(0); // remove dummy value at 0th position
    return b;
}

// Example usage:
/*** R
pK = rep(1/30, 30)  # uniform on {1,...,30}
pK = dpois(seq(0,100),5)  # K-1 ~ Poisson(5)
pK = dgeom(seq(0,1000),0.1)  # K-1 ~ Geometric(0.1)
b <- b_MFM(log(pK),1,100,50)
*/

#endif

