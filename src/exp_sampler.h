#ifndef EXP_SAMPLER_H
#define EXP_SAMPLER_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "clustering.h"
#include "gibbs.h"

// Exp(theta) with a Gamma(a,b) prior on theta.
class Exp_Params {
    public:
        double a, b;
};
class Exp_Cluster {
        int n; // number of elements in the cluster
        double sum_x;
        bool initialized;
        double a, b;

    public:
        // Must use a default constructor such as this (taking no arguments), in order for sets of clusters to work properly
        Exp_Cluster() { initialized = false; } // DO NOT MODIFY TO TAKE ARGUMENTS
        int size() { return n; }

        void reset() {
            sum_x = 0;
        }
        void initialize(Exp_Params p) {
            a = p.a; b = p.b;  // prior parameters
            n = 0;
            sum_x = 0;
            initialized = true;
        }
        void insert(double x) {
            if (!initialized) stop("(Cluster error) Attempting to insert item into uninitialized cluster.");
            n += 1;
            sum_x += x;
        }
        void remove(double x) {
            if (!initialized) stop("(Cluster error) Attempting to remove item from uninitialized cluster.");
            n -= 1;
            sum_x -= x;
        }
        // log (marginal) likelihood
        double log_lik() {
            return a*log(b) - (a+n)*log(b+sum_x) - lgamma(a) + lgamma(a+n);
        }
        // log (marginal) likelihood if x were added to the cluster
        double log_lik_w(double x) {
            return a*log(b) - (a+n+1)*log(b+sum_x+x) - lgamma(a) + lgamma(a+n+1);
        }
};


// [[Rcpp::export]]
IntegerMatrix Exp_Sampler(NumericVector data_, NumericVector A, NumericVector B, NumericVector params, int n_samples, int spacing) {
    vector<double> data = as< vector<double> >(data_);
    Exp_Params p; p.a = params[0]; p.b = params[1];
    //Exp_Params p = {params[0],params[1]};
    Clustering<Exp_Cluster,Exp_Params,double> C(data,p,false);
    IntegerMatrix z = gibbs(C,A,B,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


/*** R
n <- 5
alpha <- 1
Exp_Sampler(rexp(n),seq(1,n),rep(alpha,n),c(1,1),10)
*/



#endif


