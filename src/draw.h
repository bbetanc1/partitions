#ifndef DRAW_H
#define DRAW_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "logsumexp.h"

int randp(vector<double> p) {
    double u = unif_rand();
    int i = -1;
    double c = 0;
    do {
        i += 1;
        c += p[i];
    } while (u > c);
    return i;
}

class Draw {
        double log_s;
        vector<double> p;
        vector<double> log_p;
        vector<int> xs;
        bool p_valid;

        void compute_p() {
            for (int i=0; i<int(log_p.size()); i++) p.push_back(exp(log_p[i]-log_s));
            p_valid = true;
        }
    public:
        Draw() {
            log_s = -INFINITY;
            p_valid = false;
        }
        void reset() {
            log_s = -INFINITY;
            p.resize(0);
            log_p.resize(0);
            xs.resize(0);
            p_valid = false;
        }
        void insert(double log_px, int x) {
            log_s = logsumexp(log_s, log_px);
            log_p.push_back(log_px);
            xs.push_back(x);
        }
        int sample() {
            if (!p_valid) compute_p();
            return xs[randp(p)];
        }
};


// [[Rcpp::export]]
IntegerVector draw(int n, NumericVector p, IntegerVector x) {
    GetRNGstate();
    int k = p.size();
    Draw D;
    IntegerVector y(n);
    for (int i=0; i<k; i++) D.insert(log(p[i]),x[i]);
    for (int i=0; i<n; i++) y[i] = D.sample();
    PutRNGstate();
    return y;
}

/*** R
draw(100,c(1,2,1),c(0,1,2))
*/

#endif

