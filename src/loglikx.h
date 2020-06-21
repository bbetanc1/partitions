#ifndef LOGLIKX_H
#define LOGLIKX_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// log-likehood for model in BetancourtZanellaSteorts2020

// [[Rcpp::export]]
double loglikxSP(NumericVector betas, IntegerMatrix x, IntegerVector z , List params) {
    int nrow = x.nrow(), ncol = x.ncol();
    
    vector< vector<int> > c; // counts for each field
    vector< vector<double> > a; // params for each field
    vector<int> fsize; //size of each field
    
    int K = max(z);
    int L = ncol;
    fsize.resize(L);
    c.resize(L);
    a.resize(L);
    for (int i=0; i<L; i++) {
        NumericVector f0 = params[i];
        int M = f0.size();
        fsize[i] = M;
        c[i].assign(M,0);
        a[i].assign(M,0);
        for (int j=0; j<M; j++) { a[i][j] = f0[j];}
    }
    
    double s=0;
    for (int j=0; j<L; j++){
        for (int i=0; i<nrow; i++) {
            s += log(betas[j]) + log(a[j][x(i,j)-1]);
        }
    }
            
    double lg=0.0;
    for (int k=0; k<K; k++) {
        for (int j=0; j<L; j++){
            for (int i=0; i<fsize[j]; i++){
                c[j][i] = 0;
            }
        }
        for (int j=0; j<L; j++){
            for (int i=0; i<nrow; i++) {
                if ((z[i]-1)==k){
                    c[j][x(i,j)-1] += 1;
                }
            }
        }
        
        vector<double> s0;
        vector<double> sle;
        s0.assign(L,0); // to save second term of f(x_c, betas, thetas)
        for (int j=0; j<L; j++){
            int M = a[j].size();
            sle.assign(M,0);
            for (int m=0; m<M; m++){
                sle[m] = log(a[j][m]) + c[j][m]*(log(betas[j]*a[j][m] +
                                                     (1 - betas[j])) - log(betas[j]) - log(a[j][m]));
            }
            lg += logsumexpv(sle);
        }
        
    }
    
    lg = lg + s;
    
    return lg;
}

#endif

