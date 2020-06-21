#ifndef LOGLIKF_H
#define LOGLIKF_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

vector<vector< vector<int> > > counts(IntegerMatrix x, IntegerVector z, List params) {
    int nrow = x.nrow(), ncol = x.ncol();
    
    vector<vector< vector<int> > > c;
    vector< vector<double> > a; // params for each field
    vector<int> fsize; //vector with sizes of each field
    
    int K = max(z);
    int L = ncol;
    fsize.resize(L);

    // counts for each field
    for (int i=0; i<L; i++) {
        NumericVector f0 = params[i];
        int M = f0.size();
        fsize[i] = M;
    }
    
    c.resize(K);
    for (int i = 0; i < K; i++)
    {
        c[i].resize(L);
        for (int j = 0; j < L; j++)
        {
            c[i][j].resize(fsize[j]);
        }
    }
    
    for (int k=0; k<K; k++) {
        for (int j=0; j<L; j++){
            for (int i=0; i<fsize[j]; i++){
                c[k][j][i] = 0;
            }
        }
        for (int j=0; j<L; j++){
            for (int i=0; i<nrow; i++) {
                if ((z[i]-1)==k){
                    c[k][j][x(i,j)-1] += 1;
                }
            }
        }
    }
    
    return c;
}


double loglikspb1(double beta, int fd, IntegerVector z,
                 List params, int N, vector<vector< vector<int> > >  ccs) {
   
    int K = max(z);
    NumericVector a = params[fd];
    int M = a.size();
  
    vector<double> sle;
    double lg=0.0;
    for (int k=0; k<K; k++) {
            sle.assign(M,0);
            for (int m=0; m<M; m++){
                sle[m] = log(a[m]) + ccs[k][fd][m]*(log(beta*a[m] +
                                                          (1 - beta)) - log(beta) - log(a[m]));
                
            }
            lg += logsumexpv(sle);
    }
    
    lg = lg + N*log(beta);
    
    return lg;
}


#endif

