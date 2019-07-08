#ifndef LOGLIKX_H
#define LOGLIKX_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;


double logama(vector<double> a, vector<int> c){
    double suma=0;
    double sumac=0;
    double sumlga=0;
    double sumlgac=0;
    int M = a.size();
    for(int i = 0; i < M; ++i) {
        suma += a[i];
        sumac += a[i] + c[i];
        sumlga += lgamma(a[i]);
        sumlgac += lgamma(a[i] + c[i]);
    }
    double lgsum = sumlgac - lgamma(sumac) + lgamma(suma) - sumlga;
    return lgsum;
}

// [[Rcpp::export]]
double loglikx(IntegerMatrix x, IntegerVector z, List params) {
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
    
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
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
        for (int i=0; i<L; i++){
            lg += logama(a[i], c[i]);
        }
       
    }
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;

    return lg;
}

// log-likehood for model in SteortsHallFienberg2016

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
    
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
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
            //cout << "betas=" << betas[j];
            int M = a[j].size();
            sle.assign(M,0);
            for (int m=0; m<M; m++){
                //cout << "C=" << c[j][k];
                sle[m] = log(a[j][m]) + c[j][m]*(log(betas[j]*a[j][m] +
                                                     (1 - betas[j])) - log(betas[j]) - log(a[j][m]));
            }
            lg += logsumexpv(sle);
        }
        
    }
    
    lg = lg + s;
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;
    
    return lg;
}

#endif

