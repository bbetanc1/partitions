#ifndef LOGLIKF_H
#define LOGLIKF_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

double logamaf(vector<double> a, vector<int> c){
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
double loglikf(double conc, IntegerMatrix x, IntegerVector z, List params) {
    int nrow = x.nrow(), ncol = x.ncol();
    
    vector< vector<int> > c; // counts for each field
    vector< vector<double> > a; // params for each field
    vector<int> fsize; //vector with sizes of each field
    
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
        for (int j=0; j<M; j++) { a[i][j] = f0[j]*conc;}
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
            lg += logamaf(a[i], c[i]);
        }
       
    }
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;

    return lg;
}

// [[Rcpp::export]]
double loglikds(NumericVector vconc, IntegerMatrix x, IntegerVector z,
                List params) {
    int nrow = x.nrow(), ncol = x.ncol();
    
    vector< vector<int> > c; // counts for each field
    vector< vector<double> > a; // params for each field
    vector<int> fsize; //vector with sizes of each field
    
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
        for (int j=0; j<M; j++) { a[i][j] = f0[j]*vconc[i];}
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
            lg += logamaf(a[i], c[i]);
        }
       
    }
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;

    return lg;
}

// [[Rcpp::export]]
NumericVector loglikdsV(NumericVector vconc, IntegerMatrix x, IntegerVector z,
                List params) {
    int nrow = x.nrow(), ncol = x.ncol();
    
    vector< vector<int> > c; // counts for each field
    vector< vector<double> > a; // params for each field
    vector<int> fsize; //vector with sizes of each field
    NumericVector lg(ncol);
    
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
        for (int j=0; j<M; j++) { a[i][j] = f0[j]*vconc[i];}
    }
    
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
        for (int j=0; j<L; j++){
            for (int i=0; i<fsize[j]; i++){
                c[j][i] = 0;
            }
            for (int k=0; k<K; k++) {
                for (int i=0; i<nrow; i++) {
                    if ((z[i]-1)==k){
                        c[j][x(i,j)-1] += 1;
                    }
                }
            lg(j) = lg(j) + logamaf(a[j], c[j]);
            }
        }
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;
    
    return lg;
}

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
    
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
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
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;
    
    return c;
}

double loglikcn(NumericVector vconc, IntegerMatrix x, IntegerVector z,
                List params,  vector<vector< vector<int> > >  ccs) {
    int ncol = x.ncol();
    
    vector< vector<double> > a; // params for each field
    vector<int> fsize; //vector with sizes of each field
    //vector<int> c;
    
    int K = max(z);
    int L = ncol;
    fsize.resize(L);
    a.resize(L);
    for (int i=0; i<L; i++) {
        NumericVector f0 = params[i];
        int M = f0.size();
        fsize[i] = M;
        a[i].assign(M,0);
        for (int j=0; j<M; j++) { a[i][j] = f0[j]*vconc[i];}
    }
    
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
    double lg=0.0;
    for (int k=0; k<K; k++) {
        for (int i=0; i<L; i++){
            //c.resize(fsize[i]);
            //for (int j=0; j<fsize[i]; j++){
              //  c[j] =  ccs[k][i][j];
            //}
            lg += logamaf(a[i], ccs[k][i]);
        }
    }
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;
    
    return lg;
}

//loglik for slice sampler of distortions (betas), did ti with separate counts (ccs)
// to not have to recompute them for each field
// This is wrong, usw loglikxSP for now
double loglikspb(NumericVector betas, IntegerMatrix x, IntegerVector z,
                List params,  vector<vector< vector<int> > >  ccs) {
    int ncol = x.ncol();
    int nrow = x.nrow();
    
    vector< vector<double> > a; // params for each field
    vector<int> fsize; //vector with sizes of each field
    //vector<int> c;
    
    int K = max(z);
    int L = ncol;
    fsize.resize(L);
    a.resize(L);
    double lb0=0;
    for (int i=0; i<L; i++) {
        lb0 += log(betas[i]);
        NumericVector f0 = params[i];
        int M = f0.size();
        fsize[i] = M;
        a[i].assign(M,0);
        for (int j=0; j<M; j++) { a[i][j] = f0[j];}
    }
    
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
    vector<double> sle;
    double lg=0.0;
    double lb1=0.0;
    
    for (int j=0; j<L; j++){
        for (int i=0; i<nrow; i++){
            lb1 += log(betas[j]) + log(a[j][x(i,j)-1]);
            
        }
            int M = a[j].size();
            sle.assign(M,0);
        for (int k=0; k<K; k++) {
            for (int m=0; m<M; m++){
                sle[m] = log(a[j][m]) + ccs[k][j][m]*(log(betas[j]*a[j][m] +
                                                          (1 - betas[j])) - log(betas[j]) - log(a[j][m]));

            }
            lg += logsumexpv(sle);
        }
    }
    
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;
    lg= lg + lb1;
    
    return lg;
}

double loglikspb1(double beta, int fd, IntegerVector z,
                 List params, int N, vector<vector< vector<int> > >  ccs) {
   
    int K = max(z);
    NumericVector a = params[fd];
    int M = a.size();
    
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
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
    
    //cout << "counts" << c[0][0] << " " << c[0][1] << " " << c[0][6]  << endl;
    //cout << "lgT" << lgT << endl;
    lg = lg + N*log(beta);
    
    return lg;
}



// [[Rcpp::export]]
double loglikone(double conc, int fd, IntegerMatrix x, IntegerVector z,
                 List params) {
  int nrow = x.nrow();
    
    vector<int>  c; // counts for each field
    vector<double> a; // params for each field
    
    int K = max(z);
   
    NumericVector f0 = params[fd];
    int M = f0.size();
    c.resize(M);
    a.resize(M);
    for (int j=0; j<M; j++) { a[j] = f0[j]*conc;}
         
    //cout << "aaaa" << a[0][0] << " " << a[0][1] << endl;
    
    double lg=0.0;
    for (int k=0; k<K; k++) {

         for (int i=0; i < M; i++){
              c[i] = 0;
         }
        
         for (int i=0; i<nrow; i++) {
              if ((z[i]-1)==k){
                  c[x(i,fd)-1] += 1;
              }
         }
        
         lg = logamaf(a, c);       
    }
    
    return lg;
}

#endif

