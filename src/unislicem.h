
#ifndef MSLICE_H
#define MSLICE_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "logprior.h"

// [[Rcpp::export]]
NumericVector unislicem(NumericVector x1, int N, int Khat, double lx,
                IntegerVector Nk, NumericVector hpriorpar, double w, int m, NumericVector lo,
                NumericVector up, std::string Prior, IntegerVector samind) {
    NumericVector xx2 = x1;
    
    //cout << "size=" << x1.size() << endl;
    for (int k0 = 0; k0<x1.size(); k0++) {
        if (samind[k0]==1){
            double lower = lo[k0];
            double upper = up[k0];
            double xx0 = xx2[k0];
            
            // cout << "low-up" << lower << " "<< upper <<endl;
            
            double gx0 = logprior(xx0, k0+1, N, Khat, xx2,
                          Nk, hpriorpar, Prior) + lx;
            
            //cout << "check1-uni" << endl;
            
            double logy = gx0 - rexp(1)[0];
            
            // Find the initial interval to sample from.
            
            double u = runif(1,0,w)[0];
            double L = xx0 - u;
            double R = xx0 + (w-u);  // should guarantee that x0 is in [L,R], even with roundoff
            double g;
            
            // Expand the interval until its ends are outside the slice, or until
            // the limit on steps is reached.
            
            if (m==0)  // no limit on number of steps
            {
                //cout << "check-uni" << endl;
                while(L > lower)
                {   g = logprior(L, k0+1, N, Khat, xx2, Nk,
                        hpriorpar, Prior) + lx;
                    if (g <= logy) break;
                    L = L - w;
                }
                
                while(R < upper)
                {   g = logprior(R, k0+1, N, Khat, xx2, Nk,
                        hpriorpar, Prior) + lx;
                    if (g <= logy) break;
                    R = R + w;
                }
                
            }
            
            else if (m>1)  // limit on steps, bigger than one
            {
                int J = floor(runif(1,0,m)[0]);
                int K = (m-1) - J;
                
                while (J>0)
                { if (L<=lower) break;
                    g = logprior(L, k0+1, N, Khat, xx2, Nk,
                        hpriorpar, Prior) + lx;
                    if (g <= logy) break;
                    L = L - w;
                    J = J - 1;
                }
                
                while (K>0)
                { if (R>=upper) break;
                    g = logprior(R, k0+1, N, Khat, xx2, Nk,
                        hpriorpar, Prior) + lx;
                    if (g <= logy) break;
                    R = R + w;
                    K = K - 1;
                }
            }
            
            // Shrink interval to lower and upper bounds.
            
            if (L<lower)
            { L = lower;
            }
            if (R>upper)
            { R = upper;
            }
            
            // Sample from the interval, shrinking it on each rejection.
            double x2;
            double gx2;
            
            do
            {
                x2 = runif(1,L,R)[0];
                
                gx2 = logprior(x2, k0+1, N, Khat, xx2, Nk,
                      hpriorpar, Prior) + lx;
                
                if (gx2 < logy){
                    if (x2>xx0)
                    { R = x2;
                    }
                    else
                    { L = x2;
                    }
                }
            }while(gx2 < logy);
            
            xx2[k0] = x2;
            
        }
        
    }
     // Return new vector.
    
    return xx2;
}


// [[Rcpp::export]]
NumericVector unislicemESCD(NumericVector x1, double lx, IntegerVector Lm, NumericVector mu0, NumericVector hpriorpar, double w, int m, NumericVector lo, NumericVector up, IntegerVector samind) {
    NumericVector xx2 = x1;
    
    //cout << "size=" << x1.size() << endl;
    for (int k0 = 0; k0<x1.size(); k0++) {
        if (samind[k0]==1){
            double lower = lo[k0];
            double upper = up[k0];
            double xx0 = xx2[k0];
            
            // cout << "low-up" << lower << " "<< upper <<endl;
            
            double gx0 = logpriorESCD(xx0, k0+1, xx2, Lm, mu0, hpriorpar) + lx;
            
            //cout << "check1-uni" << endl;
            
            double logy = gx0 - rexp(1)[0];
            
            // Find the initial interval to sample from.
            
            double u = runif(1,0,w)[0];
            double L = xx0 - u;
            double R = xx0 + (w-u);  // should guarantee that x0 is in [L,R], even with roundoff
            double g;
            
            // Expand the interval until its ends are outside the slice, or until
            // the limit on steps is reached.
            
            if (m==0)  // no limit on number of steps
            {
                //cout << "check-uni" << endl;
                while(L > lower)
                {   g = logpriorESCD(L, k0+1, xx2, Lm, mu0, hpriorpar) +lx;
                    
                    if (g <= logy) break;
                    L = L - w;
                }
                
                while(R < upper)
                {   g = logpriorESCD(R, k0+1, xx2, Lm, mu0, hpriorpar) + lx;
                    if (g <= logy) break;
                    R = R + w;
                }
                
            }
            
            else if (m>1)  // limit on steps, bigger than one
            {
                int J = floor(runif(1,0,m)[0]);
                int K = (m-1) - J;
                
                while (J>0)
                { if (L<=lower) break;
                    g = logpriorESCD(L, k0+1, xx2, Lm, mu0, hpriorpar) + lx;
                    if (g <= logy) break;
                    L = L - w;
                    J = J - 1;
                }
                
                while (K>0)
                { if (R>=upper) break;
                    g = logpriorESCD(R, k0+1, xx2, Lm, mu0, hpriorpar) + lx;
                    if (g <= logy) break;
                    R = R + w;
                    K = K - 1;
                }
            }
            
            // Shrink interval to lower and upper bounds.
            
            if (L<lower)
            { L = lower;
            }
            if (R>upper)
            { R = upper;
            }
            
            // Sample from the interval, shrinking it on each rejection.
            double x2;
            double gx2;
            
            do
            {
                x2 = runif(1,L,R)[0];
                
                gx2 = logpriorESCD(x2, k0+1, xx2, Lm, mu0, hpriorpar) + lx;
                
                if (gx2 < logy){
                    if (x2>xx0)
                    { R = x2;
                    }
                    else
                    { L = x2;
                    }
                }
            }while(gx2 < logy);
            
            xx2[k0] = x2;
            
        }
        
    }
    // Return new vector.
    
    return xx2;
}

#endif

