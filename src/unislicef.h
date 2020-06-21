#ifndef SLICEF_H
#define SLICEF_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "loglikf.h"

// slice sampler for distortions (betas) in model with SteortsHallFienberg2016 likelihood

// [[Rcpp::export]]
NumericVector unislicespb1(NumericVector betas, IntegerMatrix x, IntegerVector z,
                         List params, NumericVector hpriords, double w, int m,
                         double lower, double upper, NumericVector x1, int N, int Khat,
                         IntegerVector Nk, NumericVector hpriorpar, std::string Prior) {
    
    int fields = x.ncol();
    NumericVector xx2 = betas;
    double cb = hpriords[0];
    double db = hpriords[1];
    //double lp = logpriorf(x1, N, Khat, Nk, hpriorpar, Prior); //constant
    vector<vector< vector<int> > > ccs = counts(x, z, params);
    
    for (int fd = 0; fd < fields; fd++) {
        
        double xx0 = xx2[fd];
        double gx0 = loglikspb1(xx0, fd, z, params, N, ccs) + (cb-1)*log(xx0) + (db-1)*log(1-xx0);
        
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
            {
                g = loglikspb1(L, fd, z, params, N, ccs) + (cb-1)*log(L) + (db-1)*log(1-L);
                if (g <= logy) break;
                L = L - w;
            }
            
            while(R < upper)
            {
                g = loglikspb1(R, fd, z, params, N, ccs) + (cb-1)*log(R) + (db-1)*log(1-R);
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
                g = loglikspb1(L, fd, z, params, N, ccs) + (cb-1)*log(L) + (db-1)*log(1-L);
                if (g <= logy) break;
                L = L - w;
                J = J - 1;
            }
            
            while (K>0)
            { if (R>=upper) break;
                g = loglikspb1(R, fd, z, params, N, ccs) + (cb-1)*log(R) + (db-1)*log(1-R);
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
            gx2 = loglikspb1(x2, fd, z, params, N, ccs) + (cb-1)*log(x2) + (db-1)*log(1-x2);
            
            if (gx2 < logy){
                if (x2>xx0)
                { R = x2;
                }
                else
                { L = x2;
                }
            }
        }while(gx2 < logy);
        
        xx2[fd] = x2;
        
    }   
    // Return new vector.
    return xx2; 
}

#endif

