#ifndef SLICEF_H
#define SLICEF_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "logpriorf.h"
#include "loglikf.h"

// [[Rcpp::export]]
double unislicef(double conc, IntegerMatrix x, IntegerVector z, List params,
                NumericVector hpriorconc, double w, int m, double lower,
                 double upper, NumericVector x1, int N, int Khat,
                 IntegerVector Nk, NumericVector hpriorpar, std::string Prior) {
    
    //cout << "size=" << x1.size() << endl;
    double xx0 = conc;
    double cg = hpriorconc[0];
    double dg = hpriorconc[1];
    double lp = logpriorf(x1, N, Khat, Nk, hpriorpar, Prior); //constant
    
    // cout << "low-up" << lower << " "<< upper <<endl;
    double gx0 = loglikf(xx0, x, z, params) + (cg-1)*log(xx0) - dg*xx0 + lp;
   
    
    //cout << "check1-uni" << endl;
    
    double logy = gx0 - rexp(1)[0];
    
    // Find the initial interval to sample from.
    
    double u = runif(1,0,w)[0];
    double L = xx0 - u;
    double R = xx0 + (w-u);  // should guarantee that x0 is in [L,R], even with roundoff
    double g;
    
    // Expand the interval until its ends are outside the slice, or until
    // the limit on steps is reached.
    
    
    //cout << "checkM" << " " << m<< " " << isinf(m/1.0) << endl;
    if (m==0)  // no limit on number of steps
    {
        //cout << "check-uni" << endl;
        while(L > lower)
        {   g = loglikf(L, x, z, params) + (cg-1)*log(L) - dg*L + lp;
            if (g <= logy) break;
            L = L - w;
        }
        
        while(R < upper)
        {   g = loglikf(R, x, z, params) + (cg-1)*log(R) - dg*R + lp;
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
            g = loglikf(L, x, z, params) + (cg-1)*log(L) - dg*L + lp;
            if (g <= logy) break;
            L = L - w;
            J = J - 1;
        }
        
        while (K>0)
        { if (R>=upper) break;
            g = loglikf(R, x, z, params) + (cg-1)*log(R) - dg*R + lp;
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
        
        gx2 = loglikf(x2, x, z, params) + (cg-1)*log(x2) - dg*x2 + lp;
        
        if (gx2 < logy){
            if (x2>xx0)
            { R = x2;
            }
            else
            { L = x2;
            }
        }
    }while(gx2 < logy);
    
     // Return new vector.
    
    return x2;
}


// slice sampler for distortions (deltas) in basic model with Dirichlte-Multinomial likelihood

// [[Rcpp::export]]
NumericVector unisliceds(NumericVector vconc, IntegerMatrix x, IntegerVector z,
                  List params, NumericVector hpriorconc, double w, int m,
                  double lower, double upper, NumericVector x1, int N, int Khat,
                 IntegerVector Nk, NumericVector hpriorpar, std::string Prior) {
    
    int fields = x.ncol();
    NumericVector xx2 = vconc;
    double cg = hpriorconc[0];
    double dg = hpriorconc[1];
    double lp = logpriorf(x1, N, Khat, Nk, hpriorpar, Prior); //constant
    vector<vector< vector<int> > > ccs = counts(x, z, params);

    for (int fd = 0; fd < fields; fd++) {

       double xx0 = xx2[fd];
       double gx0 = loglikcn(xx2, x, z, params, ccs) + //loglikds(xx2, x, z, params)
                    (cg-1)*log(xx0) - dg*xx0 + lp;
    
       //cout << "check1-uni" << endl;
    
       double logy = gx0 - rexp(1)[0];
    
       // Find the initial interval to sample from.
    
       double u = runif(1,0,w)[0];
       double L = xx0 - u;
       double R = xx0 + (w-u);  // should guarantee that x0 is in [L,R], even with roundoff
       double g;
    
       // Expand the interval until its ends are outside the slice, or until
       // the limit on steps is reached.
      
       //cout << "checkM" << " " << m<< " " << isinf(m/1.0) << endl;
       if (m==0)  // no limit on number of steps
       {
        //cout << "check-uni" << endl;
	  while(L > lower)
      {  xx2[fd] = L;
         g = loglikcn(xx2, x, z, params, ccs) + (cg-1)*log(L) - dg*L + lp; //loglikone(L, fd, x, z, params)
	     if (g <= logy) break;
	     L = L - w;
	  }
        
	  while(R < upper)
	  {  xx2[fd] = R;
         g = loglikcn(xx2, x, z, params, ccs) + (cg-1)*log(R) - dg*R + lp; //loglikone(R, fd, x, z, params)
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
         xx2[fd] = L;
	     g = loglikcn(xx2, x, z, params, ccs) + (cg-1)*log(L) - dg*L + lp;
	     if (g <= logy) break;
	     L = L - w;
	     J = J - 1;
	  }
        
	  while (K>0)
	  { if (R>=upper) break;
         xx2[fd] = R;
	     g = loglikcn(xx2, x, z, params, ccs) + (cg-1)*log(R) - dg*R + lp;
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
      xx2[fd] = x2;
	  gx2 = loglikcn(xx2, x, z, params, ccs)  + (cg-1)*log(x2) - dg*x2 + lp;
        
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
        
        //cout << "checkM" << " " << m << " " << isinf(m/1.0) << endl;
        if (m==0)  // no limit on number of steps
        {
            //cout << "check-uni" << endl;
            while(L > lower)
            {
                g = loglikspb1(L, fd, z, params, N, ccs) + (cb-1)*log(L) + (db-1)*log(1-L); //loglikone(L, fd, x, z, params)
                if (g <= logy) break;
                L = L - w;
            }
            
            while(R < upper)
            {
                g = loglikspb1(R, fd, z, params, N, ccs) + (cb-1)*log(R) + (db-1)*log(1-R); //loglikone(R, fd, x, z, params)
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

// [[Rcpp::export]]
NumericVector unislicespb(NumericVector betas, IntegerMatrix x, IntegerVector z,
                          List params, NumericVector hpriords, double w, int m,
                          double lower, double upper, NumericVector x1, int N, int Khat,
                          IntegerVector Nk, NumericVector hpriorpar, std::string Prior) {
    
    int fields = x.ncol();
    NumericVector xx2 = betas;
    double cb = hpriords[0];
    double db = hpriords[1];
    double lp = logpriorf(x1, N, Khat, Nk, hpriorpar, Prior); //constant
    
    for (int fd = 0; fd < fields; fd++) {
        
        double xx0 = xx2[fd];
        double gx0 = loglikxSP(xx2, x, z, params) + (cb-1)*log(xx0) + (db-1)*log(1-xx0) + lp;
        
        //cout << "check1-uni" << endl;
        
        double logy = gx0 - rexp(1)[0];
        
        // Find the initial interval to sample from.
        
        double u = runif(1,0,w)[0];
        double L = xx0 - u;
        double R = xx0 + (w-u);  // should guarantee that x0 is in [L,R], even with roundoff
        double g;
        
        // Expand the interval until its ends are outside the slice, or until
        // the limit on steps is reached.
        
        //cout << "checkM" << " " << m<< " " << isinf(m/1.0) << endl;
        if (m==0)  // no limit on number of steps
        {
            //cout << "check-uni" << endl;
            while(L > lower)
            {  xx2[fd] = L;
                g = loglikxSP(xx2, x, z, params) + (cb-1)*log(L) + (db-1)*log(1-L) + lp; //loglikone(L, fd, x, z, params)
                if (g <= logy) break;
                L = L - w;
            }
            
            while(R < upper)
            {  xx2[fd] = R;
                g = loglikxSP(xx2, x, z, params) + (cb-1)*log(R) + (db-1)*log(1-R) + lp; //loglikone(R, fd, x, z, params)
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
                xx2[fd] = L;
                g = loglikxSP(xx2, x, z, params) + (cb-1)*log(L) + (db-1)*log(1-L) + lp;
                if (g <= logy) break;
                L = L - w;
                J = J - 1;
            }
            
            while (K>0)
            { if (R>=upper) break;
                xx2[fd] = R;
                g = loglikxSP(xx2, x, z, params) + (cb-1)*log(R) + (db-1)*log(1-R) + lp;
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
            xx2[fd] = x2;
            gx2 = loglikxSP(xx2, x, z, params) + (cb-1)*log(x2) + (db-1)*log(1-x2) + lp;
            
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

