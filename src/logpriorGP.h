#ifndef LOGPRIORGP_H
#define LOGPRIORGP_H

#include <cmath>

double logpriorGP(double x0, int k0, int N, int M, int Khat, NumericVector x1,
                IntegerVector Ml, NumericVector hpriorpar) {
    NumericVector xx1 = x1;
    xx1[k0] = x0;
    Ml.push_front(0);
        double lp = 0.0;
        //double a = xx1[0];
        //double q = xx1[1];
        double alpha = xx1[2];
        double r = xx1[3];
        double p = xx1[4];
        double ag1 = hpriorpar[0];
        double bg1 = hpriorpar[1];
        double ag2 = hpriorpar[2];
        double bg2 = hpriorpar[3];
        double ab1 = hpriorpar[4];
        double bb1 = hpriorpar[5];
        double term=0.0;
        // I am using the joint distribution P(alpha, r, p|C)
        // I know a lot of terms are constant but I am keeping everything to make it easier to check
        // everything is there.
    
    for (int m = 1; m <= M ; m++) {
            double lp0 = lgamma(m + r) + (m)*log(p) + r*log(1-p) -
            lgamma(r) - lgamma(m + 1) - log(1 - pow(1-p, r));
            term = term + Ml[m]*lgamma(m+1) + lgamma(Ml[m] + alpha*exp(lp0)) -
            lgamma(alpha*exp(lp0));
    }
   
    //cout << "term=" << " " << term << endl;
    
    lp = lgamma(alpha) - lgamma(alpha + Khat) + term  + (ag1-1)*log(alpha) - bg1*alpha +
    (ag2-1)*log(r) - bg2*r + (ab1 - 1)*log(p) + (bb1 - 1)*log(1 - p);
    
    
   // cout << "lp=" << " " << lp << endl;
        /*
        for (int m=1; m < M + 1; m++) {
            double lp0 = lgamma(m - 1 + r) + (m-1)*log(p) + r*log(1-p) -
            lgamma(r) - lgamma(m);
            term = term + Ml[m]*lgamma(m+1) + lgamma(Ml[m] + alpha*exp(lp0)) -
            lgamma(alpha*exp(lp0));
        }

        lp = lgamma(a + Khat) + Khat*log(q) + a*log(1-q) + lgamma(alpha) - lgamma(a) -
        lgamma(alpha + Khat) + term  + (ag1-1)*log(alpha) - bg1*alpha + (ag2-1)*log(r) -
        bg2*r + (ab1 - 1)*log(p) + (bb1 - 1)*log(1 - p);
        */
    
       return lp;
}

#endif

