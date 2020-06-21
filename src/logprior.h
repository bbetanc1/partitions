#ifndef LOGPRIOR_H
#define LOGPRIOR_H

#include <cmath>

double logprior(double x0, int k0, int N, int Khat, NumericVector x1,
                IntegerVector Nk, NumericVector hpriorpar, std::string Prior) {
    NumericVector xx1 = x1;
    xx1[k0-1] = x0;
    double lp = 0.0;
    if (Prior=="DP"){
        double theta = xx1[0];
        double ag1 = hpriorpar[0];
        double bg1 = hpriorpar[1];
        double term=0.0;
        for (int i=0; i<Khat; i++) {
            term = term + lgamma(Nk[i]);
        }
        lp = lgamma(theta) + Khat*log(theta) - lgamma(theta + N) +
        term + (ag1-1)*log(theta) - bg1*theta;
    }
    if (Prior=="PY"){
        double theta = xx1[0];
        double delta = xx1[1];
        double ag1 = hpriorpar[0];
        double bg1 = hpriorpar[1];
        double ab1 = hpriorpar[2];
        double bb1 = hpriorpar[3];
        double term=0.0;
        IntegerVector sK=seq_len(Khat);
        for (int i=0; i<Khat; i++)  {
            term = term + log(theta + delta*(sK[i])) +
            lgamma(Nk[i] - delta) - lgamma(1 - delta);
        }
        lp = lgamma(theta + 1) - log(theta + delta*Khat) -
        lgamma(theta + N) + term +
        (ag1-1)*log(theta) - bg1*theta + (ab1 - 1)*log(delta) +
        (bb1 - 1)*log(1 - delta);
    }
    if (Prior=="ESCNB"){
        double rnb = xx1[0];
        double pnb = xx1[1];
        double ag = hpriorpar[0]; // gamma prior for r
        double bg = hpriorpar[1];
        double ab = hpriorpar[2]; // beta prior for p
        double bb = hpriorpar[3];
        double term=0.0;
        for (int i=0; i<Khat; i++) {
            term = term + lgamma(rnb + Nk[i]) - lgamma(rnb);
        }
        lp = N*log(pnb) + lgamma(Khat + 1) + Khat*(rnb*log(1-pnb) -
                                                   log(1 - pow((1-pnb),rnb))) + term +
        (ag-1)*log(rnb) - bg*rnb + (ab - 1)*log(pnb) + (bb - 1)*log(1 - pnb);
    }

       return lp;
}

double logpriorESCD(double x0, int k0, NumericVector x1, IntegerVector Lm, NumericVector mu0,
                    NumericVector hpriorpar) {
    NumericVector xx1 = x1;
    xx1[k0-1] = x0;
    double lp = 0.0;
    double alpha = xx1[0];
    double rnb = xx1[1];
    double pnb = xx1[2];
    double ag = hpriorpar[0]; // gamma prior for r
    double bg = hpriorpar[1];
    double ab = hpriorpar[2]; // beta prior for p
    double bb = hpriorpar[3];
    int M = Lm.size();
    double term=0.0;
    for (int i=0; i<M; i++) {
        term = term + Lm[i]*lgamma(i+1) + lgamma(Lm[i] + alpha*mu0[i]) - lgamma(alpha*mu0[i]);
    }
    lp = term + (ag-1)*log(rnb) - bg*rnb + (ab - 1)*log(pnb) + (bb - 1)*log(1 - pnb);
    
    return lp;
}



#endif

