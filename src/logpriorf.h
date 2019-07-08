#ifndef LOGPRIORF_H
#define LOGPRIORF_H

#include <cmath>

double logpriorf(NumericVector x1, int N, int Khat,
                IntegerVector Nk, NumericVector hpriorpar, std::string Prior) {
    NumericVector xx1 = x1;
    double lp = 0.0;
    if (Prior=="DP"){ //Prior == "DP" str1.compare(str2) != 0
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
    if (Prior=="PERPS"){
        double alpha = xx1[0];
        double lambda = xx1[1];
        double ag1 = hpriorpar[0];
        double bg1 = hpriorpar[1];
        double ag2 = hpriorpar[2];
        double bg2 = hpriorpar[3];
        lp = N*log(lambda) + Khat*log(alpha) - alpha - Khat*lambda +
        alpha*exp(-lambda) +
        (ag1-1)*log(alpha) - bg1*alpha + (ag2-1)*log(lambda) - bg2*lambda;
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
        //sum(log(theta + delta*(1:Khat)) + lgamma(Nk - delta) - lgamma(1 - delta))
        lp = lgamma(theta + 1) - log(theta + delta*Khat) -
             lgamma(theta + N) + term +
             (ag1-1)*log(theta) - bg1*theta + (ab1 - 1)*log(delta) +
             (bb1 - 1)*log(1 - delta);
    }
    if (Prior=="NBNBO"){
        double anb = xx1[0];
        double qnb = xx1[1];
        double rnb = xx1[2];
        double pnb = xx1[3];
        double ag1 = hpriorpar[0];
        double bg1 = hpriorpar[1];
        double ag2 = hpriorpar[2];
        double bg2 = hpriorpar[3];
        double ab1 = hpriorpar[4];
        double bb1 = hpriorpar[5];
        double ab2 = hpriorpar[6];
        double bb2 = hpriorpar[7];
        double term=0.0;
        for (int i=0; i<Khat; i++) {
            term = term + lgamma(rnb + Nk[i]) - lgamma(rnb);
        }
        lp = N*log(pnb) + lgamma(anb + Khat) - lgamma(anb) +
             anb*log(1 - qnb) + Khat*(log(qnb) + rnb*log(1-pnb)) -
             (anb + Khat)*log(1 - qnb*pow((1-pnb),rnb)) +
             term + (ag1-1)*log(anb) - bg1*anb + (ag2-1)*log(rnb) -
             bg2*rnb + (ab1 - 1)*log(qnb) + (bb1 - 1)*log(1 - qnb) +
             (ab2 - 1)*log(pnb) + (bb2 - 1)*log(1 - pnb);
    }
    if (Prior=="NBNB"){
        double anb = xx1[0];
        double qnb = xx1[1];
        double rnb = xx1[2];
        double pnb = xx1[3];
        double ag1 = hpriorpar[0];
        double bg1 = hpriorpar[1];
        double ag2 = hpriorpar[2];
        double bg2 = hpriorpar[3];
        double ab1 = hpriorpar[4];
        double bb1 = hpriorpar[5];
        double ab2 = hpriorpar[6];
        double bb2 = hpriorpar[7];
        double term=0.0;
        for (int i=0; i<Khat; i++) {
            term = term + lgamma(rnb + Nk[i] - 1) - lgamma(rnb);
        }
        lp = lgamma(anb + Khat) - lgamma(anb) + anb*log(1 - qnb) +
        Khat*(log(qnb) + rnb*log(1-pnb)) + (N-Khat)*log(pnb) +
        term + (ag1-1)*log(anb) - bg1*anb + (ag2-1)*log(rnb) -
        bg2*rnb + (ab1 - 1)*log(qnb) +  (bb1 - 1)*log(1 - qnb) +
        (ab2 - 1)*log(pnb) + (bb2 - 1)*log(1 - pnb);
    }
    if (Prior=="NBNBF"){
        double anb = xx1[0];
        double qnb = xx1[1];
        double rnb = xx1[2];
        double pnb = xx1[3];
        double ag1 = hpriorpar[0];
        double bg1 = hpriorpar[1];
        double ag2 = hpriorpar[2];
        double bg2 = hpriorpar[3];
        double ab1 = hpriorpar[4];
        double bb1 = hpriorpar[5];
        double ab2 = hpriorpar[6];
        double bb2 = hpriorpar[7];
        double term=0.0;
        for (int i=0; i<Khat; i++) {
            term = term + lgamma(rnb + Nk[i]) - lgamma(rnb);
        }
        lp = N*log(pnb) + lgamma(anb + Khat) + Khat*(log(qnb) +
        rnb*log(1-pnb) - log(1 - pow((1-pnb),rnb))) +
        term + anb*log(1-qnb) - log(1 - pow((1-qnb),anb)) -
	lgamma(anb) + (ag1-1)*log(anb) - bg1*anb +
        (ag2-1)*log(rnb) - bg2*rnb +  (ab1 - 1)*log(qnb) +
        (bb1 - 1)*log(1 - qnb) + (ab2 - 1)*log(pnb) +
        (bb2 - 1)*log(1 - pnb);
    }
    return lp;
}

#endif

