#ifndef LOGSUMEXP_H
#define LOGSUMEXP_H

#include <cmath>

double Log(double x) {
    if (x==0) { return -INFINITY; }
    return log(x);
}

// [[Rcpp::export]]
double logsumexp(double a, double b) {
    double m = fmax(a,b);
    if (m==-INFINITY) { return -INFINITY; }
    return Log(exp(a-m) + exp(b-m)) + m;
}

double maxi(vector<double> q)
{
    double maxq;
    maxq = q[1];
    
    for(int l = 0; l< int(q.size()) ;l++){
        if(q[l]>maxq){
            maxq = q[l];
        }
    }
    return maxq;
}

double logsumexpv(vector<double> q) {
    double m = maxi(q);
    if (m==-INFINITY) {
        return -INFINITY;
    }else{
        double sle=0.0;
        for(int l = 0; l< int(q.size()) ;l++){
            sle += exp(q[l]-m);
        }
        return Log(sle) + m;
    }
}

#endif

