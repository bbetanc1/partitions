#ifndef GIBBS_H
#define GIBBS_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "clustering.h"
#include "draw.h"

template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs(Clustering<cluster_type,param_type,data_type>& C, NumericVector A, NumericVector B, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    // GetRNGstate();
    int n = C.data.size(); // number of data points
    IntegerMatrix z(n_samples,n); // record of assignments
    Set<cluster_type>* cs = &(C.clusters); // pointer to set of clusters
    Draw D; // used for drawing Multinomial(1,p) samples

    for (int r=0; r<n_samples; r++) {
        for (int s=0; s<spacing; s++) {
            for (int i=0; i<n; i++) {
                C.remove(i); // remove i from its current cluster
                
                // compute reassignment probabilities
                D.reset();
                for (int id=cs->first(); id!=cs->flag; id=cs->next(id)) {
                    cluster_type* c = cs->item(id);
                    double log_p = log(A[c->size()-1]) + c->log_lik_w(C.data[i]) - c->log_lik();
                    D.insert(log_p,id);
                }
                double log_p = log(B[cs->size()-1]) + C.empty.log_lik_w(C.data[i]);
                D.insert(log_p,cs->flag);

                int id = D.sample(); // choose a new cluster for i

                // assign i to its new cluster
                if (id==cs->flag) {
                    if (cs->size()==B.size()) stop("Number of clusters exceeds range of B vector.");
                    id = C.insert_new(i);
                } else {
                    C.insert(i,id);
                }
            }
        }
        // record assignments
        for (int i=0; i<n; i++) z(r,i) = C.assignments[i];
    }
    // PutRNGstate();
    return z;
}



#endif


