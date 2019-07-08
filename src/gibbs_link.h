#ifndef GIBBS_LINK_H
#define GIBBS_LINK_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "linkage.h"


template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs_link(Linkage<cluster_type,param_type,data_type>& L, NumericVector logA, NumericVector logB, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    int n = L.data.size(); // number of data points
    IntegerMatrix z(n_samples,n); // record of assignments

    for (int r=0; r<n_samples; r++) {
        for (int s=0; s<spacing; s++) {
            // choose a pair of records
            int i1 = ceil(unif_rand()*n)-1;
            int i2 = ceil(unif_rand()*(n-1))-1; i2 += (i2>=i1);

            // get the corresponding clusters
            int id1 = L.get_cluster_id(i1);
            int id2 = L.get_cluster_id(i2);
            cluster_type* c1 = L.clusters.item(id1);
            cluster_type* c2 = L.clusters.item(id2);
            int n1 = c1->size();
            int n2 = c2->size();

            // if c1==c2, redefine c2 to be an empty cluster.
            if (c1==c2) {
                id2 = L.create_cluster();
                c2 = L.clusters.item(id2);
                n2 = 0;
            }
            int t = L.clusters.size()-1;
            if (t==logB.size()) stop("Number of clusters exceeds range of B vector.");
            
            for (int db=0; db<L.ndb; db++) { // for database db, consider swapping the records in that db
                int j1 = c1->slots[db];
                int j2 = c2->slots[db];

                // check for forbidden swaps
                if ((j1==NIL) && (j2==NIL)) continue;
                if ((j1==NIL) && (j2==i2) && (n2>1)) continue;
                if ((j2==NIL) && (j1==i1) && (n1>1)) continue;
                
                // Gibbs swap move
                if (j1==NIL) {
                    L.remove(j2); n2 -= 1;
                    double log_p1 = c1->log_lik_w(L.data[j2]) - c1->log_lik() + (n1>0? logA[n1] : logB[t]);
                    double log_p2 = c2->log_lik_w(L.data[j2]) - c2->log_lik() + (n2>0? logA[n2] : logB[t]);
                    double p1 = 1/(1+exp(log_p2 - log_p1));
                    if (unif_rand() < p1) { L.insert(j2,id1); n1 += 1; }
                    else { L.insert(j2,id2); n2 += 1; }
                } else if (j2==NIL) {
                    L.remove(j1); n1 -= 1;
                    double log_p1 = c1->log_lik_w(L.data[j1]) - c1->log_lik() + (n1>0? logA[n1] : logB[t]);
                    double log_p2 = c2->log_lik_w(L.data[j1]) - c2->log_lik() + (n2>0? logA[n2] : logB[t]);
                    double p1 = 1/(1+exp(log_p2 - log_p1));
                    if (unif_rand() < p1) { L.insert(j1,id1); n1 += 1; }
                    else { L.insert(j1,id2); n2 += 1; }
                } else {
                    double log_lik_stay = c1->log_lik() + c2->log_lik();
                    double log_lik_swap = c1->log_lik_wwo(L.data[j2],L.data[j1]) + c2->log_lik_wwo(L.data[j1],L.data[j2]);
                    double p_swap = 1/(1+exp(log_lik_stay - log_lik_swap));
                    if (unif_rand() < p_swap) L.swap(j1,j2);
                }
            }
            // if either cluster is empty, remove it
            if (n1==0) L.remove_cluster(id1);
            if (n2==0) L.remove_cluster(id2);
        }
        // record assignments
        for (int i=0; i<n; i++) z(r,i) = L.get_cluster_id(i);
    }
    return z;
}



#endif


