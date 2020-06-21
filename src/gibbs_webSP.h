#ifndef GIBBS_WEB_H
#define GIBBS_WEB_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "web.h"

template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs_web(Web<cluster_type,param_type,data_type>& W, NumericVector logA, NumericVector logB, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    int n = W.data.size(); // number of data points
    vector<int> I; // indices of items in restricted Gibbs
    IntegerMatrix z(n_samples,n); // record of assignments

    for (int r=0; r<n_samples; r++) {
        for (int s=0; s<spacing; s++) {
            // randomly choose a pair of anchor elements (husbands)
            
            int i1 = ceil(unif_rand()*n)-1;
            int i2 = ceil(unif_rand()*(n-1))-1; i2 += (i2>=i1);

            // get the corresponding clusters of the husbands
            int id1 = W.get_cluster_id(i1);
            int id2 = W.get_cluster_id(i2);
            if (id1==id2) { id2 = W.create_cluster(); } // if i1 and i2 belong to the same cluster, make an empty cluster for c2.
            cluster_type* c1 = W.clusters.item(id1);
            cluster_type* c2 = W.clusters.item(id2);
            int t = W.clusters.size()-1;
            if (t==logB.size()) stop("Number of clusters exceeds range of B vector.");

            // get the indices of all elements of husbands and wives in the two clusters, and store them all in I.
            I = c1->items; I.insert(I.end(), c2->items.begin(), c2->items.end());

            for (int rep=0; rep<2; rep++) {  // two restricted Gibbs sweeps (ATTN: maybe try doing more than 2?)
                for (int j=0; j<int(I.size()); j++) { // for each element in each of the two clusters
                    int i = I[j];
                    // check for forbidden moves
                    if ((c1->size()==0) || (c2->size()==0)) { // if both husbands are together
                        if ((i!=i1) && (i!=i2)) continue; // if i is a wife, then she can't move
                        // (if i is a husband, he can move)
                    } else { // otherwise, the two husbands are not together
                        if ((i==i1) && (W.cluster_of(i1)->size()>1)) continue; // if i is a husband with >0 wives, he can't move
                        if ((i==i2) && (W.cluster_of(i2)->size()>1)) continue; //  "  "  "  "
                        // (if i is a wife, she can move)
                    }

                    // remove i from its current cluster
                    W.remove(i);

                    // compute probability of assignment for each cluster
                    // p1 is unnormalized probability of being in cluster 1, likewise for p2
                    //cout << "p_c1=" << c1->log_lik_wwoSP(W.data[i]) ;
                    double log_p1 = c1->log_lik_wwoSP(W.data[i]) + (c1->size()>0? logA[c1->size()] : logB[t]);
                    double log_p2 = c2->log_lik_wwoSP(W.data[i]) + (c2->size()>0? logA[c2->size()] : logB[t]);
                    double p1 = 1/(1+exp(log_p2 - log_p1));  // normalize p1

                    // reassign
                    if (unif_rand() < p1) { W.insert(i,id1); }
                    else { W.insert(i,id2); }
                }
            }
            // if either cluster is empty, remove it
            if (c1->size()==0) W.remove_cluster(id1);
            if (c2->size()==0) W.remove_cluster(id2);
        }
        // record assignments
        for (int i=0; i<n; i++) z(r,i) = W.get_cluster_id(i);
    }
    return z;
}


template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs_web_fbl(Web<cluster_type,param_type,data_type>& W, IntegerVector bl, NumericVector logA, NumericVector logB, int n_samples, int spacing) {
  RNGScope scope; // initialize random number generator
  int n = W.data.size(); // number of data points
  vector<int> I; // indices of items in restricted Gibbs
  IntegerMatrix z(n_samples,n); // record of assignments
  int nb = bl.size();
  
  for (int r=0; r<n_samples; r++) {
    for (int s=0; s<spacing; s++) {
      
      // randomly choose anchors in the block
      int ib1 = ceil(unif_rand()*nb)-1;
      int ib2 = ceil(unif_rand()*nb)-1;
      //ib2 += (ib2>=ib1);
      
      int i1 = bl[ib1]-1;
      int i2 = bl[ib2]-1;
      
      // get the corresponding clusters of the husbands
      int id1 = W.get_cluster_id(i1);
      int id2 = W.get_cluster_id(i2);
      if (id1==id2) { id2 = W.create_cluster(); } // if i1 and i2 belong to the same cluster, make an empty cluster for c2.
      cluster_type* c1 = W.clusters.item(id1);
      cluster_type* c2 = W.clusters.item(id2);
      int t = W.clusters.size()-1;
      if (t==logB.size()) stop("Number of clusters exceeds range of B vector.");
      
      // get the indices of all elements of husbands and wives in the two clusters, and store them all in I.
      I = c1->items; I.insert(I.end(), c2->items.begin(), c2->items.end());
      
      for (int rep=0; rep<2; rep++) {  // two restricted Gibbs sweeps (ATTN: maybe try doing more than 2?)
        for (int j=0; j<int(I.size()); j++) { // for each element in each of the two clusters
          int i = I[j];
          // check for forbidden moves
          if ((c1->size()==0) || (c2->size()==0)) { // if both husbands are together
            if ((i!=i1) && (i!=i2)) continue; // if i is a wife, then she can't move
            // (if i is a husband, he can move)
          } else { // otherwise, the two husbands are not together
            if ((i==i1) && (W.cluster_of(i1)->size()>1)) continue; // if i is a husband with >0 wives, he can't move
            if ((i==i2) && (W.cluster_of(i2)->size()>1)) continue; //  "  "  "  "
            // (if i is a wife, she can move)
          }
          
          // remove i from its current cluster
          W.remove(i);
          
          // compute probability of assignment for each cluster
          // p1 is unnormalized probability of being in cluster 1, likewise for p2
          //cout << "p_c1=" << c1->log_lik_wwoSP(W.data[i]) ;
          double log_p1 = c1->log_lik_wwoSP(W.data[i]) + (c1->size()>0? logA[c1->size()] : logB[t]);
          double log_p2 = c2->log_lik_wwoSP(W.data[i]) + (c2->size()>0? logA[c2->size()] : logB[t]);
          double p1 = 1/(1+exp(log_p2 - log_p1));  // normalize p1
          
          // reassign
          if (unif_rand() < p1) { W.insert(i,id1); }
          else { W.insert(i,id2); }
        }
      }
      // if either cluster is empty, remove it
      if (c1->size()==0) W.remove_cluster(id1);
      if (c2->size()==0) W.remove_cluster(id2);
    }
    // record assignments
    for (int i=0; i<n; i++) z(r,i) = W.get_cluster_id(i);
  }
  return z;
}



#endif


