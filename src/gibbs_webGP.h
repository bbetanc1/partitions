#ifndef GIBBS_WEBGP_H
#define GIBBS_WEBGP_H

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "web2.h"


template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs_webGP(Web<cluster_type,param_type,data_type>& W, NumericVector alphas, IntegerVector Ml,
    NumericVector logA, NumericVector logB, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    int n = W.data.size(); // number of data points
    //int maxk = alphas.size()-1; //max cluster size
    vector<int> I; // indices of items in restricted Gibbs
    //vector<int> Ml; // number of clusters of size l
    //Ml.resize(maxk);
    IntegerMatrix z(n_samples,n); // record of assignments
    
    for (int r=0; r<n_samples; r++) { //n_samples=1
        for (int s=0; s<spacing; s++) {// spacing=20000
            // randomly choose a pair of anchor elements (husbands)
            int i1 = ceil(unif_rand()*n)-1;
            int i2 = ceil(unif_rand()*(n-1))-1; i2 += (i2>=i1);
            
            // get the corresponding clusters of the husbands
            int id1 = W.get_cluster_id(i1);
            int id2 = W.get_cluster_id(i2);
            if (id1==id2) { id2 = W.create_cluster(); } // if i1 and i2 belong to the same cluster, make an empty cluster for c2.
            cluster_type* c1 = W.clusters.item(id1);
            cluster_type* c2 = W.clusters.item(id2);
            int t = W.clusters.size()-1;//number of active clusters t=|C|
            //cout << "n_active=" << t << endl;
            if (t==logB.size()) stop("Number of clusters exceeds range of B vector.");
            
            // get the indices of all elements of husbands and wives in the two clusters, and store them all in I.
            I = c1->items; I.insert(I.end(), c2->items.begin(), c2->items.end());
            
            for (int rep=0; rep<2; rep++) {  // two restricted Gibbs sweeps (ATTN: maybe try doing more than 2?)
                for (int j=0; j<int(I.size()); j++) { // for each element in each of the two clusters
                    int i = I[j];
                    
                    int n1 = c1->size();
                    int n2 = c2->size();
                    
                    // check for forbidden moves, the continue statement ignores everything ahead within the loop
                    if ((n1==0) || (n2==0)) { // if both husbands are together
                        if ((i!=i1) && (i!=i2)) continue; // if i is a wife, then she can't move
                        // (if i is a husband, he can move)
                    } else { // otherwise, the two husbands are not together
                        if ((i==i1) && (W.cluster_of(i1)->size()>1)) continue; // if i is a husband with >0 wives, he can't move
                        if ((i==i2) && (W.cluster_of(i2)->size()>1)) continue; //  "  "  "  "
                        // (if i is a wife, she can move)
                    }
                    
                    // remove i from its current cluster
                    W.remove(i);
                    
                    //sizes of corresponding clusters of chaperones after removing i
                    int n_1 = c1->size();
                    int n_2 = c2->size();
                    
                    // compute probability of assignment for each cluster
                    // p1 is unnormalized probability of being in cluster 1, likewise for p2
                    
                    double log_p1 = c1->log_lik_wwo(W.data[i]); //log-likelihood contribution of
                    double log_p2 = c2->log_lik_wwo(W.data[i]);
                    double logB1 = 0.0;
                    
                    //cout << "log_p1=" << log_p1 << endl;
                    //cout << "log_p2=" << log_p2 << endl;
                    
                    // Note here I am doing M=N.
                    // Just for reference I am using a random start with the following number of clusters
                    // of sizes 1 to 5 for Syria2000.
                    // 1   2   3   4   5
                    //729 383 114  32   7
                    
                    // Ml is a vector of size N + 1 with the number of clusters of each size (starts at 0)
                    // alphas = alpha0*exp(lp0) where alpha0, r0, p0 are the current values of the params and
                    // lp0 = lgamma(seq(N) - 1 + r0) + (seq(N)-1)*log(p0) + r0*log(1-p0) -
                    // lgamma(r0) - lgamma(seq(N)).
                    
                    // Ml and alphas have a zero in position zero, ok because of condition n_1>0
                    // logA = log(seq(N)+1), logB=log(((seq(N) + a)/(seq(N) + alpha)) * q)
                    // logA and logB also have an extra position at 0 with value 0
                    
                    
                        if(n_1 > 0){
                        // n1 = n_1 if i wasn't in cluster 1 (c1)
                        // If i was in c1 substract 1 from clusters of size n1 and add 1 to clusters of size n_1
                            if((n1 - n_1)==1){ Ml[n1]--; Ml[n_1]++;};
                            double logA1 = logA[n_1] + log(Ml[n_1+1] + alphas[n_1+1]) -
                            log(Ml[n_1] + alphas[n_1] - 1);
                            log_p1 = log_p1 + logA1;
                        }else{
                            // if i was in cluster 1 by itself substract 1 from the singleton clusters
                            if((n1 - n_1)==1) Ml[1]--;
                            logB1 = logB[t] + log(Ml[1] + alphas[1]);
                            log_p1 = log_p1 + logB1;
                        }
                        
                        if(n_2 > 0){
                        // n2 = n_2 if i wasn't in cluster 2 (c2)
                        // If i was in c2 substract 1 from clusters of size n2 and add 1 to clusters of size n_2
                            if((n2 - n_2)==1){ Ml[n2]--; Ml[n_2]++;};
                            double logA2 = logA[n_2] + log(Ml[n_2+1] + alphas[n_2+1]) -
                            log(Ml[n_2] + alphas[n_2] - 1);
                            log_p2 = log_p2 + logA2;
                        }else{
                            // if i was in cluster 2 by itself substract 1 from the singleton clusters
                            if((n2 - n_2)==1) Ml[1]--;
                            logB1 = logB[t] + log(Ml[1] + alphas[1]);
                            log_p2 = log_p2 + logB1;
                        }
                    
                        double p1 = 1/(1+exp(log_p2 - log_p1));  // normalize p1
                    
                        // reassign
                        if (unif_rand() < p1) {
                            // insert i in cluster 1 with label id1
                            W.insert(i,id1);
                            // add 1 to clusters of size n_1 + 1 and substract 1 from sizes n_1
                            Ml[n_1 + 1]++;
                            if(n_1 > 0) Ml[n_1]--;
                        }else {
                            // insert i in cluster 2 with label id2
                            W.insert(i,id2);
                            // add 1 to clusters of size n_2 + 1 and substract 1 from sizes n_2
                            Ml[n_2 + 1]++;
                            if(n_2 > 0) Ml[n_2]--;
                        }
                } //close j
            } //close rep
            // if either cluster is empty, remove it
            if (c1->size()==0) W.remove_cluster(id1);
            if (c2->size()==0) W.remove_cluster(id2);
        }
        //cout << "Singles=" << Ml[1] << endl;
        // record assignments
       for (int i=0; i<n; i++) z(r,i) = W.get_cluster_id(i);
       cout << "Ml=" << Ml[1] << " " << Ml[2]<< "  "<<  Ml[3] << endl;
 

    }
    return z;
}



#endif


