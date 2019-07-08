#ifndef GIBBS_WEB_S
#define GIBBS_WEB_S

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "web.h"


template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs_web_scale(Web<cluster_type,param_type,data_type>& W, IntegerMatrix Index, NumericVector logA, NumericVector logB, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    int n = W.data.size(); // number of data points
    vector<int> I; // indices of items in restricted Gibbs
    IntegerMatrix z(n_samples,n); // record of assignments
    int npairs = Index.nrow(); // number of pairs of data points

    for (int r=0; r<n_samples; r++) {
        //double rate=0.0;
        vector<int> c1O;
        vector<int> c2O;
        for (int s=0; s<spacing; s++) {
            
            // randomly choose a pair of anchor elements (husbands)
            int ij = ceil(unif_rand()*npairs)-1;
            
            int i1 = Index(ij,0)-1;
            int i2 = Index(ij,1)-1;

            // get the corresponding clusters of the husbands
            int id1 = W.get_cluster_id(i1);
            int id2 = W.get_cluster_id(i2);
            if (id1==id2) { id2 = W.create_cluster(); } // if i1 and i2 belong to the same cluster, make an empty cluster for c2.
            cluster_type* c1 = W.clusters.item(id1);
            cluster_type* c2 = W.clusters.item(id2);
            c1O = c1->items;
            c2O = c2->items;

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
            
            // (c1->size() == c1O.size()) && (c2->size() == c2O.size())
            //if(equal(c1O.begin(), c1O.end(), c1->items.begin()) &&
             //  equal(c2O.begin(), c2O.end(), c2->items.begin())){
             //   rate+=0.0;
            //}else{
              //  rate+=1.0/spacing;
           // }
            
            //cout << "size c1=" << c1->size() << " "<< c1O.size() << endl;
            //cout << "size c2=" << c2->size() << " "<< c2O.size() << endl;
            
            // if either cluster is empty, remove it
            if (c1->size()==0) W.remove_cluster(id1);
            if (c2->size()==0) W.remove_cluster(id2);
            
            //equal(c1O->items.begin(), c1O->items.end(), c1->items.begin()) &&
            //equal(c2O->items.begin(), c2O->items.end(), c2->items.begin())
            
        }
        // record assignments
        for (int i=0; i<n; i++) z(r,i) = W.get_cluster_id(i);
       // cout << "rate=" << rate << endl;

    }
    return z;
}

// More efficient than giving all blocked pairs in advance
template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs_web_scaleBL(Web<cluster_type,param_type,data_type>& W, List zblock, IntegerVector Idblock, NumericVector logA, NumericVector logB, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    int n = W.data.size(); // number of data points
    vector<int> I; // indices of items in restricted Gibbs
    IntegerMatrix z(n_samples,n); // record of assignments
    int nb = max(Idblock);
    //cout << "nb=" << nb << endl;
    
    for (int r=0; r<n_samples; r++) {
        //double rate=0.0;
        vector<int> c1O;
        vector<int> c2O;
        for (int s=0; s<spacing; s++) {
            
            // randomly choose a pair of anchor elements (husbands)
            
            /* this samples a block first */
            int bl = ceil(unif_rand()*nb)-1;
            //cout << "bl=" << bl << endl;
            IntegerVector zbl = zblock[bl];
            int nbl = zbl.size();
            //cout << "nbl=" << nbl << endl;
            
            
            while(nbl==1){
                bl = ceil(unif_rand()*nb)-1;
                zbl = zblock[bl];
                nbl = zbl.size();
            }

            
            int idr1 = ceil(unif_rand()*nbl)-1;
            //cout << "idr1=" << idr1 << endl;
            int i1 = zbl[idr1]-1;
            //cout << "i1=" << i1 << endl;
            int idr2 = ceil(unif_rand()*(nbl-1))-1; idr2 += (idr2>=idr1);
            //cout << "idr2=" << idr2 << endl;
            int i2 = zbl[idr2]-1;
            //cout << "i2=" << i2 << endl;

            
            /*this samples a record first
            int i1 = ceil(unif_rand()*n)-1;
            int bl = Idblock[i1];
            IntegerVector zbl = zblock[bl-1];
            int nbl = zbl.size();
            
            int idr = ceil(unif_rand()*nbl)-1;
            int i2 = zbl[idr]-1;
            */
            /*
            int nbl = 1;
            IntegerVector zbl;
            int i1;
            int bl;
            
            do{
              i1 = ceil(unif_rand()*n)-1;
              bl = Idblock[i1];
              zbl = zblock[bl-1];
              nbl = zbl.size();
            }while(nbl==1);
            
            int idr = ceil(unif_rand()*nbl)-1;
            int i2 = zbl[idr]-1;
            */
            
             //while(i1==i2){
               // idr = ceil(unif_rand()*nbl)-1;
               // i2 = zbl[idr]-1;
             //}
            
            //cout << "bl=" << bl << " "<< "nbl=" << nbl << " "<< "i1=" << i1 << " " <<"i2=" << i2 << endl;
            
            // get the corresponding clusters of the husbands
            int id1 = W.get_cluster_id(i1);
            int id2 = W.get_cluster_id(i2);
            if (id1==id2) { id2 = W.create_cluster(); } // if i1 and i2 belong to the same cluster, make an empty cluster for c2.
            cluster_type* c1 = W.clusters.item(id1);
            cluster_type* c2 = W.clusters.item(id2);
            c1O = c1->items;
            c2O = c2->items;
            
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
            
            // (c1->size() == c1O.size()) && (c2->size() == c2O.size())
            //if(equal(c1O.begin(), c1O.end(), c1->items.begin()) &&
            //  equal(c2O.begin(), c2O.end(), c2->items.begin())){
            //   rate+=0.0;
            //}else{
            //  rate+=1.0/spacing;
            // }
            
            //cout << "size c1=" << c1->size() << " "<< c1O.size() << endl;
            //cout << "size c2=" << c2->size() << " "<< c2O.size() << endl;
            
            // if either cluster is empty, remove it
            if (c1->size()==0) W.remove_cluster(id1);
            if (c2->size()==0) W.remove_cluster(id2);
            
            //equal(c1O->items.begin(), c1O->items.end(), c1->items.begin()) &&
            //equal(c2O->items.begin(), c2O->items.end(), c2->items.begin())
            
        }
        // record assignments
        for (int i=0; i<n; i++) z(r,i) = W.get_cluster_id(i);
        // cout << "rate=" << rate << endl;
        
    }
    return z;
}


template<class cluster_type, class param_type, class data_type>
IntegerMatrix gibbs_web_scaleMix(Web<cluster_type,param_type,data_type>& W, IntegerMatrix index, NumericVector logA, NumericVector logB, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    int n = W.data.size(); // number of data points
    vector<int> I; // indices of items in restricted Gibbs
    IntegerMatrix z(n_samples,n); // record of assignments
    // int npairs = n;
  //if((Rf_isMatrix(Index)==TRUE)){
        int npairs = index.nrow(); // number of pairs of data points
  // }
    
    for (int r=0; r<n_samples; r++) {
        for (int s=0; s<spacing; s++) {
            
            int ij, i1, i2;
            // randomly choose a pair of anchor elements (husbands)
            // s % 3 implies 2 iterations from best pairs and 1 from random pairs
            if((s % 3)!=0){
                ij = ceil(unif_rand()*npairs)-1;
                i1 = index(ij,0)-1;
                i2 = index(ij,1)-1;
            }else{
                i1 = ceil(unif_rand()*n)-1;
                i2 = ceil(unif_rand()*(n-1))-1; i2 += (i2>=i1);
            }
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
IntegerMatrix gibbs_web_scaleBlock(Web<cluster_type,param_type,data_type>& W, IntegerMatrix index, NumericVector logA, NumericVector logB, IntegerVector block, bool mix, int n_samples, int spacing) {
    RNGScope scope; // initialize random number generator
    int n = W.data.size(); // number of data points
    vector<int> I; // indices of items in restricted Gibbs
    IntegerMatrix z(n_samples,n+1); // record of assignments
    // int npairs = n;
    //if((Rf_isMatrix(Index)==TRUE)){
    int npairs = index.nrow(); // number of pairs of data points
    // }
    int ij, i1, i2;
    int nb = block.size(); //number of blocking fields
    int neq;
    
    for (int r=0; r<n_samples; r++) {
        int rate=0;
        vector<int> c1O;
        vector<int> c2O;
        for (int s=0; s<spacing; s++) {
            // randomly choose a pair of anchor elements (husbands)
            // s % 3 implies 2 iterations from best pairs and 1 from random pairs
            if(mix==TRUE){
                if(((s % 3)!=0)){
                    ij = ceil(unif_rand()*npairs)-1;
                    i1 = index(ij,0)-1;
                    i2 = index(ij,1)-1;
                }else{
                    do {
                        i1 = ceil(unif_rand()*n)-1;
                        i2 = ceil(unif_rand()*(n-1))-1; i2 += (i2>=i1);
                        neq = 0;
                        for(int b=0; b<nb; b++) {
                            neq += (W.data[i1][block[b]-1]==W.data[i2][block[b]-1]);
                        }
                        //cout << "neq=" << neq << endl;
                    }while(neq < nb);
                }
            }else{
                do {
                    i1 = ceil(unif_rand()*n)-1;
                    i2 = ceil(unif_rand()*(n-1))-1; i2 += (i2>=i1);
                    neq = 0;
                    for(int b=0; b<nb; b++) {
                        neq += (W.data[i1][block[b]-1]==W.data[i2][block[b]-1]);
                    }
                    //cout << "neq=" << neq << endl;
                }while(neq < nb);
            }
            // get the corresponding clusters of the chaperones
            int id1 = W.get_cluster_id(i1);
            int id2 = W.get_cluster_id(i2);
            if (id1==id2) { id2 = W.create_cluster(); } // if i1 and i2 belong to the same cluster, make an empty cluster for c2.
            cluster_type* c1 = W.clusters.item(id1);
            cluster_type* c2 = W.clusters.item(id2);
            c1O = c1->items;
            c2O = c2->items;

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
            
            /*
            // Checking whether the clusters I started with changed at all
            // This isn't a good criterion, better check if the chaperones end up together
            sort(c1O.begin(), c1O.end());
            sort(c2O.begin(), c2O.end());
            sort(c1->items.begin(),c1->items.end());
            sort(c2->items.begin(),c2->items.end());
            
            //equal(c1O.begin(), c1O.end(), c1->items.begin()) &&
            //equal(c2O.begin(), c2O.end(), c2->items.begin())
            if((c1O != c1->items) || (c2O != c2->items))
            {
                rate+=1;
            }else{
                rate+=0;
            }
            */
            //cout << "size c1=" << c1->size() << " "<< c1O.size() << endl;
            //cout << "size c2=" << c2->size() << " "<< c2O.size() << endl;
            
            // Checking whether the chaperones ended up in the same cluster
            int id1n = W.get_cluster_id(i1);
            int id2n = W.get_cluster_id(i2);
            if((id1!=id2) && (id1n==id2n))
            {
                rate+=1;
            }else{
                rate+=0;
            }
            
            
            // if either cluster is empty, remove it
            if (c1->size()==0) W.remove_cluster(id1);
            if (c2->size()==0) W.remove_cluster(id2);
            
            //equal(c1O->items.begin(), c1O->items.end(), c1->items.begin()) &&
            //equal(c2O->items.begin(), c2O->items.end(), c2->items.begin())
            
                    }
        // record assignments
        for (int i=0; i<n; i++) z(r,i) = W.get_cluster_id(i);
        cout << "rate=" << rate << endl;
        z(r,n) = rate;
    }
    return z;
}




#endif


