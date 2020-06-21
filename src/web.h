#ifndef WEB_H
#define WEB_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "set2.h"
#define UNASSIGNED -86


template<class cluster_type, class param_type, class data_type>
class Web {
        param_type params;
        vector<int> assignments; // cluster id associated with each data point

    public:
        vector<data_type> data;
        Set<cluster_type> clusters; // set of clusters
        int nclus; // number of clusters in assignments

        // Create new Web object
        Web(vector<data_type> data_, vector<int> assignments_, param_type params_, int nclus_) {
            data = data_; // copy data into Web object
            params = params_; // copy params into Web object
            assignments = assignments_; // copy initial assignment into Web object
            nclus = nclus_; // copy initialnumber of clusters into Web object
            // Note this is vector/struct copy, which automatically allocates memory and copies entries into it.
            int n = data.size();
            if (n==0) stop("(Web error) Data vector must have at least one entry.");
            
            // put each element in their respective cluster
            for (int i=0; i<n; i++) {
                int id = create_cluster1(i);
                //cout<< id << " " << assignments[i]<<endl;
                if(id < nclus){
                    for (int j=0; j<n; j++) {
                        if(assignments[j] == id){
                            insert(j, id);
                            //cout<< id << " " << assignments[j]<<endl;
                        }
                    }
                } else {
                    remove_cluster(id);
                }
            }
        }
        
        // Return id of cluster containing element i.
        int get_cluster_id(int i) {
            return assignments[i];
        }

        // Return cluster associated with a given element.
        cluster_type* cluster_of(int i) {
            return clusters.item(assignments[i]);
        }

        // Create an empty cluster.
        int create_cluster() {
            int id = clusters.insert(); // make a new cluster (this insert is defined in set2.h)
            clusters.item(id)->initialize(params); // initialize it
            // (the arrow -> dereferences pointer and accesses the function simultaneously)
            return id;
        }
    
        // Create an empty cluster.
        int create_cluster1(int i) {
            int id = clusters.insert1(i); // make a new cluster (this insert is defined in set2.h)
            clusters.item(id)->initialize(params); // initialize it
            // (the arrow -> dereferences pointer and accesses the function simultaneously)
            return id;
        }
    

        // Remove cluster with given id.
        void remove_cluster(int id) {
            cluster_type *c = clusters.item(id); // * means c is a pointer to cluster_type object
            if (c->size()!=0) stop("(Web error) Attempting to remove a non-empty cluster.");
            c->reset();
            clusters.remove(id);
        }

        // Remove element i from its cluster.
        void remove(int i) {
            clusters.item(assignments[i])->remove(data[i],i); // (this remove is in web_sampler.h)
            assignments[i] = UNASSIGNED;
        }

        // Add element i to the cluster with given id.
        void insert(int i, int id) {
            clusters.item(id)->insert(data[i],i); // (this insert is in web_sampler.h)
            assignments[i] = id;
        }
};



#endif













