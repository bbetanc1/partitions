#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "set.h"

template<class cluster_type, class param_type, class data_type>
class Clustering {
        static const int unassigned = -1; // special value indicating unnassigned elements
        param_type params;

    public:
        vector<int> assignments; // cluster id associated with each data point
        vector<data_type> data;
        Set<cluster_type> clusters; // set of clusters
        cluster_type empty; // an empty cluster
        
        // Create new clustering object
        Clustering(vector<data_type> data_, param_type params_, bool init_separate) {
            data = data_; // vector copy (automatically allocates memory and copies entries into it)
            params = params_; // struct copy (automatically allocates memory and copies entries into it)
            int n = data.size();
            if (n==0) stop("(Clustering error) Data vector must have at least one entry.");
            if (init_separate) {
                // put each element in a cluster by itself
                for (int i=0; i<n; i++) {
                    int id = clusters.insert(); // make a new cluster
                    clusters.item(id)->initialize(params); // initialize it
                    clusters.item(id)->insert(data[i]);
                    assignments.push_back(id);
                }
            } else {
                int id = clusters.insert(); // make a new cluster
                clusters.item(id)->initialize(params); // initialize it
                // put all elements in it
                for (int i=0; i<n; i++) clusters.item(id)->insert(data[i]);
                assignments.resize(n,id);
            }
            // initialize empty cluster
            empty.initialize(params);
        }
        
        // Remove element i from its cluster.
        void remove(int i) {
            cluster_type *c = clusters.item(assignments[i]);
            c->remove(data[i]);
            if (c->size()==0) {
                clusters.remove(assignments[i]);
                c->reset();
            }
            assignments[i] = unassigned;
        }
        
        // Add element i to the cluster with given id.
        void insert(int i, int id) {
            assignments[i] = id;
            clusters.item(id)->insert(data[i]);
        }
        
        // Add element i to a new cluster.
        int insert_new(int i) {
            int id = clusters.insert(); // make a new cluster
            clusters.item(id)->initialize(params); // initialize it
            insert(i,id);
            return id;
        }
};



#endif













