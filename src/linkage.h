#ifndef LINKAGE_H
#define LINKAGE_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "set.h"
#define NIL -1
#define UNASSIGNED -86

template<class cluster_type, class param_type, class data_type>
class Linkage {
        param_type params;
        vector<int> assignments; // cluster id associated with each data point
        vector<int> dbindex; // index of database containing each datapoint

    public:
        vector<data_type> data;
        Set<cluster_type> clusters; // set of clusters
        int ndb; // number of databases

        // Create new Linkage object
        Linkage(vector<data_type> data_, vector<int> dbindex_, param_type params_, int ndb_) {
            data = data_; // vector/struct copy (automatically allocates memory and copies entries into it)
            dbindex = dbindex_;
            params = params_;
            ndb = ndb_;
            int n = data.size();
            if (n==0) stop("(Linkage error) Data vector must have at least one entry.");
            assignments.resize(n,UNASSIGNED);
            // put each element in a cluster by itself
            for (int i=0; i<n; i++) {
                int id = create_cluster();
                insert(i,id);
            }
        }

        // Create an empty cluster.
        int create_cluster() {
            int id = clusters.insert(); // make a new cluster
            clusters.item(id)->initialize(params,ndb); // initialize it
            return id;
        }

        // Remove cluster with given id.
        void remove_cluster(int id) {
            cluster_type *c = clusters.item(id);
            if (c->size()!=0) stop("(Linkage error) Attempting to remove a non-empty cluster.");
            c->reset();
            clusters.remove(id);
        }

        // Return id of cluster containing element i.
        int get_cluster_id(int i) {
            return assignments[i];
        }

        // Remove element i from its cluster.
        void remove(int i) {
            clusters.item(assignments[i])->remove(data[i],i,dbindex[i]);
            assignments[i] = UNASSIGNED;
        }

        // Add element i to the cluster with given id.
        void insert(int i, int id) {
            clusters.item(id)->insert(data[i],i,dbindex[i]);
            assignments[i] = id;
        }

        // Swap elements i1 and i2 between their respective clusters.
        void swap(int i1, int i2) {
            int id1 = assignments[i1];
            int id2 = assignments[i2];
            cluster_type *c1 = clusters.item(id1);
            cluster_type *c2 = clusters.item(id2);
            c1->remove(data[i1],i1,dbindex[i1]);
            c2->remove(data[i2],i2,dbindex[i2]);
            c1->insert(data[i2],i2,dbindex[i2]);
            c2->insert(data[i1],i1,dbindex[i1]);
            assignments[i1] = id2;
            assignments[i2] = id1;
        }
};



#endif













