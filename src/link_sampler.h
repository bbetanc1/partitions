#ifndef LINK_SAMPLER_H
#define LINK_SAMPLER_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "linkage.h"
#include "gibbs_link.h"


// Multiple independent Dirichlet-Multinomials
class Link_Params {
    public:
        vector< vector<double> > a;
};
class Link_Cluster {
        int n; // number of elements in the cluster
        int L; // number of fields
        double ll; // log-likelihood
        vector< vector<int> > c; // counts for each field
        vector< vector<double> > a; // params for each field
        vector<double> a0; // sum of params, for each field
        bool initialized;
        
    public:
        vector<int> slots; // indices of datapoints in the cluster, for each DB

        // Must use a default constructor such as this (taking no arguments), in order for sets of clusters to work properly
        Link_Cluster() { initialized = false; } // DO NOT MODIFY TO TAKE ARGUMENTS
        int size() { return n; }

        void reset() {
            ll = 0;
        }
        void initialize(Link_Params p, int ndb) {
            if (!initialized) {
                a = p.a;  // prior parameters
                L = a.size();
                // initialize statistics
                n = 0;
                ll = 0;
                c.resize(L);
                a0.assign(L,0);
                for (int i=0; i<L; i++) {
                    int M = a[i].size();
                    c[i].assign(M,0);
                    for (int j=0; j<M; j++) a0[i] += a[i][j];
                }
                // initialize slots
                slots.assign(ndb,NIL);
                initialized = true;
            }
        }
        void insert(vector<int> x, int i, int db) {
            //if (!initialized) stop("(Cluster error) Attempting to insert item into uninitialized cluster.");
            // update statistics
            for (int j=0; j<L; j++) {
                ll += log(a[j][x[j]] + c[j][x[j]]) - log(a0[j] + n);
                c[j][x[j]] += 1;
            }
            slots[db] = i;
            n += 1;
        }
        void remove(vector<int> x, int i, int db) {
            //if (!initialized) stop("(Cluster error) Attempting to remove item from uninitialized cluster.");
            n -= 1;
            for (int j=0; j<L; j++) {
                c[j][x[j]] -= 1;
                ll -= log(a[j][x[j]] + c[j][x[j]]) - log(a0[j] + n);
            }
            slots[db] = NIL;
        }
        // log (marginal) likelihood
        double log_lik() { return ll; }

        // log (marginal) likelihood if x were added to the cluster
        double log_lik_w(vector<int> x) {
            double llw = ll;
            for (int j=0; j<L; j++) llw += log(a[j][x[j]] + c[j][x[j]]) - log(a0[j] + n);
            return llw;
        }

        // log (marginal) likelihood if x were added to the cluster and y were removed
        double log_lik_wwo(vector<int> x, vector<int> y) {
            double llwwo = ll;
            for (int j=0; j<L; j++) llwwo += log(a[j][x[j]] + c[j][x[j]]) - log(a[j][y[j]] + c[j][y[j]]);
            return llwwo;
        }
};
void validate(vector< vector<int> > data, Link_Params params) {
    vector< vector<double> > a = params.a;
    int L = a.size(); // number of fields
    for (int i=0; i<int(data.size()); i++) {
        if (int(data[i].size()) != L) stop("Data point %d has %d entries, but there are %d fields in params",i+1,data[i].size(),L);
        for (int j=0; j<L; j++) {
            if ((data[i][j]<0) || (data[i][j]>=int(a[j].size()))) 
                stop("Data element data[%d][%d] has value %d, outside acceptable values: {1,...,%d}",i+1,j+1,data[i][j]+1,a[j].size());
        }
    }
}


// [[Rcpp::export]]
IntegerMatrix Link_Sampler(IntegerMatrix data_, IntegerVector dbindex,
        NumericVector A, NumericVector B, List params, int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Link_Params
    Link_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    validate(data,p);
    // create linkage object
    int ndb = max(dbindex);
    vector<int> dbindex_ = as< vector<int> >(increment(dbindex,-1));
    Linkage< Link_Cluster,Link_Params,vector<int> > L(data,dbindex_,p,ndb);
    // run sampler and return the results
    IntegerMatrix z = gibbs_link(L,logA,logB,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


/*** R
n <- 10  # number of data points
dbindex <- c(rep(1,5),rep(2,5))
gamma <- 1  # PPP parameter
a <- list(rep(1,3),rep(.5,5),rep(2,2))  # Dirichlet parameters
rcat <- function(n,p) sample.int(length(p),n,TRUE,p)  # sample from categorical distn
data <- cbind(rcat(n,a[[1]]), rcat(n,a[[2]]), rcat(n,a[[3]]))  # sample some test data
Link_Sampler(data,dbindex,rep(1,n),rep(gamma,n),a,20,4)  # run Gibbs sampler for PPP model
*/
// For a DPM, the last line would instead be:
// Link_Sampler(data,dbindex,seq(1,n),rep(alpha,n),a,20,4)



#endif


