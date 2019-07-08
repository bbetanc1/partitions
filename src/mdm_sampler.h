#ifndef MDM_SAMPLER_H
#define MDM_SAMPLER_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "clustering.h"
#include "gibbs.h"

// Multiple independent Dirichlet-Multinomials
class MDM_Params {
    public:
        vector< vector<double> > a;
};
class MDM_Cluster {
        int n; // number of elements in the cluster
        int L; // number of fields
        double ll; // log-likelihood
        vector< vector<int> > c; // counts for each field
        vector< vector<double> > a; // params for each field
        vector<double> a0; // sum of params, for each field
        bool initialized;
        
    public:
        // Must use a default constructor such as this (taking no arguments), in order for sets of clusters to work properly
        MDM_Cluster() { initialized = false; } // DO NOT MODIFY TO TAKE ARGUMENTS
        int size() { return n; }

        void reset() {
            ll = 0;
        }
        void initialize(MDM_Params p) {
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
                initialized = true;
            }
        }
        void insert(vector<int> x) {
            //if (!initialized) stop("(Cluster error) Attempting to insert item into uninitialized cluster.");
            // update statistics
            for (int i=0; i<L; i++) {
                ll += log(a[i][x[i]] + c[i][x[i]]) - log(a0[i] + n);
                c[i][x[i]] += 1;
            }
            n += 1;
        }
        void remove(vector<int> x) {
            //if (!initialized) stop("(Cluster error) Attempting to remove item from uninitialized cluster.");
            n -= 1;
            for (int i=0; i<L; i++) {
                c[i][x[i]] -= 1;
                ll -= log(a[i][x[i]] + c[i][x[i]]) - log(a0[i] + n);
            }
        }
        // log (marginal) likelihood
        double log_lik() { return ll; }

        // log (marginal) likelihood if x were added to the cluster
        double log_lik_w(vector<int> x) {
            double llw = ll;
            for (int i=0; i<L; i++) llw += log(a[i][x[i]] + c[i][x[i]]) - log(a0[i] + n);
            return llw;
        }
};
void validate(vector< vector<int> > data, MDM_Params params) {
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
IntegerMatrix MDM_Sampler(IntegerMatrix data_, NumericVector A, NumericVector B, List params, int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    // convert parameter list to MDM_Params
    MDM_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    validate(data,p);
    // create clustering object
    Clustering< MDM_Cluster,MDM_Params,vector<int> > C(data,p,true);
    // run sampler and return the results
    IntegerMatrix z = gibbs(C,A,B,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


/*** R
n <- 10  # number of data points
gamma <- 1  # PPP parameter
a <- list(rep(1,3),rep(.5,5),rep(2,2))  # Dirichlet parameters
rcat <- function(n,p) sample.int(length(p),n,TRUE,p)  # sample from categorical distn
data <- cbind(rcat(n,a[[1]]), rcat(n,a[[2]]), rcat(n,a[[3]]))  # sample some test data
MDM_Sampler(data,rep(1,n),rep(gamma,n),a,20)  # run Gibbs sampler for PPP model
*/
// For a DPM, the last line would instead be:
// MDM_Sampler(data,seq(1,n),rep(alpha,n),a,20)



#endif


