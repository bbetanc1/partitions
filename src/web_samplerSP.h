#ifndef WEB_SAMPLER_H
#define WEB_SAMPLER_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
using namespace std;

#include "web.h"
#include "gibbs_webSP.h"
#include "logsumexp.h"
#include "gibbs_webSP_Prior.h"
#include "gibbs_webSPscale.h"


// Multiple independent Dirichlet-Multinomials
class Web_Params {
    public:
        vector< vector<double> > a;
        vector<double> betas;
};
class Web_Cluster {
        int n; // number of elements in the cluster
        int L; // number of fields
        vector< vector<int> > c; // counts for each field
        vector< vector<double> > a; // params for each field
        vector<double> a0; // sum of params, for each field
        vector<double> a1;
        vector<double> betas; // distortion probability for each field
        bool initialized;
        
    public:
        vector<int> items; // list of elements belonging to this cluster

        // Must use a default constructor such as this (taking no arguments), in order for sets of clusters to work properly
        Web_Cluster() { initialized = false; } // DO NOT MODIFY TO TAKE ARGUMENTS
        int size() { return n; }

        void reset() { }

        void initialize(Web_Params p) {
            if (!initialized) {
                a = p.a;  // prior parameters
                betas = p.betas; // distortion probabilities
                L = a.size();
                // initialize statistics
                n = 0;
                c.resize(L);
                a0.assign(L,0);
                //a1.assign(L,0);
                for (int i=0; i<L; i++) {
                    int M = a[i].size();
                    c[i].assign(M,0);
                    for (int j=0; j<M; j++) a0[i] += a[i][j];
                }
                items.resize(0); // initialize items
                initialized = true;
            }
        }
        void insert(vector<int> x, int i) {
            //if (!initialized) stop("(Cluster error) Attempting to insert item into uninitialized cluster.");
            // update statistics
            for (int j=0; j<L; j++) {
                c[j][x[j]] += 1;
            }
            items.push_back(i); // add i to items
            n += 1;
        }
        void remove(vector<int> x, int i) {
            //if (!initialized) stop("(Cluster error) Attempting to remove item from uninitialized cluster.");
            n -= 1;
            for (int j=0; j<L; j++) {
                c[j][x[j]] -= 1;
            }
            items.erase(std::remove(items.begin(), items.end(), i), items.end()); // remove i from items
            // (standard way to remove an element with a given value from a vector)
        }
        // Basic Dirichlet-Multinomial likelihood
        // log (marginal) likelihood if x were added to the cluster minus old log (marginal) likelihood w.o x
        double log_lik_wwo(vector<int> x) {
            double llwwo = 0;
            for (int j=0; j<L; j++) llwwo += log(a[j][x[j]] + c[j][x[j]]) - log(a0[j] + n);
            return llwwo;
        }
    
    // Dirichlet-Multinomial and distortion model as in SteortsHallFienberg(2016)
    // log (marginal) likelihood if x were added to the cluster minus old log (marginal) likelihood w.o x
    double log_lik_wwoSP(vector<int> x) {
        double llwwo = 0;
        vector<double> slenew;
        vector<double> sleold;
        for (int j=0; j<L; j++){
            //cout << "betas=" << betas[j];
            int M = a[j].size();
            slenew.assign(M,0);
            sleold.assign(M,0);
            for (int m=0; m<M; m++){
                    //cout << "C=" << c[j][k];
                slenew[m] = log(a[j][m]) + (c[j][m] + 1*(x[j]==m))*(log(betas[j]*a[j][m] +
                                                        (1 - betas[j])) - log(betas[j]) - log(a[j][m]));
                sleold[m] = log(a[j][m]) + c[j][m]*(log(betas[j]*a[j][m] +
                                                        (1 - betas[j])) - log(betas[j]) - log(a[j][m]));
            }
            // The last two terms cancel out in chaperones when normalizing the probabilities
            llwwo += logsumexpv(slenew) - logsumexpv(sleold) + log(betas[j]) + log(a[j][x[j]]);
        }
        return llwwo;
    }
};

void validate(vector< vector<int> > data, Web_Params params) {
    vector< vector<double> > a = params.a;
    vector<double> betas = params.betas;
    int L = a.size(); // number of fields
    for (int j=0; j<L; j++) {
        if (betas[j] < 0 || betas[j] > 1) stop("Distortion parameter %d has value %f, outside acceptable values: [0,1]",j+1,betas[j]);
    }
    for (int i=0; i<int(data.size()); i++) {
        if (int(data[i].size()) != L) stop("Data point %d has %d entries, but there are %d fields in params",i+1,data[i].size(),L);
        for (int j=0; j<L; j++) {
            if ((data[i][j]<0) || (data[i][j]>=int(a[j].size()))) 
                stop("Data element data[%d][%d] has value %d, outside acceptable values: {1,...,%d}",i+1,j+1,data[i][j]+1,a[j].size());
        }
    }
}


// [[Rcpp::export]]
IntegerMatrix Web_SamplerSP(IntegerMatrix data_, IntegerVector assignments, NumericVector A,
                           NumericVector B, NumericVector distortions, List params, int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Web_Params
    Web_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    for (int i=0; i<int(distortions.size()); i++) p.betas.push_back(distortions[i]);
    validate(data,p);
    // create Web object
    int nclus = max(assignments);
    vector<int> assignments_ = as< vector<int> >(increment(assignments,-1));
    Web< Web_Cluster,Web_Params,vector<int> > W(data,assignments_,p,nclus);
    // run sampler and return the results
    IntegerMatrix z = gibbs_web(W,logA,logB,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}

// [[Rcpp::export]]
IntegerMatrix Web_SamplerSP_Par(IntegerMatrix data_, IntegerVector assignments, int tclus, NumericVector A,
                            NumericVector B, NumericVector distortions, List params, int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Web_Params
    Web_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    for (int i=0; i<int(distortions.size()); i++) p.betas.push_back(distortions[i]);
    validate(data,p);
    // create Web object
    int nclus = max(assignments);
    vector<int> assignments_ = as< vector<int> >(increment(assignments,-1));
    Web< Web_Cluster,Web_Params,vector<int> > W(data,assignments_,p,nclus);
    // run sampler and return the results
    IntegerMatrix z = gibbs_web_Par(W, tclus,logA,logB,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


// [[Rcpp::export]]
IntegerMatrix Web_SamplerSP_Prior(IntegerMatrix data_, IntegerVector assignments, NumericVector A,
                            NumericVector B, NumericVector distortions, List params,int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Web_Params
    Web_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    for (int i=0; i<int(distortions.size()); i++) p.betas.push_back(distortions[i]);
    validate(data,p);
    // create Web object
    int nclus = max(assignments);
    vector<int> assignments_ = as< vector<int> >(increment(assignments,-1));
    Web< Web_Cluster,Web_Params,vector<int> > W(data,assignments_,p,nclus);
    // run sampler and return the results
    IntegerMatrix z = gibbs_web_prior(W,logA,logB,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


// [[Rcpp::export]]
IntegerMatrix Web_SamplerSP_Scale(IntegerMatrix data_, IntegerVector assignments, IntegerMatrix Index, NumericVector A, NumericVector B, NumericVector distortions, List params, int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Web_Params
    Web_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    for (int i=0; i<int(distortions.size()); i++) p.betas.push_back(distortions[i]);
    validate(data,p);
    // create Web object
    int nclus = max(assignments);
    vector<int> assignments_ = as< vector<int> >(increment(assignments,-1));
    Web< Web_Cluster,Web_Params,vector<int> > W(data,assignments_,p,nclus);
    // run sampler and return the results
    IntegerMatrix z = gibbs_web_scale(W,Index,logA,logB,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


// [[Rcpp::export]]
IntegerMatrix Web_SamplerSP_BL(IntegerMatrix data_, IntegerVector assignments, List zblock, IntegerVector Idblock,NumericVector A, NumericVector B, NumericVector distortions, List params, int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Web_Params
    Web_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    for (int i=0; i<int(distortions.size()); i++) p.betas.push_back(distortions[i]);
    validate(data,p);
    // create Web object
    int nclus = max(assignments);
    vector<int> assignments_ = as< vector<int> >(increment(assignments,-1));
    Web< Web_Cluster,Web_Params,vector<int> > W(data,assignments_,p,nclus);
    // run sampler and return the results
    IntegerMatrix z = gibbs_web_scaleBL(W,zblock,Idblock,logA,logB,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}

// [[Rcpp::export]]
IntegerMatrix Web_SamplerSP_ScaleMix(IntegerMatrix data_, IntegerVector assignments, IntegerMatrix index, NumericVector A, NumericVector B, NumericVector distortions, List params, int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Web_Params
    Web_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    for (int i=0; i<int(distortions.size()); i++) p.betas.push_back(distortions[i]);
    validate(data,p);
    // create Web object
    int nclus = max(assignments);
    vector<int> assignments_ = as< vector<int> >(increment(assignments,-1));
    Web< Web_Cluster,Web_Params,vector<int> > W(data,assignments_,p,nclus);
    // run sampler and return the results
    IntegerMatrix z = gibbs_web_scaleMix(W,index,logA,logB,n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


// [[Rcpp::export]]
IntegerMatrix Web_SamplerSP_ScaleBlock(IntegerMatrix data_, IntegerVector assignments, IntegerMatrix index, NumericVector A, NumericVector B, NumericVector distortions, List params, IntegerVector block, bool mix,
                                       int n_samples, int spacing) {
    // convert data from matrix (with 1-based values) to vector of vectors (with 0-based values)
    vector< vector<int> > data;
    data.resize(data_.nrow());
    for (int i=0; i<data_.nrow(); i++) {
        data[i].resize(data_.ncol());
        for (int j=0; j<data_.ncol(); j++) data[i][j] = data_(i,j)-1;
    }
    NumericVector logA = log(A); logA.push_front(0);
    NumericVector logB = log(B); logB.push_front(0);
    // convert parameter list to Web_Params
    Web_Params p;
    for (int i=0; i<int(params.size()); i++) p.a.push_back(as< vector<double> >(params[i]));
    for (int i=0; i<int(distortions.size()); i++) p.betas.push_back(distortions[i]);
    validate(data,p);
    // create Web object
    int nclus = max(assignments);
    vector<int> assignments_ = as< vector<int> >(increment(assignments,-1));
    Web< Web_Cluster,Web_Params,vector<int> > W(data,assignments_,p,nclus);
    // run sampler and return the results
    IntegerMatrix z = gibbs_web_scaleBlock(W,index,logA,logB,block, mix, n_samples,spacing);
    return increment(z,1); // switch back to 1-based values
}


/*** R
n <- 10  # number of data points
gamma <- 1  # PPP parameter
a <- list(rep(1,3),rep(.5,5),rep(2,2))  # Dirichlet parameters
rcat <- function(n,p) sample.int(length(p),n,TRUE,p)  # sample from categorical distn
data <- cbind(rcat(n,a[[1]]), rcat(n,a[[2]]), rcat(n,a[[3]]))  # sample some test data
Web_Sampler(data,rep(1,n),rep(gamma,n),a,20,4)  # run Gibbs sampler for PPP model
*/
// For a DPM, the last line would instead be:
// Web_Sampler(data,seq(1,n),rep(alpha,n),a,20,4)



#endif


