
#include <Rcpp.h>
using namespace Rcpp;

IntegerMatrix increment(IntegerMatrix x, int a) {
    for (int i=0; i<x.nrow(); i++) for (int j=0; j<x.ncol(); j++) x(i,j) += a;
    return x;
}
IntegerVector increment(IntegerVector x, int a) {
    for (int i=0; i<x.size(); i++) x[i] += a;
    return x;
}

#include "web_samplerSP.h"
#include "loglikx.h"
#include "unislicem.h"
#include "loglikf.h"
#include "unislicef.h"

