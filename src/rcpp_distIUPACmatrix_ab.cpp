
#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

//' @useDynLib distIUPAC, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @export rcpp_distIUPACmatrix_ab
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_distIUPACmatrix_ab( std::string a, std::string b, int nsites, Rcpp::NumericMatrix scoreMatrix ) {
  std::unordered_map<std::string, double> iupac_dist;
  int nCols = scoreMatrix.ncol();
  int nRows = scoreMatrix.nrow();
  CharacterVector Colnames = colnames(scoreMatrix);
  CharacterVector Rownames = rownames(scoreMatrix);
  for( int is = 0; is < nCols; is++ ){
    for( int js = 0; js < nRows; js++){
      std::string isName = "";
      isName = Colnames(is);
      std::string jsName = "";
      jsName = Rownames(js);
      iupac_dist[isName+jsName] = scoreMatrix(is,js);
    }
  }
  double eqnum = 0;
  int ab_n = nsites;
  for( int s=0; s < nsites; s++){
    std::string as;
    std::string bs;
    as = a[s];
    bs = b[s];
    double ab_dist;
    ab_dist = iupac_dist[as+bs];
    if(ab_dist >= 0.0){
      eqnum = eqnum + ab_dist;
    } else {
      ab_n = ab_n -1;
    };
  }
  return Rcpp::NumericVector::create(Rcpp::Named("distIUPAC") = eqnum, Rcpp::Named("sitesUsed") = ab_n);
}
