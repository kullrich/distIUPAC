#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

//' @useDynLib distIUPAC, .registration = TRUE
//' @import Rcpp
//' @export distIUPACmatrix
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List distIUPACmatrix( Rcpp::StringVector dnavector, Rcpp::NumericMatrix scoreMatrix ) {
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
  int n = dnavector.size();
  Rcpp::NumericMatrix distMatrix(n, n);
  CharacterVector dnavectornames = dnavector.attr("names");
  colnames(distMatrix) = dnavectornames;
  rownames(distMatrix) = dnavectornames;
  Rcpp::NumericMatrix sitesMatrix(n, n);
  colnames(sitesMatrix) = dnavectornames;
  rownames(sitesMatrix) = dnavectornames;
  int nsites = dnavector[1].size();
  for( int i=0; i < n; i++ ){
    for( int j=i; j < n; j++ ){
      double eqnum = 0;
      int ij_n = nsites;
      for( int s=0; s < nsites; s++){
        std::string is;
        std::string js;
        is = dnavector[i][s];
        js = dnavector[j][s];
        double ij_dist;
        ij_dist = iupac_dist[is+js];
        if(ij_dist >= 0.0){
          eqnum = eqnum + ij_dist;
        } else {
          ij_n = ij_n -1;
        };
      }
      distMatrix(i,j) = eqnum / ij_n;
      distMatrix(j,i) = eqnum / ij_n;
      sitesMatrix(i,j) = ij_n;
      sitesMatrix(j,i) = ij_n;
    }
  }
  return Rcpp::List::create(Rcpp::Named("distIUPAC") = distMatrix, Rcpp::Named("sitesUsed") = sitesMatrix);
}
