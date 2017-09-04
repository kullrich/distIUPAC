#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

//' @export distIUPACmatrix
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::NumericMatrix distIUPACmatrix( Rcpp::StringVector myvector, Rcpp::NumericMatrix scoreMatrix ) {
  std::map<std::string, double> iupac_dist;
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
  int n = myvector.size();
  Rcpp::NumericMatrix distMatrix(n, n);
  int nsites = myvector[1].size();
  for( int i=0; i < n - 1; i++ ){
    for( int j=1; j < n; j++ ){
      double eqnum = 0;
      int ij_n = nsites;
      for( int s=0; s < nsites; s++){
        std::string is;
        std::string js;
        is = myvector[i][s];
        js = myvector[j][s];
        double ij_dist;
        ij_dist = iupac_dist[is+js];
        if(ij_dist >= 0.0){
          eqnum = eqnum + ij_dist;
        } else {
          ij_n = ij_n -1;
        };
        distMatrix(i,j) = eqnum / ij_n;
        distMatrix(j,i) = eqnum / ij_n;
      }
    }
  }
  return distMatrix;
}
