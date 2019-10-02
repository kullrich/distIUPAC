#include <Rcpp.h>
#include <string.h>
#include <unordered_map>
using namespace Rcpp;

//' @useDynLib distIUPAC, .registration = TRUE
//' @import Rcpp
//' @export rcpp_iupacString2diploidString
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::StringVector rcpp_iupacString2diploidString( std::string a, int nsites, std::string name_a = "a1", std::string name_b = "a2") {
  std::unordered_map<std::string, std::string> iupac_code_map;
  iupac_code_map["A"]="AA";
  iupac_code_map["B"]="NN";
  iupac_code_map["C"]="CC";
  iupac_code_map["D"]="NN";
  iupac_code_map["G"]="GG";
  iupac_code_map["H"]="NN";
  iupac_code_map["K"]="GT";
  iupac_code_map["M"]="AC";
  iupac_code_map["N"]="NN";
  iupac_code_map["R"]="AG";
  iupac_code_map["S"]="CG";
  iupac_code_map["T"]="TT";
  iupac_code_map["V"]="NN";
  iupac_code_map["W"]="AT";
  iupac_code_map["X"]="NN";
  iupac_code_map["Y"]="CT";
  iupac_code_map["-"]="NN";
  iupac_code_map["+"]="NN";
  iupac_code_map["."]="NN";
  std::string b = a;
  for( int s=0; s < nsites; s++){
    std::string as;
    std::string ar;
    std::string br;
    as = a[s];
    ar = iupac_code_map[as];
    br = iupac_code_map[as];
    ar = ar[0];
    br = br[1];
    a.replace(s,1,ar);
    b.replace(s,1,br);
  }
  return Rcpp::StringVector::create(Rcpp::Named(name_a) = a, Rcpp::Named(name_b) = b);
}
