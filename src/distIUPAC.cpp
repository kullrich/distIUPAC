#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

//' @export distIUPAC
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List distIUPAC( Rcpp::StringVector myvector ) {
  std::map<std::string, double> iupac_dist;
  iupac_dist["AA"]=0.0;iupac_dist["AC"]=1.0;iupac_dist["AG"]=1.0;iupac_dist["AT"]=1.0;iupac_dist["AR"]=0.5;iupac_dist["AY"]=1.0;iupac_dist["AS"]=1.0;iupac_dist["AW"]=0.5;iupac_dist["AK"]=1.0;iupac_dist["AM"]=0.5;iupac_dist["AB"]=-1.0;iupac_dist["AD"]=-1.0;iupac_dist["AH"]=-1.0;iupac_dist["AV"]=-1.0;iupac_dist["A."]=-1.0;iupac_dist["A-"]=-1.0;iupac_dist["AN"]=-1.0;iupac_dist["AX"]=-1.0;
  iupac_dist["CA"]=1.0;iupac_dist["CC"]=0.0;iupac_dist["CG"]=1.0;iupac_dist["CT"]=1.0;iupac_dist["CR"]=1.0;iupac_dist["CY"]=0.5;iupac_dist["CS"]=0.5;iupac_dist["CW"]=1.0;iupac_dist["CK"]=1.0;iupac_dist["CM"]=0.5;iupac_dist["CB"]=-1.0;iupac_dist["CD"]=-1.0;iupac_dist["CH"]=-1.0;iupac_dist["CV"]=-1.0;iupac_dist["C."]=-1.0;iupac_dist["C-"]=-1.0;iupac_dist["CN"]=-1.0;iupac_dist["CX"]=-1.0;
  iupac_dist["GA"]=1.0;iupac_dist["GC"]=1.0;iupac_dist["GG"]=0.0;iupac_dist["GT"]=1.0;iupac_dist["GR"]=0.5;iupac_dist["GY"]=1.0;iupac_dist["GS"]=0.5;iupac_dist["GW"]=1.0;iupac_dist["GK"]=0.5;iupac_dist["GM"]=1.0;iupac_dist["GB"]=-1.0;iupac_dist["GD"]=-1.0;iupac_dist["GH"]=-1.0;iupac_dist["GV"]=-1.0;iupac_dist["G."]=-1.0;iupac_dist["G-"]=-1.0;iupac_dist["GN"]=-1.0;iupac_dist["GX"]=-1.0;
  iupac_dist["TA"]=1.0;iupac_dist["TC"]=1.0;iupac_dist["TG"]=1.0;iupac_dist["TT"]=0.0;iupac_dist["TR"]=1.0;iupac_dist["TY"]=0.5;iupac_dist["TS"]=1.0;iupac_dist["TW"]=0.5;iupac_dist["TK"]=0.5;iupac_dist["TM"]=1.0;iupac_dist["TB"]=-1.0;iupac_dist["TD"]=-1.0;iupac_dist["TH"]=-1.0;iupac_dist["TV"]=-1.0;iupac_dist["T."]=-1.0;iupac_dist["T-"]=-1.0;iupac_dist["TN"]=-1.0;iupac_dist["TX"]=-1.0;
  iupac_dist["RA"]=0.5;iupac_dist["RC"]=1.0;iupac_dist["RG"]=0.5;iupac_dist["RT"]=1.0;iupac_dist["RR"]=0.0;iupac_dist["RY"]=1.0;iupac_dist["RS"]=1.0;iupac_dist["RW"]=1.0;iupac_dist["RK"]=1.0;iupac_dist["RM"]=1.0;iupac_dist["RB"]=-1.0;iupac_dist["RD"]=-1.0;iupac_dist["RH"]=-1.0;iupac_dist["RV"]=-1.0;iupac_dist["R."]=-1.0;iupac_dist["R-"]=-1.0;iupac_dist["RN"]=-1.0;iupac_dist["RX"]=-1.0;
  iupac_dist["YA"]=1.0;iupac_dist["YC"]=0.5;iupac_dist["YG"]=1.0;iupac_dist["YT"]=0.5;iupac_dist["YR"]=1.0;iupac_dist["YY"]=0.0;iupac_dist["YS"]=1.0;iupac_dist["YW"]=1.0;iupac_dist["YK"]=1.0;iupac_dist["YM"]=1.0;iupac_dist["YB"]=-1.0;iupac_dist["YD"]=-1.0;iupac_dist["YH"]=-1.0;iupac_dist["YV"]=-1.0;iupac_dist["Y."]=-1.0;iupac_dist["Y-"]=-1.0;iupac_dist["YN"]=-1.0;iupac_dist["YX"]=-1.0;
  iupac_dist["SA"]=1.0;iupac_dist["SC"]=0.5;iupac_dist["SG"]=0.5;iupac_dist["ST"]=1.0;iupac_dist["SR"]=1.0;iupac_dist["SY"]=1.0;iupac_dist["SS"]=0.0;iupac_dist["SW"]=1.0;iupac_dist["SK"]=1.0;iupac_dist["SM"]=1.0;iupac_dist["SB"]=-1.0;iupac_dist["SD"]=-1.0;iupac_dist["SH"]=-1.0;iupac_dist["SV"]=-1.0;iupac_dist["S."]=-1.0;iupac_dist["S-"]=-1.0;iupac_dist["SN"]=-1.0;iupac_dist["SX"]=-1.0;
  iupac_dist["WA"]=0.5;iupac_dist["WC"]=1.0;iupac_dist["WG"]=1.0;iupac_dist["WT"]=0.5;iupac_dist["WR"]=1.0;iupac_dist["WY"]=1.0;iupac_dist["WS"]=1.0;iupac_dist["WW"]=0.0;iupac_dist["WK"]=1.0;iupac_dist["WM"]=1.0;iupac_dist["WB"]=-1.0;iupac_dist["WD"]=-1.0;iupac_dist["WH"]=-1.0;iupac_dist["WV"]=-1.0;iupac_dist["W."]=-1.0;iupac_dist["W-"]=-1.0;iupac_dist["WN"]=-1.0;iupac_dist["WX"]=-1.0;
  iupac_dist["KA"]=1.0;iupac_dist["KC"]=1.0;iupac_dist["KG"]=0.5;iupac_dist["KT"]=0.5;iupac_dist["KR"]=1.0;iupac_dist["KY"]=1.0;iupac_dist["KS"]=1.0;iupac_dist["KW"]=1.0;iupac_dist["KK"]=0.0;iupac_dist["KM"]=1.0;iupac_dist["KB"]=-1.0;iupac_dist["KD"]=-1.0;iupac_dist["KH"]=-1.0;iupac_dist["KV"]=-1.0;iupac_dist["K."]=-1.0;iupac_dist["K-"]=-1.0;iupac_dist["KN"]=-1.0;iupac_dist["KX"]=-1.0;
  iupac_dist["MA"]=0.5;iupac_dist["MC"]=0.5;iupac_dist["MG"]=1.0;iupac_dist["MT"]=1.0;iupac_dist["MR"]=1.0;iupac_dist["MY"]=1.0;iupac_dist["MS"]=1.0;iupac_dist["MW"]=1.0;iupac_dist["MK"]=1.0;iupac_dist["MM"]=0.0;iupac_dist["MB"]=-1.0;iupac_dist["MD"]=-1.0;iupac_dist["MH"]=-1.0;iupac_dist["MV"]=-1.0;iupac_dist["M."]=-1.0;iupac_dist["M-"]=-1.0;iupac_dist["MN"]=-1.0;iupac_dist["MX"]=-1.0;
  iupac_dist["BA"]=-1.0;iupac_dist["BC"]=-1.0;iupac_dist["BG"]=-1.0;iupac_dist["BT"]=-1.0;iupac_dist["BR"]=-1.0;iupac_dist["BY"]=-1.0;iupac_dist["BS"]=-1.0;iupac_dist["BW"]=-1.0;iupac_dist["BK"]=-1.0;iupac_dist["BM"]=-1.0;iupac_dist["BB"]=-1.0;iupac_dist["BD"]=-1.0;iupac_dist["BH"]=-1.0;iupac_dist["BV"]=-1.0;iupac_dist["B."]=-1.0;iupac_dist["B-"]=-1.0;iupac_dist["BN"]=-1.0;iupac_dist["BX"]=-1.0;
  iupac_dist["DA"]=-1.0;iupac_dist["DC"]=-1.0;iupac_dist["DG"]=-1.0;iupac_dist["DT"]=-1.0;iupac_dist["DR"]=-1.0;iupac_dist["DY"]=-1.0;iupac_dist["DS"]=-1.0;iupac_dist["DW"]=-1.0;iupac_dist["DK"]=-1.0;iupac_dist["DM"]=-1.0;iupac_dist["DB"]=-1.0;iupac_dist["DD"]=-1.0;iupac_dist["DH"]=-1.0;iupac_dist["DV"]=-1.0;iupac_dist["D."]=-1.0;iupac_dist["D-"]=-1.0;iupac_dist["DN"]=-1.0;iupac_dist["DX"]=-1.0;
  iupac_dist["HA"]=-1.0;iupac_dist["HC"]=-1.0;iupac_dist["HG"]=-1.0;iupac_dist["HT"]=-1.0;iupac_dist["HR"]=-1.0;iupac_dist["HY"]=-1.0;iupac_dist["HS"]=-1.0;iupac_dist["HW"]=-1.0;iupac_dist["HK"]=-1.0;iupac_dist["HM"]=-1.0;iupac_dist["HB"]=-1.0;iupac_dist["HD"]=-1.0;iupac_dist["HH"]=-1.0;iupac_dist["HV"]=-1.0;iupac_dist["H."]=-1.0;iupac_dist["H-"]=-1.0;iupac_dist["HN"]=-1.0;iupac_dist["HX"]=-1.0;
  iupac_dist["VA"]=-1.0;iupac_dist["VC"]=-1.0;iupac_dist["VG"]=-1.0;iupac_dist["VT"]=-1.0;iupac_dist["VR"]=-1.0;iupac_dist["VY"]=-1.0;iupac_dist["VS"]=-1.0;iupac_dist["VW"]=-1.0;iupac_dist["VK"]=-1.0;iupac_dist["VM"]=-1.0;iupac_dist["VB"]=-1.0;iupac_dist["VD"]=-1.0;iupac_dist["VH"]=-1.0;iupac_dist["VV"]=-1.0;iupac_dist["V."]=-1.0;iupac_dist["V-"]=-1.0;iupac_dist["VN"]=-1.0;iupac_dist["VX"]=-1.0;
  iupac_dist[".A"]=-1.0;iupac_dist[".C"]=-1.0;iupac_dist[".G"]=-1.0;iupac_dist[".T"]=-1.0;iupac_dist[".R"]=-1.0;iupac_dist[".Y"]=-1.0;iupac_dist[".S"]=-1.0;iupac_dist[".W"]=-1.0;iupac_dist[".K"]=-1.0;iupac_dist[".M"]=-1.0;iupac_dist[".B"]=-1.0;iupac_dist[".D"]=-1.0;iupac_dist[".H"]=-1.0;iupac_dist[".V"]=-1.0;iupac_dist[".."]=-1.0;iupac_dist[".-"]=-1.0;iupac_dist[".N"]=-1.0;iupac_dist[".X"]=-1.0;
  iupac_dist["-A"]=-1.0;iupac_dist["-C"]=-1.0;iupac_dist["-G"]=-1.0;iupac_dist["-T"]=-1.0;iupac_dist["-R"]=-1.0;iupac_dist["-Y"]=-1.0;iupac_dist["-S"]=-1.0;iupac_dist["-W"]=-1.0;iupac_dist["-K"]=-1.0;iupac_dist["-M"]=-1.0;iupac_dist["-B"]=-1.0;iupac_dist["-D"]=-1.0;iupac_dist["-H"]=-1.0;iupac_dist["-V"]=-1.0;iupac_dist["-."]=-1.0;iupac_dist["--"]=-1.0;iupac_dist["-N"]=-1.0;iupac_dist["-X"]=-1.0;
  iupac_dist["NA"]=-1.0;iupac_dist["NC"]=-1.0;iupac_dist["NG"]=-1.0;iupac_dist["NT"]=-1.0;iupac_dist["NR"]=-1.0;iupac_dist["NY"]=-1.0;iupac_dist["NS"]=-1.0;iupac_dist["NW"]=-1.0;iupac_dist["NK"]=-1.0;iupac_dist["NM"]=-1.0;iupac_dist["NB"]=-1.0;iupac_dist["ND"]=-1.0;iupac_dist["NH"]=-1.0;iupac_dist["NV"]=-1.0;iupac_dist["N."]=-1.0;iupac_dist["N-"]=-1.0;iupac_dist["NN"]=-1.0;iupac_dist["NX"]=-1.0;
  iupac_dist["XA"]=-1.0;iupac_dist["XC"]=-1.0;iupac_dist["XG"]=-1.0;iupac_dist["XT"]=-1.0;iupac_dist["XR"]=-1.0;iupac_dist["XY"]=-1.0;iupac_dist["XS"]=-1.0;iupac_dist["XW"]=-1.0;iupac_dist["XK"]=-1.0;iupac_dist["XM"]=-1.0;iupac_dist["XB"]=-1.0;iupac_dist["XD"]=-1.0;iupac_dist["XH"]=-1.0;iupac_dist["XV"]=-1.0;iupac_dist["X."]=-1.0;iupac_dist["X-"]=-1.0;iupac_dist["XN"]=-1.0;iupac_dist["XX"]=-1.0;
  int n = myvector.size();
  Rcpp::NumericMatrix distMatrix(n, n);
  CharacterVector myvectornames = myvector.attr("names");
  colnames(distMatrix) = myvectornames;
  rownames(distMatrix) = myvectornames;
  Rcpp::NumericMatrix sitesMatrix(n, n);
  colnames(sitesMatrix) = myvectornames;
  rownames(sitesMatrix) = myvectornames;
  int nsites = myvector[1].size();
  for( int i=0; i < n; i++ ){
    sitesMatrix(i,i) = nsites;
    for( int j=i; j < n; j++ ){
      sitesMatrix(j,j) = nsites;
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
      }
      distMatrix(i,j) = eqnum / ij_n;
      distMatrix(j,i) = eqnum / ij_n;
      sitesMatrix(i,j) = ij_n;
      sitesMatrix(j,i) = ij_n;
    }
  }
  return Rcpp::List::create(Rcpp::Named("distIUPAC") = distMatrix, Rcpp::Named("sitesUsed") = sitesMatrix);
}
