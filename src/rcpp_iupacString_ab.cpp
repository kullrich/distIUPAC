#include <Rcpp.h>
#include <string.h>
#include <unordered_map>
using namespace Rcpp;

//' @useDynLib distIUPAC, .registration = TRUE
//' @import Rcpp
//' @export rcpp_iupacString_ab
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::StringVector rcpp_iupacString_ab( std::string a, std::string b, int nsites, std::string name = "iupacString") {
  std::unordered_map<std::string, std::string> iupac_code_map;
  iupac_code_map["AA"]="A";iupac_code_map["AC"]="M";iupac_code_map["AG"]="R";iupac_code_map["AT"]="W";iupac_code_map["AR"]="R";iupac_code_map["AY"]="H";iupac_code_map["AS"]="V";iupac_code_map["AW"]="W";iupac_code_map["AK"]="D";iupac_code_map["AM"]="M";iupac_code_map["AB"]="N";iupac_code_map["AD"]="D";iupac_code_map["AH"]="H";iupac_code_map["AV"]="V";iupac_code_map["A."]="N";iupac_code_map["A-"]="N";iupac_code_map["AN"]="N";iupac_code_map["AX"]="N";
  iupac_code_map["CA"]="M";iupac_code_map["CC"]="C";iupac_code_map["CG"]="S";iupac_code_map["CT"]="Y";iupac_code_map["CR"]="V";iupac_code_map["CY"]="Y";iupac_code_map["CS"]="S";iupac_code_map["CW"]="H";iupac_code_map["CK"]="B";iupac_code_map["CM"]="M";iupac_code_map["CB"]="B";iupac_code_map["CD"]="N";iupac_code_map["CH"]="H";iupac_code_map["CV"]="V";iupac_code_map["C."]="N";iupac_code_map["C-"]="N";iupac_code_map["CN"]="N";iupac_code_map["CX"]="N";
  iupac_code_map["GA"]="R";iupac_code_map["GC"]="S";iupac_code_map["GG"]="G";iupac_code_map["GT"]="K";iupac_code_map["GR"]="R";iupac_code_map["GY"]="B";iupac_code_map["GS"]="S";iupac_code_map["GW"]="D";iupac_code_map["GK"]="K";iupac_code_map["GM"]="V";iupac_code_map["GB"]="B";iupac_code_map["GD"]="D";iupac_code_map["GH"]="N";iupac_code_map["GV"]="V";iupac_code_map["G."]="N";iupac_code_map["G-"]="N";iupac_code_map["GN"]="N";iupac_code_map["GX"]="N";
  iupac_code_map["TA"]="W";iupac_code_map["TC"]="Y";iupac_code_map["TG"]="K";iupac_code_map["TT"]="T";iupac_code_map["TR"]="D";iupac_code_map["TY"]="Y";iupac_code_map["TS"]="B";iupac_code_map["TW"]="W";iupac_code_map["TK"]="K";iupac_code_map["TM"]="H";iupac_code_map["TB"]="B";iupac_code_map["TD"]="D";iupac_code_map["TH"]="H";iupac_code_map["TV"]="N";iupac_code_map["T."]="N";iupac_code_map["T-"]="N";iupac_code_map["TN"]="N";iupac_code_map["TX"]="N";
  iupac_code_map["RA"]="R";iupac_code_map["RC"]="V";iupac_code_map["RG"]="R";iupac_code_map["RT"]="D";iupac_code_map["RR"]="R";iupac_code_map["RY"]="N";iupac_code_map["RS"]="V";iupac_code_map["RW"]="D";iupac_code_map["RK"]="D";iupac_code_map["RM"]="H";iupac_code_map["RB"]="N";iupac_code_map["RD"]="D";iupac_code_map["RH"]="N";iupac_code_map["RV"]="V";iupac_code_map["R."]="N";iupac_code_map["R-"]="N";iupac_code_map["RN"]="N";iupac_code_map["RX"]="N";
  iupac_code_map["YA"]="H";iupac_code_map["YC"]="Y";iupac_code_map["YG"]="B";iupac_code_map["YT"]="Y";iupac_code_map["YR"]="N";iupac_code_map["YY"]="Y";iupac_code_map["YS"]="B";iupac_code_map["YW"]="H";iupac_code_map["YK"]="B";iupac_code_map["YM"]="H";iupac_code_map["YB"]="B";iupac_code_map["YD"]="N";iupac_code_map["YH"]="H";iupac_code_map["YV"]="N";iupac_code_map["Y."]="N";iupac_code_map["Y-"]="N";iupac_code_map["YN"]="N";iupac_code_map["YX"]="N";
  iupac_code_map["SA"]="V";iupac_code_map["SC"]="S";iupac_code_map["SG"]="S";iupac_code_map["ST"]="B";iupac_code_map["SR"]="V";iupac_code_map["SY"]="B";iupac_code_map["SS"]="S";iupac_code_map["SW"]="N";iupac_code_map["SK"]="B";iupac_code_map["SM"]="V";iupac_code_map["SB"]="B";iupac_code_map["SD"]="N";iupac_code_map["SH"]="N";iupac_code_map["SV"]="V";iupac_code_map["S."]="N";iupac_code_map["S-"]="N";iupac_code_map["SN"]="N";iupac_code_map["SX"]="N";
  iupac_code_map["WA"]="W";iupac_code_map["WC"]="H";iupac_code_map["WG"]="D";iupac_code_map["WT"]="W";iupac_code_map["WR"]="D";iupac_code_map["WY"]="H";iupac_code_map["WS"]="N";iupac_code_map["WW"]="W";iupac_code_map["WK"]="D";iupac_code_map["WM"]="H";iupac_code_map["WB"]="N";iupac_code_map["WD"]="D";iupac_code_map["WH"]="H";iupac_code_map["WV"]="N";iupac_code_map["W."]="N";iupac_code_map["W-"]="N";iupac_code_map["WN"]="N";iupac_code_map["WX"]="N";
  iupac_code_map["KA"]="D";iupac_code_map["KC"]="B";iupac_code_map["KG"]="K";iupac_code_map["KT"]="K";iupac_code_map["KR"]="D";iupac_code_map["KY"]="B";iupac_code_map["KS"]="B";iupac_code_map["KW"]="D";iupac_code_map["KK"]="K";iupac_code_map["KM"]="N";iupac_code_map["KB"]="B";iupac_code_map["KD"]="D";iupac_code_map["KH"]="N";iupac_code_map["KV"]="N";iupac_code_map["K."]="N";iupac_code_map["K-"]="N";iupac_code_map["KN"]="N";iupac_code_map["KX"]="N";
  iupac_code_map["MA"]="M";iupac_code_map["MC"]="M";iupac_code_map["MG"]="V";iupac_code_map["MT"]="H";iupac_code_map["MR"]="V";iupac_code_map["MY"]="H";iupac_code_map["MS"]="V";iupac_code_map["MW"]="H";iupac_code_map["MK"]="N";iupac_code_map["MM"]="M";iupac_code_map["MB"]="N";iupac_code_map["MD"]="N";iupac_code_map["MH"]="H";iupac_code_map["MV"]="V";iupac_code_map["M."]="N";iupac_code_map["M-"]="N";iupac_code_map["MN"]="N";iupac_code_map["MX"]="N";
  iupac_code_map["BA"]="N";iupac_code_map["BC"]="B";iupac_code_map["BG"]="B";iupac_code_map["BT"]="B";iupac_code_map["BR"]="N";iupac_code_map["BY"]="B";iupac_code_map["BS"]="B";iupac_code_map["BW"]="N";iupac_code_map["BK"]="B";iupac_code_map["BM"]="N";iupac_code_map["BB"]="B";iupac_code_map["BD"]="N";iupac_code_map["BH"]="N";iupac_code_map["BV"]="N";iupac_code_map["B."]="N";iupac_code_map["B-"]="N";iupac_code_map["BN"]="N";iupac_code_map["BX"]="N";
  iupac_code_map["DA"]="D";iupac_code_map["DC"]="N";iupac_code_map["DG"]="D";iupac_code_map["DT"]="D";iupac_code_map["DR"]="D";iupac_code_map["DY"]="N";iupac_code_map["DS"]="N";iupac_code_map["DW"]="D";iupac_code_map["DK"]="D";iupac_code_map["DM"]="N";iupac_code_map["DB"]="N";iupac_code_map["DD"]="D";iupac_code_map["DH"]="N";iupac_code_map["DV"]="N";iupac_code_map["D."]="N";iupac_code_map["D-"]="N";iupac_code_map["DN"]="N";iupac_code_map["DX"]="N";
  iupac_code_map["HA"]="H";iupac_code_map["HC"]="H";iupac_code_map["HG"]="N";iupac_code_map["HT"]="H";iupac_code_map["HR"]="N";iupac_code_map["HY"]="H";iupac_code_map["HS"]="N";iupac_code_map["HW"]="H";iupac_code_map["HK"]="N";iupac_code_map["HM"]="H";iupac_code_map["HB"]="N";iupac_code_map["HD"]="N";iupac_code_map["HH"]="H";iupac_code_map["HV"]="N";iupac_code_map["H."]="N";iupac_code_map["H-"]="N";iupac_code_map["HN"]="N";iupac_code_map["HX"]="N";
  iupac_code_map["VA"]="V";iupac_code_map["VC"]="V";iupac_code_map["VG"]="V";iupac_code_map["VT"]="N";iupac_code_map["VR"]="V";iupac_code_map["VY"]="N";iupac_code_map["VS"]="N";iupac_code_map["VW"]="N";iupac_code_map["VK"]="N";iupac_code_map["VM"]="V";iupac_code_map["VB"]="N";iupac_code_map["VD"]="N";iupac_code_map["VH"]="N";iupac_code_map["VV"]="V";iupac_code_map["V."]="N";iupac_code_map["V-"]="N";iupac_code_map["VN"]="N";iupac_code_map["VX"]="N";
  iupac_code_map[".A"]="N";iupac_code_map[".C"]="N";iupac_code_map[".G"]="N";iupac_code_map[".T"]="N";iupac_code_map[".R"]="N";iupac_code_map[".Y"]="N";iupac_code_map[".S"]="N";iupac_code_map[".W"]="N";iupac_code_map[".K"]="N";iupac_code_map[".M"]="N";iupac_code_map[".B"]="N";iupac_code_map[".D"]="N";iupac_code_map[".H"]="N";iupac_code_map[".V"]="N";iupac_code_map[".."]="N";iupac_code_map[".-"]="N";iupac_code_map[".N"]="N";iupac_code_map[".X"]="N";
  iupac_code_map["-A"]="N";iupac_code_map["-C"]="N";iupac_code_map["-G"]="N";iupac_code_map["-T"]="N";iupac_code_map["-R"]="N";iupac_code_map["-Y"]="N";iupac_code_map["-S"]="N";iupac_code_map["-W"]="N";iupac_code_map["-K"]="N";iupac_code_map["-M"]="N";iupac_code_map["-B"]="N";iupac_code_map["-D"]="N";iupac_code_map["-H"]="N";iupac_code_map["-V"]="N";iupac_code_map["-."]="N";iupac_code_map["--"]="N";iupac_code_map["-N"]="N";iupac_code_map["-X"]="N";
  iupac_code_map["NA"]="N";iupac_code_map["NC"]="N";iupac_code_map["NG"]="N";iupac_code_map["NT"]="N";iupac_code_map["NR"]="N";iupac_code_map["NY"]="N";iupac_code_map["NS"]="N";iupac_code_map["NW"]="N";iupac_code_map["NK"]="N";iupac_code_map["NM"]="N";iupac_code_map["NB"]="N";iupac_code_map["ND"]="N";iupac_code_map["NH"]="N";iupac_code_map["NV"]="N";iupac_code_map["N."]="N";iupac_code_map["N-"]="N";iupac_code_map["NN"]="N";iupac_code_map["NX"]="N";
  iupac_code_map["XA"]="N";iupac_code_map["XC"]="N";iupac_code_map["XG"]="N";iupac_code_map["XT"]="N";iupac_code_map["XR"]="N";iupac_code_map["XY"]="N";iupac_code_map["XS"]="N";iupac_code_map["XW"]="N";iupac_code_map["XK"]="N";iupac_code_map["XM"]="N";iupac_code_map["XB"]="N";iupac_code_map["XD"]="N";iupac_code_map["XH"]="N";iupac_code_map["XV"]="N";iupac_code_map["X."]="N";iupac_code_map["X-"]="N";iupac_code_map["XN"]="N";iupac_code_map["XX"]="N";
  for( int s=0; s < nsites; s++){
    std::string as;
    std::string bs;
    as = a[s];
    bs = b[s];
    if(as!=bs){
      a.replace(s,1,iupac_code_map[as+bs]);
    }
  }
  return Rcpp::StringVector::create(Rcpp::Named(name) = a);
}
