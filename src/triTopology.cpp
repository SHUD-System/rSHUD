#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

//' Determine the topological relation among triangles. 
//' @param tri triangle defination, m x 3
//' @return topological relation of triangles. m x 4; cols = c(ID, Nabor1, Nabor2, Nabor3)
//' @export
// [[Rcpp::export]]
NumericMatrix triTopology(NumericMatrix tri) {
  long n = tri.nrow();
  int flag;
  int nid[4] = {0,1,2,0};
  long ed0[2];
  long ed1[2];
  long i, j, ii,jj;
  NumericMatrix out(n, 4);
  for(i=0; i < n; i++){
    for(j=0; j < 3; j++){
      ed0[0]=  tri(i, nid[j]);
      ed0[1] =  tri(i, nid[j+1]);
      flag = 0;
      for(ii=0; ii < n && !flag; ii++){
        if(i == ii){
          continue;
        }
        for(jj=0;jj<3 && !flag;jj++){
          ed1[0]=  tri(ii, nid[jj]);
          ed1[1] =  tri(ii, nid[jj+1]);
          if(ed0[0] ==  ed1[0]){
            if(ed0[1]== ed1[1]){
              flag = ii + 1;
            }
          }else if(ed0[1] ==  ed1[0]){
            if(ed0[0] ==  ed1[1]){
              flag = ii + 1 ;
            }
          }
        }
      }
      out(i, j+1) = flag;
    } // for of three edges.
    out(i, 0) = i;
  }// for of nrow
  return  out;
}