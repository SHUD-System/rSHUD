#include <Rcpp.h>
#include <iostream>
#include <math.h>
using namespace Rcpp;

//' Calculate the Polygon Area;
//' @param X required
//' @param Y required
//' @return Area
//' @export
// [[Rcpp::export]]
double polygonArea(NumericVector X, NumericVector Y)  {
  int n=X.size();
  double area = 0.;
  
  // Calculate value of shoelace formula
  int j = n - 1;
  for (int i = 0; i < n; i++)
  {
    area += (X(j) + X(i) ) * (Y(j) - Y(i));
    j = i;  // j is previous vertex to i
  }
  if (area > 0){
    return 0.5 * area;
  }else{
    return -0.5 * area;
  }
}
