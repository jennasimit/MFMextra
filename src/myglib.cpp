#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec getEigenValues(arma::mat M) {
  arma::eig_sym(M);
}

// [[Rcpp::export]]
double sumlogeigen(arma::mat M) {
  arma::vec eM = arma::eig_sym(M);
  //eM.print();
  arma::vec lM = log(eM);
  // lM.print();
  // return(wrap(lM));
  int n = M.n_rows;
  double ret=lM(0);
  for(int i=1; i<n; i++)
    ret += lM(i);
  return ret;
  // return Rcpp::List::create(lM,ret,lM(0));
}

// [[Rcpp::export]]
arma::mat sigmult(arma::mat A, arma::mat B, arma::mat C, arma::mat D) {
  return( A - B * inv(C) * D ) ;
}

// [[Rcpp::export]]
arma::mat ApBmult(arma::mat A, double phi, double psi) {
  int n = A.n_rows;
  arma::vec v(n);
  v.fill(pow(phi,2));
  v(0,0) = pow(psi,2);
  arma::mat B = diagmat(v);
  // A.print("A");
  // B.print("B");
  return( A * B * A.t() ) ;
}

// [[Rcpp::export]]
SEXP E1helper(arma::mat pvarj, arma::mat V, arma::mat x) {
  arma::mat A = inv(V);
  arma::mat B = inv(pvarj);
  arma::mat C = inv(A+B);
  double logdetW = sumlogeigen(pvarj);
  int n = C.n_rows;
  arma::mat I = arma::eye<arma::mat>(n,n);
  arma::mat icb = I - C*B;
  arma::mat F1 = B*C*A*C*B + icb.t() * B * icb;
  arma::mat ret = -1.0*x.t() * F1 * x - logdetW - sumlogeigen(A + B);
  return(wrap(ret));
}

// [[Rcpp::export]]
double wsum(NumericVector x, NumericVector w) {
  int n=x.size();
  double ret=0.0;
  for(int i=0; i<n; i++)
    ret += x[i]*w[i];
  return ret;
}
  

// [[Rcpp::export]]
SEXP armqr(arma::mat M) {
  arma::mat Q, R;
  arma::qr(Q,R,M);
  
  return Rcpp::List::create(Rcpp::Named("Q")=Q,Rcpp::Named("R")=R);
}


// [[Rcpp::export]]
arma::mat invxxt(arma::mat x) {
  return( inv(x * x.t()) ) ;
}
// arma::mat F1(arma::mat A, arma::mat B, arma::mat C, arma::mat I) {
//   arma::mat A = 1/V;
//   arma::vec = 1/pvar;
//   arma::mat 

// }

// [[Rcpp::export]]
arma::mat cpdiff(arma::mat x, arma::vec w) {
  arma::mat x1 = cross(x,w) ;
  arma::mat x2 = pow(x,2);
  arma::mat x3 = cross(x2,w) ;
  return x1 - x3;
}
