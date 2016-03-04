// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat matOmega1(arma::mat W, double rho) {
  // Omega1 - SAR(1)
  int n = W.n_rows;
  arma::mat I = arma::eye<arma::mat>(n, n);
  arma::mat S1 = I - rho * W;
  return arma::inv_sympd(S1.t() * S1);
}

// [[Rcpp::export]]
arma::mat matOmega2(int nTime, double rho) {
  // Omega2 - AR(1)
  arma::mat Ome2(nTime, nTime);
  Ome2.fill(0.0);
  for(int i = 0; i < nTime; ++i) {
    Ome2.diag(i) += pow(rho, i);
  }
  Ome2 += Ome2.t();
  Ome2.diag() *= 0.0;
  Ome2.diag() += 1;
  return 1/(1-pow(rho, 2)) * Ome2;
}

// [[Rcpp::export]]
arma::mat matBlockDiagonal(arma::mat X, int n) {
  arma::mat XX(X.n_cols*n, X.n_cols*n);
  XX.fill(0.0);
  for(int r = 0; r < n; ++r) {
    XX.submat(arma::span(r * X.n_cols, (r + 1) * X.n_cols - 1), arma::span(r * X.n_cols, (r + 1) * X.n_cols - 1)) =
      X;
  }
  return XX;
}

// [[Rcpp::export]]
arma::mat matA(double sigma2, arma::mat Ome2, int nDomains, arma::colvec sigmaSamplingError) {
  // Some Matrix A - used to speed up VInv
  arma::mat Ad = sigma2 * Ome2;
  arma::mat A = matBlockDiagonal(Ad, nDomains);
  A.diag() += sigmaSamplingError;
  return A;
}

// [[Rcpp::export]]
Rcpp::List matVInvST(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  // Inverse of variance covariance matrix of spatio-temporal model
  arma::mat Ome2 = matOmega2(Z1.n_rows / W.n_rows, rho2);
  arma::mat Ome1 = matOmega1(W, rho1);
  arma::mat Ad = sigma2 * Ome2;
  arma::mat A = matBlockDiagonal(Ad, W.n_rows);
  A.diag() += sigmaSamplingError;
  arma::mat Ainv = A;

  for(int r = 0; r < W.n_rows; ++r) {
    Ainv.submat(arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1), arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1)) =
      arma::inv(Ainv.submat(arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1), arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1)));
  }

  arma::mat V = sigma1 * Z1 * Ome1 * Z1.t() + A;
  arma::mat AinvZ1 = Ainv * Z1;
  arma::mat Ome1inv = arma::inv(sigma1 * Ome1);
  arma::mat Vinv = Ainv - AinvZ1 * arma::inv(Ome1inv + Z1.t() * AinvZ1) * AinvZ1.t();
  return Rcpp::List::create(Rcpp::Named("V", V),
                            Rcpp::Named("Vinv", Vinv));
}

// [[Rcpp::export]]
arma::mat matVST(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  // Variance-Covariance of spatio-temporal model
  arma::mat Ome2 = matOmega2(Z1.n_rows / W.n_rows, rho2);
  arma::mat A = matA(sigma2, Ome2, W.n_cols, sigmaSamplingError);
  arma::mat Ome1 = matOmega1(W, rho1);
  return sigma1 * Z1 * Ome1 * Z1.t() + A;
}

// [[Rcpp::export]]
arma::mat matVDerS1(arma::mat Ome1, arma::mat Z1) {
  // derivative of V with respect to sigma1
  return Z1 * Ome1 * Z1.t();
}

// [[Rcpp::export]]
arma::mat matVDerS2(arma::mat Ome2, int nDomains) {
  // derivative of V with respect to sigma2
  return matBlockDiagonal(Ome2, nDomains);
}

// [[Rcpp::export]]
arma::mat matVDerR1(double rho1, double sigma1, arma::mat Z1, arma::mat Ome1, arma::mat W) {
  // derivative of V with respect to rho1
  return -sigma1 * Z1 * Ome1 * (-W-W.t() + 2 * rho1 * W.t() * W) * Ome1 * Z1.t();
}

// [[Rcpp::export]]
arma::mat matVDerR2(double rho2, double sigma2, arma::mat Ome2, int nDomains) {
  // derivative of V with respect to rho2
  arma::mat ome2derR2(Ome2.n_cols, Ome2.n_cols);
  ome2derR2.fill(0.0);
  for(int i = 1; i < Ome2.n_cols; ++i) {
    ome2derR2.diag(i) += i * pow(rho2, i-1);
  }
  ome2derR2 += ome2derR2.t();
  ome2derR2 = 1/(1-pow(rho2, 2)) * (ome2derR2 + 2 * rho2 * Ome2);
  return sigma2 * matBlockDiagonal(ome2derR2, nDomains);
}
