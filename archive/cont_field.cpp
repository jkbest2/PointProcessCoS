#include <TMB.hpp>

// Construct the precision matrix for a stationary GMRF constructed using the
// SPDE approximation from R-INLA. This is specific to alpha = 2, which
// corresponds to a Matern with nu = 1 in R^2
template<class Type>
Eigen::SparseMatrix<Type> Q_stat(Type tau, Type kappa2,
                                 Eigen::SparseMatrix<Type> C0,
                                 Eigen::SparseMatrix<Type> G1,
                                 Eigen::SparseMatrix<Type> G2) {
  int n = C0.cols();
  Eigen::SparseMatrix<Type> Q(n, n);
  Q = pow(tau, 2) * (pow(kappa2, 2) * C0 + Type(2.0) * kappa2 * G1 + G2);
  return(Q);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(N_vertices);        // Number of mesh vertices
  DATA_INTEGER(N_obs);             // Number of observations
  DATA_VECTOR(obs_value);          // Value of the process; locations defined by
                                   // A matrix

  // Matrices for constructing precision matrix Q
  DATA_SPARSE_MATRIX(C0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // `A` matrix for projecting mesh to observation locations
  DATA_SPARSE_MATRIX(A);

  PARAMETER(mu);                   // Mean log-density
  PARAMETER(sigma);                // Observation standard deviation
  PARAMETER_VECTOR(spat);          // Spatial random effect

  PARAMETER(log_tau);              // Precision parameter
  Type tau = exp(log_tau);
  PARAMETER(log_kappa2);           // Precision & range parameter
  Type kappa2 = exp(log_kappa2);

  vector<Type> nll(2);
  nll.setZero();

  // Get contribution of spatial random effects
  Eigen::SparseMatrix<Type> Q(N_vertices, N_vertices);
  Q = Q_stat(tau, kappa2, C0, G1, G2);
  density::GMRF_t<Type> gmrf(Q);
  nll(1) += gmrf(spat);

  // Calculate intensity function
  vector<Type> pred_value(N_obs);
  pred_value = A * spat + mu;

  // Likelihood
  nll(0) -= sum(dnorm(obs_value, pred_value, sigma, true));

  Type sigma_spat;
  sigma_spat = sqrt(1 / (4 * PI * kappa2 * tau * tau));
  Type rho_spat;
  rho_spat = sqrt(8 / kappa2);

  ADREPORT(mu);
  ADREPORT(sigma_spat);
  ADREPORT(rho_spat);

  return nll.sum();
}
