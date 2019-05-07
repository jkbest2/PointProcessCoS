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
  DATA_INTEGER(N_vertices);        // Vertices in mesh

  // Observed counts
  DATA_INTEGER(N_cells);          // Number of cells for point process
  DATA_VECTOR(pp_counts);          // Number of points in each cell
  DATA_SCALAR(pp_cell_area);       // Area of a cell
  // DATA_VECTOR(count_quads);

  // Matrices for constructing precision matrix Q
  DATA_SPARSE_MATRIX(C0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // `A` matrix for projecting mesh to observation locations
  DATA_SPARSE_MATRIX(A);
  // A matrix of R rows x N_quads columns; specifies rows of A to use for each
  // quadrat observation
  // DATA_IMATRIX(A_obs);
  // int N_quad_obs;                  // Number of fine cells in a quadrat
  // N_quad_obs = A_obs.rows();
  // int N_quads;                     // Number of quadrats observed
  // N_quads = A_obs.cols();

  PARAMETER(mu);                        // Mean log-density
  PARAMETER_VECTOR(spat);               // Spatial random effect

  PARAMETER(log_tau);                   // Precision parameter
  Type tau = exp(log_tau);
  PARAMETER(log_kappa2);                // Precision & range parameter
  Type kappa2 = exp(log_kappa2);

  vector<Type> nll(3);
  nll.setZero();

  // Get contribution of spatial random effects
  Eigen::SparseMatrix<Type> Q(N_vertices, N_vertices);
  Q = Q_stat(tau, kappa2, C0, G1, G2);
  density::GMRF_t<Type> gmrf(Q);
  nll(2) += gmrf(spat);

  // Calculate intensity function
  vector<Type> log_lambda_spat(N_cells);
  log_lambda_spat = A * spat;

  // Likelihood of total count
  vector<Type> lambda(N_cells);
  for(int i = 0; i < N_cells; i++) {
    lambda(i) = pp_cell_area * exp(mu + log_lambda_spat(i));
    nll(1) -= dpois(pp_counts(i), lambda(i), true);
  }

  //Likelihood of quadrat counts
  // vector<Type> quad_lambda(N_quads);
  // for (int quad = 0; quad < N_quads; quad++) {
  //   vector<Type> temp_lambda(N_quad_obs);
  //   for (int i = 0; i < N_quad_obs; i++) {
  //     temp_lambda(i) += exp(log_lambda(A_obs(i, quad)));
  //   }
  //   vector<Type> quad_lambda(N_quads);
  //   quad_lambda(quad) = temp_lambda.sum();
  //   nll(0) -= dpois(count_quads(quad), quad_lambda(quad), true);
  // }

  REPORT(lambda);
  // REPORT(quad_lambda);
  // REPORT(log_lambda);

  Type sigma;
  sigma = sqrt(1 / (4 * PI * kappa2));
  Type rho;
  rho = sqrt(8 / kappa2);
  ADREPORT(sigma);
  ADREPORT(rho);

  return nll.sum();
}
