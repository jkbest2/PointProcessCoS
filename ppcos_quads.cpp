#include <TMB.hpp>

// Log-normal pdf, defined as in R
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// Objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(N_vertices);        // Vertices in mesh

  // Observed point process
  DATA_VECTOR(quadrat_count);      // Number of points in the intersection of
                                   // each dual mesh cell and quadrat
  DATA_VECTOR(quadrat_exposure);   // Area of a cell/quadrat intersection

  // Areal observations
  DATA_INTEGER(use_areal);         // Use areal (1) or not (0)?
  DATA_VECTOR(dual_exposure);      // Area of each dual mesh polygon
  DATA_SCALAR(area_est);           // Estimate of total area abundance
  // DATA_SCALAR(area_sd);            // Std. deviation of area abundance estimate
  // DATA_INTEGER(N_pmarg);           // Number of values to use in marginalizing
  //                                  // Poisson
  // DATA_VECTOR(pmarg_range);        // Vector of values to use to integrate out
                                   // Poisson total abundance

  // Pass spatial precision matrix as data in order to simplify and speed up
  // estimation
  DATA_SPARSE_MATRIX(Q);

  PARAMETER(mu);                   // Mean log-density
  PARAMETER_VECTOR(spat);          // Spatial random effect

  // Initialize negative log-likelihood vector, where nll(0) is the areal
  // observation likelihood, nll(1) is the point process likelihood, and nll(2)
  // is the (expected) log-density of the latent spatial process
  vector<Type> nll(3);
  nll.setZero();

  // Get contribution of spatial random effects. With a constant precision
  // matrix, the only non-constant contribution to the log-likelihood is the
  // quadratic form. Because we have a precision matrix, we don't even need to
  // do a solve! Note that the PARAMETER_VECTOR macro apparently constructs the
  // object as an array rather than a vector, so it is necessary to use the
  // `.matrix()` method to perform matrix-vector products. Using the `.value()`
  // method at the end converts the 1×1 matrix to a scalar value.
  nll(0) = (spat.matrix().transpose() * Q * spat.matrix()).value();

  // Likelihood of partially observed (within quadrats) point process
  vector<Type> quad_lambda(N_vertices);
  for (int i = 0; i < N_vertices; i++) {
    // Calculate the expected Poisson rate within each quadrat/dual mesh
    // intersected polygon, then add the effect to the likelihood of it was
    // observed (intersected area greater than zero).
    quad_lambda(i) = quadrat_exposure(i) * exp(mu + spat(i));
    if (quadrat_exposure(i) > 0) {
      nll(1) -= dpois(quadrat_count(i), quad_lambda(i), true);
    }
  }

  // Calculate the total intensity for each dual mesh polygon.
  vector<Type> lambda(N_vertices);
  for (int i = 0; i < N_vertices; i++) {
    lambda(i) = dual_exposure(i) * exp(mu + spat(i));
  }
  Type tot_intens;
  // Use a simple zeroth-order numerical integration to get the likelihood of
  // the total abundance estimate given the current intensity map.
  tot_intens = lambda.sum();
  if (use_areal) {
    nll(2) -= dpois(area_est, tot_intens);
  }

  // Below uses a log-normal error distribution of total abundance estimate,
  // manually marginalizing out Poisson distribution over possible actual
  // abundances.
  // for (int j = 0; i < N_pmarg; i++) {
  //   nll(2) -= dpois(pmarg_range(j), tot_intens, true) +
  //     dlnorm(area_est, log(pmarg_range(j)), area_sd, true);
  // }

  REPORT(mu);
  REPORT(spat);
  REPORT(lambda);
  REPORT(nll);

  ADREPORT(tot_intens);

  return nll.sum();
}
