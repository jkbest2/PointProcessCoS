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

  // Covariate values at each mesh vertex
  DATA_MATRIX(X)

  // Observed point process
  DATA_VECTOR(quadrat_count);      // Number of points in the intersection of
                                   // each dual mesh cell and quadrat
  DATA_VECTOR(quadrat_exposure);   // Area of a cell/quadrat intersection

  // Areal observations
  DATA_VECTOR(dual_exposure);      // Area of each dual mesh polygon
  DATA_SCALAR(area_est);           // Estimate of total area abundance

  PARAMETER_VECTOR(beta);          // Mean log-density

  // Initialize negative log-likelihood vector, where nll(0) is the areal
  // observation likelihood, nll(1) is the point process likelihood, and nll(2)
  // is the (expected) log-density of the latent spatial process
  vector<Type> nll(2);
  nll.setZero();

  // Get the mean vector using the design matrix and the betas
  vector<Type> mu(N_vertices);
  mu = X * beta;

  // Likelihood of partially observed (within quadrats) point process
  vector<Type> quad_lambda(N_vertices);
  for (int i = 0; i < N_vertices; i++) {
    // Calculate the expected Poisson rate within each quadrat/dual mesh
    // intersected polygon, then add the effect to the likelihood of it was
    // observed (intersected area greater than zero).
    quad_lambda(i) = quadrat_exposure(i) * exp(mu(i));
    if (quadrat_exposure(i) > 0) {
      nll(0) -= dpois(quadrat_count(i), quad_lambda(i), true);
    }
  }

  // Calculate the total intensity for each dual mesh polygon.
  vector<Type> lambda(N_vertices);
  for (int i = 0; i < N_vertices; i++) {
    lambda(i) = dual_exposure(i) * exp(mu(i));
  }
  Type tot_intens;
  // Use a simple zeroth-order numerical integration to get the likelihood of
  // the total abundance estimate given the current intensity map.
  tot_intens = lambda.sum();
  nll(1) -= dpois(area_est, tot_intens);

  REPORT(mu);
  REPORT(lambda);
  REPORT(nll);

  ADREPORT(tot_intens);

  return nll.sum();
}
