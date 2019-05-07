#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(N_vertices);        // Vertices in mesh
  DATA_INTEGER(N_areaobs);         // Number of areal observations

  // Observed counts
  DATA_VECTOR(pp_counts);          // Number of points in each cell
  DATA_VECTOR(pp_cell_area);       // Area of a cell

  // Areal observations
  DATA_VECTOR(areaobs);            // Areal observations

  // Pass spatial precision matrix as data in order to simplify and speed up
  // estimation
  DATA_SPARSE_MATRIX(Q);

  PARAMETER(mu);                        // Mean log-density
  PARAMETER_VECTOR(spat);               // Spatial random effect

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

  // Likelihood of point process
  vector<Type> lambda(N_vertices);
  for (int i = 0; i < N_vertices; i++) {
    // Note that `pp_cell_area == 0` outside the boundary. TMB's `dpois`
    // function doesn't handle Poisson distributions with rate 0 (a NaN is
    // returned), so here we branch on the area.
    lambda(i) = pp_cell_area(i) * exp(mu + spat(i));
    if (pp_cell_area(i) > 0) {
      nll(1) -= dpois(pp_counts(i), lambda(i), true);
    } else {
      nll(1) = 0;
    }
  }

  Type tot_intens;
  tot_intens = lambda.sum();
  for (int j = 0; j < N_areaobs; j++) {
    nll(2) -= dpois(areaobs(j), tot_intens, true);
  }

  REPORT(mu);
  REPORT(spat);
  REPORT(lambda);
  REPORT(nll);

  ADREPORT(tot_intens);

  return nll.sum();
}
