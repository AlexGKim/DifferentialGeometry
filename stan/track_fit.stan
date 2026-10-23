// One template track, fit in magnitude space by the truncated model: the
// mean functions b, d and the placement, with the departure modes off.
// This is how the orders of Section 2.3 are measured --- what it takes to
// reproduce a track to 0.02 mag is a question about the light curve, not
// about the log rate or the tangent angle separately --- and how the
// SALT2-placed start of Section 5 is refined.  Run with optimize.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  real z;
  int<lower=0> N;
  vector[N] dt;
  vector[N] f;
  vector<lower=0>[N] sig_f;
  array[N] int<lower=1, upper=2> band;
  // A ridge of the form of Eq. (coefprior), lambda k^{2 alpha} c_k^2: modes
  // the track does not constrain (arclength no epoch reaches) go to zero
  // instead of taking whatever the unconstrained fit leaves there.
  real<lower=0> ridge;
  real<lower=0> ridge_alpha;
}
transformed data {
#include dgsn_grid_tdata.stan
  vector[Kq - 1] a = rep_vector(0, Kq - 1);
  vector[Kpsi - 1] c = rep_vector(0, Kpsi - 1);
  vector[Kq - 1] ridge_q;
  vector[Kpsi - 1] ridge_psi;
  for (k in 1:(Kq - 1)) ridge_q[k] = ridge * pow(k, 2 * ridge_alpha);
  for (k in 1:(Kpsi - 1)) ridge_psi[k] = ridge * pow(k, 2 * ridge_alpha);
}
parameters {
  real m_peak;
  real colr;
  real x_s;
  real delta_t;
  vector[Kq - 1] b;
  vector[Kpsi - 1] d;
}
model {
  matrix[N, 4] lf = log_flux(dt, band, z, m_peak, colr, x_s, 0, 0, delta_t,
                             a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                             Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                             edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u,
                             xg, wg, P0, T);
  target += -0.5 * neg2_loglike(f, sig_f, exp(lf[:, 1]));
  target += -0.5 * (dot_product(ridge_q, square(b)) + dot_product(ridge_psi, square(d)));
}
