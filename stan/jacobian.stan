// The predicted log flux of every epoch, as a function of the raw coefficients
// the sampler moves on.  Run with fixed_param, once per perturbed coefficient,
// it gives the Jacobian of the forward map -- and from it the Fisher
// information, which is what preconditions the hierarchy.
//
// The Hessian cannot: measured at the SALT2 start, at the posterior mean and
// at both optimize modes it is indefinite, the posterior being nowhere
// locally Gaussian.  J' Sigma^-1 J is positive semi-definite by construction,
// so with the prior precision added it is always a usable metric.
//
// The gauge is applied here rather than in Python so that it is the same code
// the fit uses: a preconditioner built from a different gauge map would be
// silently wrong.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  int<lower=1> N_sn;
  int<lower=0> N;
  vector[N] dt;
  array[N] int<lower=1, upper=2> band;
  array[N_sn] int<lower=1> start;
  array[N_sn] int<lower=1> stop;
  vector[N_sn] z;
  // the raw coefficients, as hier.stan's parameters block holds them
  vector[Kq - 1] z_a;
  vector[Kq - 1] b_til;
  vector[Kpsi - 1] z_c;
  vector[Kpsi - 1] d_til;
  real<lower=0> alpha_q_m;
  real<lower=0> alpha_psi_m;
  vector[N_sn] m_peak;
  vector[N_sn] colr;
  vector[N_sn] x_s;
  vector[N_sn] x_p;
  vector[N_sn] x_psi;
  vector[N_sn] delta_t;
}
transformed data {
#include dgsn_grid_tdata.stan
  vector[Kq - 1] kq;
  vector[Kpsi - 1] kpsi;
  for (k in 1:(Kq - 1)) kq[k] = k;
  for (k in 1:(Kpsi - 1)) kpsi[k] = k;
}
generated quantities {
  vector[N] log_flux_out;
  {
    // the gauge of hier.stan, Eq. (gauge)
    vector[Kq - 1] at = pow(kq, -(alpha_q_m + 0.5)) .* z_a;
    vector[Kpsi - 1] ct = pow(kpsi, -(alpha_psi_m + 0.5)) .* z_c;
    vector[Kq - 1] a = at / sqrt(dot_self(at));
    vector[Kq - 1] b = b_til - dot_product(b_til, a) * a;
    vector[Kpsi - 1] c = ct / sqrt(quad_form(Gs, ct));
    vector[Kpsi - 1] d = d_til - (d_til' * Gs * c) * c;
    for (n in 1:N_sn) {
      int lo = start[n];
      int hi = stop[n];
      matrix[hi - lo + 1, 4] lf = log_flux(dt[lo:hi], band[lo:hi], z[n],
          m_peak[n], colr[n], x_s[n], x_p[n], x_psi[n], delta_t[n],
          a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
          Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
          edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
      log_flux_out[lo:hi] = lf[:, 1];
    }
  }
}
