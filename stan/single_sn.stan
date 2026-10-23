// Stage 2 of Section 5.2: one supernova, six parameters, the global
// functions supplied as data.  Also the fitter for the stage-1 dilation and
// translation tests and the floors, through optimize with prior_on = 0.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  vector[Kq - 1] a;               // gauge-fixed global coefficients
  vector[Kq - 1] b;
  vector[Kpsi - 1] c;
  vector[Kpsi - 1] d;
  real z;
  int<lower=0> N;                 // the epoch set of Section 4
  vector[N] dt;                   // mjd - t0
  vector[N] f;
  vector<lower=0>[N] sig_f;
  array[N] int<lower=1, upper=2> band;
  // Priors of Section 5, with the population hyperparameters as data.
  int<lower=0, upper=1> prior_on;
  real mu_p;
  real mu_psi;
  real mu_delta;
  real<lower=0> s_delta;
  real xs_bar;
  real m_bright;
}
transformed data {
#include dgsn_grid_tdata.stan
}
parameters {
  real m_peak;
  real colr;
  real x_s;
  real x_p;
  real x_psi;
  real delta_t;
}
model {
  matrix[N, 4] lf = log_flux(dt, band, z, m_peak, colr, x_s, x_p, x_psi, delta_t,
                             a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                             Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                             edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u,
                             xg, wg, P0, T);
  target += -0.5 * neg2_loglike(f, sig_f, exp(lf[:, 1]));
  if (prior_on) {
    x_p ~ normal(mu_p, 1);
    x_psi ~ normal(mu_psi, 1);
    delta_t ~ normal(mu_delta, s_delta);
    x_s ~ normal(xs_bar, 1);
    m_peak ~ normal(m_bright, 1);
    colr ~ normal(0, 1);
  }
}
generated quantities {
  vector[N] f_hat;
  real s_pmin;                    // the interval check of Section 3
  real s_pmax;
  real kappa0 = turning_deriv_row(0, Kpsi, ulo, uhi, T) * (d + x_psi * c) / delta;
  real frac_tapered;              // epochs with taper below 0.99
  real frac_clamped;              // epochs moved by more than eps_s of the interval
  {
    matrix[N, 4] lf = log_flux(dt, band, z, m_peak, colr, x_s, x_p, x_psi, delta_t,
                               a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                               Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                               edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u,
                               xg, wg, P0, T);
    f_hat = exp(lf[:, 1]);
    frac_tapered = 0;
    frac_clamped = 0;
    for (i in 1:N) {
      frac_tapered += (lf[i, 2] < log(0.99));
      frac_clamped += (abs(lf[i, 4] - lf[i, 3]) > eps_s * (smax - smin));
    }
    frac_tapered /= N;
    frac_clamped /= N;
  }
  {
    // The data window on the sampled t_max: s at [-15, +40] days.
    vector[2] ends = arclength([pmin + 5, pmax - 5]', x_s, b + x_p * a, Kq, vlo, vhi,
                               delta_p, m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
    s_pmin = ends[1];
    s_pmax = ends[2];
  }
}
