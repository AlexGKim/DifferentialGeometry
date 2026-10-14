// Stages 3 and 4 of Section 5.2: the hierarchy.  Everything single_sn.stan
// held as data is sampled here: the four raw coefficient vectors, their
// scales and decays, the scatters, the population means, and six parameters
// per supernova.  Stage 4 is this program on the full training sample.
//
// The likelihood is a reduce_sum over supernovae.  Nothing per-epoch is
// saved: at stage 4 a flux per epoch per draw would be gigabytes, and a
// fixed_param pass on saved draws produces fluxes where wanted.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  int<lower=1> N_sn;
  int<lower=0> N;                  // all epochs, ragged by supernova
  vector[N] dt;
  vector[N] f;
  vector<lower=0>[N] sig_f;
  array[N] int<lower=1, upper=2> band;
  array[N_sn] int<lower=1> start;  // 1-based, inclusive
  array[N_sn] int<lower=1> stop;
  vector[N_sn] z;
  vector[N_sn] m_bright;
  real xs_bar;
  int<lower=1> grainsize;
}
transformed data {
#include dgsn_grid_tdata.stan
  vector[Kq - 1] kq;
  vector[Kpsi - 1] kpsi;
  // The prior scale is log sd_k = log tau - alpha log k: an intercept read at
  // k = 1, where no mode of either series lives, and a slope against log k.
  // Intercept and slope correlate exactly as they do in any uncentred
  // regression -- measured at +0.55 for the traversal pair, +0.81 on the
  // training sample -- so log k is centred on its mean, which is what the
  // prior's own information matrix picks out, and tau becomes the scale of a
  // mid-range mode rather than of a mode nothing constrains.
  real kbar_q;
  real kbar_psi;
  array[N_sn] int sn_index;
  for (k in 1:(Kq - 1)) kq[k] = k;
  for (k in 1:(Kpsi - 1)) kpsi[k] = k;
  kbar_q = exp(mean(log(kq)));
  kbar_psi = exp(mean(log(kpsi)));
  for (n in 1:N_sn) sn_index[n] = n;
}
parameters {
  // The departure coefficients non-centred, a_tilde = tau_q k^-alpha_q z_a:
  // the gauge normalizes the functions they carry, so their scale is the
  // prior's alone.  The mean coefficients centred: the data pin them, and
  // the first hierarchical fit showed the non-centred form trading tau
  // against ||z_b|| at correlation -0.6, with 88 percent of transitions at
  // maximum treedepth.
  vector[Kq - 1] z_a;
  vector[Kq - 1] b_til;
  vector[Kpsi - 1] z_c;
  vector[Kpsi - 1] d_til;
  real<lower=0> tau_q;
  real<lower=0> tau_psi;
  real<lower=0> alpha_q_m;         // alpha - 1/2, Eq. (coefprior)
  real<lower=0> alpha_psi_m;
  real<lower=0> sigma_g;
  real<lower=0> sigma_r;
  real mu_p;
  real mu_psi;
  real mu_delta;
  real<lower=0> s_delta;
  vector[N_sn] m_peak;
  vector[N_sn] colr;
  vector[N_sn] x_s;
  // The amplitudes with sampled locations, as deviations about them: the
  // first fit's longest directions were the common modes of the per-SN
  // block riding with their locations (corr(mu_psi, mean x_psi) = +1.00)
  // and the s_delta funnel (corr(log s_delta, log rms(delta_t - mu_delta))
  // = +0.95); written this way each common mode is one coordinate.
  vector[N_sn] u_p;                // x_p = mu_p + u_p
  vector[N_sn] u_psi;              // x_psi = mu_psi + u_psi
  vector[N_sn] u_del;              // delta_t = mu_delta + s_delta u_del
}
transformed parameters {
  vector[N_sn] x_p = mu_p + u_p;
  vector[N_sn] x_psi = mu_psi + u_psi;
  vector[N_sn] delta_t = mu_delta + s_delta * u_del;
  // Eq. (gauge): a, b under the flat measure in v; c, d under ds through G_s.
  vector[Kq - 1] a;
  vector[Kq - 1] b;
  vector[Kpsi - 1] c;
  vector[Kpsi - 1] d;
  {
    vector[Kq - 1] at = pow(kq, -(alpha_q_m + 0.5)) .* z_a;
    vector[Kpsi - 1] ct = pow(kpsi, -(alpha_psi_m + 0.5)) .* z_c;
    a = at / sqrt(dot_self(at));
    b = b_til - dot_product(b_til, a) * a;
    c = ct / sqrt(quad_form(Gs, ct));
    d = d_til - (d_til' * Gs * c) * c;
  }
}
model {
  z_a ~ std_normal();
  z_c ~ std_normal();
  // Eq. (coefprior) on the centred pair; tau_q appears nowhere else for a,
  // whose scale the gauge removes, so it is b_til and d_til that inform it.
  b_til ~ normal(0, tau_q * pow(kq / kbar_q, -(alpha_q_m + 0.5)));
  d_til ~ normal(0, tau_psi * pow(kpsi / kbar_psi, -(alpha_psi_m + 0.5)));
  alpha_q_m ~ exponential(1);      // Section 5: alpha - 1/2 ~ Exp(1)
  alpha_psi_m ~ exponential(1);
  // tau, sigma, mu, s_delta flat on their supports
  u_p ~ std_normal();              // x_p ~ N(mu_p, 1)
  u_psi ~ std_normal();            // x_psi ~ N(mu_psi, 1)
  u_del ~ std_normal();            // delta_t ~ N(mu_delta, s_delta)
  x_s ~ normal(xs_bar, 1);
  m_peak ~ normal(m_bright, 1);
  colr ~ normal(0, 1);
  target += reduce_sum(partial_sum, sn_index, grainsize,
                       dt, f, sig_f, band, start, stop, z,
                       m_peak, colr, x_s, x_p, x_psi, delta_t, sigma_g, sigma_r,
                       a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                       Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                       edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
}
generated quantities {
  real alpha_q = alpha_q_m + 0.5;
  real alpha_psi = alpha_psi_m + 0.5;
  // The two sign conventions of Section 2.3, as joint relabelings.
  // int_0^{v_max} q_1 dv on the fixed nodes right of zero: the panels past
  // edge m+1 in the two-sided layout, with the plain weights.
  int flip_q = (dot_product(w_v[(m * J + 1):(2 * m * J)],
                            (Xv_nodes * a)[(m * J + 1):(2 * m * J)]) < 0);
  int flip_psi = (turning_row(guhi, Kpsi, ulo, uhi, P0, T) * c < 0);
  vector[Kq - 1] a_rel = flip_q ? -a : a;
  vector[Kpsi - 1] c_rel = flip_psi ? -c : c;
  vector[N_sn] x_p_rel = flip_q ? -x_p : x_p;
  vector[N_sn] x_psi_rel = flip_psi ? -x_psi : x_psi;
  real mu_p_rel = flip_q ? -mu_p : mu_p;
  real mu_psi_rel = flip_psi ? -mu_psi : mu_psi;
  real kappa0 = turning_deriv_row(0, Kpsi, ulo, uhi, T) * d / delta;
  // Per draw, per supernova: log-likelihood, the interval check, and the
  // softening diagnostics of Section 5.2 item 3.
  vector[N_sn] lp_sn;
  vector[N_sn] s_pmin;
  vector[N_sn] s_pmax;
  real frac_tapered = 0;
  real frac_clamped = 0;
  for (n in 1:N_sn) {
    int lo = start[n];
    int hi = stop[n];
    matrix[hi - lo + 1, 4] lf = log_flux(dt[lo:hi], band[lo:hi], z[n],
        m_peak[n], colr[n], x_s[n], x_p[n], x_psi[n], delta_t[n],
        a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
        Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
        edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
    lp_sn[n] = -0.5 * neg2_loglike(f[lo:hi], sig_f[lo:hi], band[lo:hi], exp(lf[:, 1]),
                                   sigma_g, sigma_r);
    for (i in 1:(hi - lo + 1)) {
      frac_tapered += (lf[i, 2] < log(0.99));
      frac_clamped += (abs(lf[i, 4] - lf[i, 3]) > eps_s * (smax - smin));
    }
    {
      vector[2] ends = arclength([pmin + 5, pmax - 5]', x_s[n], b + x_p[n] * a, Kq, vlo, vhi,
                                 delta_p, m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      s_pmin[n] = ends[1];
      s_pmax[n] = ends[2];
    }
  }
  frac_tapered /= N;
  frac_clamped /= N;
}
