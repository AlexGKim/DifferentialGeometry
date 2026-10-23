// TEMPORARY, for the floor check only: the pre-gauge hier_common, kept so the
// constrained optimum can be compared with the free one on the same data.
// The hierarchy with both departure modes removed: the turning of
// Eq. (curvlin) reduced to psi = -pi/2 + Psi_0 and the traversal of
// Eq. (travlog) to log ds/dp = x_s + q_0, so every supernova walks one curve
// at one rate profile and differs only in brightness, colour, overall speed
// and timing.  Four parameters an object instead of six.
//
// Why this program exists.  hier_nomode.stan tests whether x_psi Psi_1 earns
// its place; this tests the same of x_p q_1.  Table 1 of the note calls x_p
// "timing: departure from dilation", so dropping it asserts that the rate
// profile has one common shape and objects differ only in how fast they run
// along it -- x_s multiplying the whole of exp(q_0).
//
// What goes with the two terms.  The gauge of Eq. (gauge) becomes vacuous on
// both pairs: with no q_1 there is nothing to normalize and nothing for q_0
// to be orthogonal to, so b = b_til, and likewise d = d_til.  Both sign
// freedoms go with their modes --- (q_1, x_p) -> (-q_1, -x_p) and
// (Psi_1, x_psi) -> (-Psi_1, -x_psi) --- so this program has no discrete
// symmetry and needs no relabeling, where hier.stan has two and
// hier_nomode.stan one.  mu_p and u_p go with x_p.
//
// tau_q, alpha_q_m, tau_psi and alpha_psi_m all stay: they are the prior
// scale and decay of b_til and d_til, which are q_0's and Psi_0's own
// expansions and have nothing to do with the modes.
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
  int<lower=1> grainsize;
  // Over the shorter block: b_til and d_til alone, no z_a and no z_c.
  matrix[(Kq - 1) + (Kpsi - 1), (Kq - 1) + (Kpsi - 1)] R;
}
transformed data {
#include dgsn_grid_tdata.stan
  vector[Kq - 1] kq;
  vector[Kpsi - 1] kpsi;
  real kbar_q;
  real kbar_psi;
  array[N_sn] int sn_index;
  // The modes and their amplitudes, held at zero as constants: they cost
  // nothing per draw and no dead column reaches the output.
  vector[Kq - 1] a_zero = rep_vector(0, Kq - 1);
  vector[Kpsi - 1] c_zero = rep_vector(0, Kpsi - 1);
  vector[N_sn] x_p_zero = rep_vector(0, N_sn);
  vector[N_sn] x_psi_zero = rep_vector(0, N_sn);
  for (k in 1:(Kq - 1)) kq[k] = k;
  for (k in 1:(Kpsi - 1)) kpsi[k] = k;
  kbar_q = exp(mean(log(kq)));
  kbar_psi = exp(mean(log(kpsi)));
  for (n in 1:N_sn) sn_index[n] = n;
}
parameters {
  vector[(Kq - 1) + (Kpsi - 1)] w_coef;
  real<lower=0> tau_q;
  real<lower=0> tau_psi;
  real<lower=0> alpha_q_m;         // alpha - 1/2, Eq. (coefprior)
  real<lower=0> alpha_psi_m;
  real mu_delta;
  real<lower=0> s_delta;
  vector[N_sn] m_peak;
  vector[N_sn] colr;
  // Eq. (xsmean): the mean log rate named, the deviations summing to zero.
  real s_bar;
  sum_to_zero_vector[N_sn] x_s_dev;
  vector[N_sn] u_del;              // delta_t = mu_delta + s_delta u_del
}
transformed parameters {
  vector[N_sn] x_s = s_bar + x_s_dev;
  vector[N_sn] delta_t = mu_delta + s_delta * u_del;
  vector[Kq - 1] b_til;
  vector[Kpsi - 1] d_til;
  // No gauge: Eq. (gauge)'s four conditions on the traversal pair and two on
  // the turning pair all involve a mode, and there are none.  q_0 and Psi_0
  // are functions the data determine outright.
  vector[Kq - 1] b;
  vector[Kpsi - 1] d;
  {
    vector[(Kq - 1) + (Kpsi - 1)] raw = R * w_coef;
    b_til = raw[1:(Kq - 1)];
    d_til = raw[Kq:(Kq - 1 + Kpsi - 1)];
  }
  b = b_til;
  d = d_til;
}
model {
  // Eq. (coefprior) on both centred vectors, at the geometric-mean pivot.
  b_til ~ normal(0, tau_q * pow(kq / kbar_q, -(alpha_q_m + 0.5)));
  d_til ~ normal(0, tau_psi * pow(kpsi / kbar_psi, -(alpha_psi_m + 0.5)));
  alpha_q_m ~ exponential(1);      // Section 5: alpha - 1/2 ~ Exp(1)
  alpha_psi_m ~ exponential(1);
  u_del ~ std_normal();            // delta_t ~ N(mu_delta, s_delta)
  // Propriety, as Section 5 requires: with a flat prior the likelihood tends
  // to a positive constant as x_s -> -infinity, where the supernova freezes
  // at peak.  Neither of these is a population model.
  s_bar ~ normal(0, 10);
  x_s_dev ~ normal(0, 1);
  m_peak ~ normal(m_bright, 1);
  colr ~ normal(0, 1);
  target += reduce_sum(partial_sum, sn_index, grainsize,
                       dt, f, sig_f, band, start, stop, z,
                       m_peak, colr, x_s, x_p_zero, x_psi_zero, delta_t,
                       a_zero, b, c_zero, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                       Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                       edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
}
generated quantities {
  real alpha_q = alpha_q_m + 0.5;
  real alpha_psi = alpha_psi_m + 0.5;
  // No flip_q, no flip_psi: both sign freedoms went with their modes, so
  // every quantity below is already in one branch and needs no relabeling.
  real kappa0 = turning_deriv_row(0, Kpsi, ulo, uhi, T) * d / delta;
  vector[N_sn] lp_sn;
  vector[N_sn] s_pmin;
  vector[N_sn] s_pmax;
  real frac_tapered = 0;
  real frac_clamped = 0;
  for (n in 1:N_sn) {
    int lo = start[n];
    int hi = stop[n];
    matrix[hi - lo + 1, 4] lf = log_flux(dt[lo:hi], band[lo:hi], z[n],
        m_peak[n], colr[n], x_s[n], 0, 0, delta_t[n],
        a_zero, b, c_zero, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
        Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
        edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
    lp_sn[n] = -0.5 * neg2_loglike(f[lo:hi], sig_f[lo:hi], exp(lf[:, 1]));
    for (i in 1:(hi - lo + 1)) {
      frac_tapered += (lf[i, 2] < log(0.99));
      frac_clamped += (abs(lf[i, 4] - lf[i, 3]) > eps_s * (smax - smin));
    }
    {
      vector[2] ends = arclength([pmin + 5, pmax - 5]', x_s[n], b, Kq, vlo, vhi,
                                 delta_p, m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      s_pmin[n] = ends[1];
      s_pmax[n] = ends[2];
    }
  }
  frac_tapered /= N;
  frac_clamped /= N;
}
