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
  int<lower=1> grainsize;
  // A fixed rotation of the coefficient coordinates, identity unless a
  // preconditioner has been built.  All of the posterior's stiffness lives
  // in these: of the 226 coordinates of a thirty-object fit the seven
  // stiffest directions put 96 to 100 percent of their weight here, and
  // whitening them takes the condition number from 2.4e5 to 7.0e3 -- the
  // steps a diagonal metric needs from 485 to 84 -- while whitening the
  // per-object blocks changes nothing.  R is constant, so its Jacobian is
  // constant and the posterior is untouched: this is a reparameterization.
  matrix[2 * (Kq - 1) + 2 * (Kpsi - 1), 2 * (Kq - 1) + 2 * (Kpsi - 1)] R;
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
  vector[2 * (Kq - 1) + 2 * (Kpsi - 1)] w_coef;
  real<lower=0> tau_q;
  real<lower=0> tau_psi;
  real<lower=0> alpha_q_m;         // alpha - 1/2, Eq. (coefprior)
  real<lower=0> alpha_psi_m;
  real mu_p;
  real mu_psi;
  // The three placement parameters whose population means were fitted
  // Gaussians' are written as a mean and deviations summing to zero, which
  // is Eq. (xsmean)'s device and assumes nothing about how the objects are
  // distributed.  See hier_common.stan for the full argument and the
  // measurements; in brief: m_peak's prior read the brightest good g epoch,
  // which is the same photometry the likelihood reads, and the box replaces
  // it with the zp = 30 convention that zp and m_peak are degenerate in;
  // colr's N(0,1) asserted a mean anchor colour of zero where three fits put
  // it at +0.060; and <Delta> = 0 is a gauge, a common shift of Delta being
  // a shift of the phase origin that q_0 reshapes to absorb, measured at
  // 2.19 nats over 1568 epochs at (K_psi, K_q) = (24, 16).
  vector<lower=5, upper=30>[N_sn] m_peak;
  real<lower=-5, upper=5> c_bar;
  sum_to_zero_vector[N_sn] c_dev;
  // The mean log rate as a named global, and the per-supernova part summing
  // to zero.  The traversal basis drops the constant Legendre term, so
  // qbar_0 = 0 identically and q_0 cannot carry the mean: it has to live in
  // x_s, and Section 5 used to centre x_s's prior on a \salt\ track's value,
  // -3.0.  Measured, that prior was both misplaced and idle --- the
  // population mean came out -2.85 and -2.70 at the two orders, 2.3 and 3.1
  // sigma away, while the realized spread of x_s was 0.12 to 0.25 against a
  // prior sd of 1, and no object of thirty drew as much as a tenth of its
  // precision from it.  Naming the mean and constraining the deviations to
  // sum to zero removes the template constant and assumes nothing about how
  // the population is distributed, where a fitted mean with a normal prior
  // would assume it is Gaussian.  sum_to_zero_vector is an orthogonal
  // transform carrying N_sn - 1 unconstrained coordinates, so it adds no
  // degeneracy of its own.
  // The mean log rate as q_0's own k = 0 coefficient, and the per-supernova
  // rate as deviations summing to zero.  b_0 restores the degree-zero
  // Legendre term that traversal_row's [2:Kq] slice drops, P_0 = 1/sqrt(vhi -
  // vlo), and is held outside the rotated block because the same row builds
  // q_1, whose gauge requires qbar_1 = 0.  (b_0, x_s) carry an exact gauge --
  // only their sum reaches the exponent -- so sum(x_s) = 0 fixes it at no
  // cost.  The bound is written on the mean log rate so it survives a change
  // of phase window, and it is needed: both tails are flat, at -4.83e6 frozen
  // and -2.93e6 clamped.  See hier_common.stan for the full argument.
  real<lower = -20 * sqrt(vhi - vlo), upper = 10 * sqrt(vhi - vlo)> b_0;
  sum_to_zero_vector[N_sn] x_s;
  // The amplitudes with sampled locations, as deviations about them: the
  // first fit's longest direction was the common mode of the per-SN block
  // riding with its location (corr(mu_psi, mean x_psi) = +1.00), and
  // written this way each common mode is one coordinate.  The second of
  // those directions was the s_delta funnel, which went with the population
  // of Delta; delta_t below is constrained instead, so there is no scale
  // left to fall into.
  vector[N_sn] u_p;                // x_p = mu_p + u_p
  vector[N_sn] u_psi;              // x_psi = mu_psi + u_psi
  // Delta flat on a box, per object, with no parameter for the mean; see
  // hier_common.stan for the measurements.  In brief: <Delta> is a near
  // gauge worth 1.96 nats over 1568 epochs, the box confines it rather than
  // deleting it, the bound is what makes the flat prior proper -- both tails
  // are flat, the phase clamp saturating the model's dependence on Delta
  // entirely -- and +-3 d is 8.6 times the largest t0 uncertainty DR2
  // reports over the thirty.  The cost is that the common mode, the long
  // direction of this block at a ratio of 3.3 to 4.0 over the deviations at
  // N = 30, is left mixed across all thirty coordinates where a diagonal
  // metric cannot scale it on its own.
  vector<lower=-3, upper=3>[N_sn] delta_t;
}
transformed parameters {
  // what reaches the exponent: the deviation plus q_0's degree-zero term
  vector[N_sn] x_s_tot = x_s + b_0 * inv_sqrt(vhi - vlo);
  vector[N_sn] x_p = mu_p + u_p;
  vector[N_sn] x_psi = mu_psi + u_psi;
  vector[N_sn] colr = c_bar + c_dev;
  // the coefficients the priors are written on, read out of the rotation
  vector[Kq - 1] z_a;
  vector[Kq - 1] b_til;
  vector[Kpsi - 1] z_c;
  vector[Kpsi - 1] d_til;
  // Eq. (gauge): a, b under the flat measure in v; c, d under ds through G_s.
  vector[Kq - 1] a;
  vector[Kq - 1] b;
  vector[Kpsi - 1] c;
  vector[Kpsi - 1] d;
  {
    vector[2 * (Kq - 1) + 2 * (Kpsi - 1)] raw = R * w_coef;
    int n1 = Kq - 1;
    int n3 = Kpsi - 1;
    z_a = raw[1:n1];
    b_til = raw[(n1 + 1):(2 * n1)];
    z_c = raw[(2 * n1 + 1):(2 * n1 + n3)];
    d_til = raw[(2 * n1 + n3 + 1):(2 * n1 + 2 * n3)];
  }
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
  // tau, sigma and the two mu flat on their supports
  u_p ~ std_normal();              // x_p ~ N(mu_p, 1)
  u_psi ~ std_normal();            // x_psi ~ N(mu_psi, 1)
  // Propriety, not population structure.  Section 5 records that with a flat
  // prior the likelihood tends to a positive constant as x_s -> -infinity,
  // where the supernova freezes at peak, so some proper prior in that
  // direction is required however uninformative it is in the bulk.  Neither of
  // these is a population model: no distribution over the x_s is inferred, and
  // the widths are deliberately far looser than the data -- s_bar is pinned to
  // about 0.07 and the realized spread of the deviations is 0.12 to 0.25, so
  // the normal below is four to eight times too wide to shape either.  What
  // they replace is a prior centred on a \salt\ track's -3.0, which was 2.3
  // to 3.1 sigma from where the data put the mean.
  // b_0 and x_s are flat; x_s needs a wall because both ends of a hyperplane
  // ray saturate into constant light curves, the density flat at -2.47e5.  A
  // declared bound is impossible on a sum_to_zero_vector, so this is a
  // sampling statement, contributing exactly 0 inside.
  x_s ~ uniform(-20, 20);
  // c_dev has no prior and needs none: colour is the only placement direction
  // the forward map does not clamp, so a runaway inside the sum-to-zero
  // hyperplane drives a predicted r flux to infinity where the Gaussian-in-
  // flux likelihood kills it, and the flat prior is proper on its own.  See
  // hier_common.stan for the measurement.  m_peak, c_bar and delta_t are flat
  // on their boxes.
  target += reduce_sum(partial_sum, sn_index, grainsize,
                       dt, f, sig_f, band, start, stop, z,
                       m_peak, colr, x_s_tot, x_p, x_psi, delta_t,
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
        m_peak[n], colr[n], x_s_tot[n], x_p[n], x_psi[n], delta_t[n],
        a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
        Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
        edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
    lp_sn[n] = -0.5 * neg2_loglike(f[lo:hi], sig_f[lo:hi], exp(lf[:, 1]));
    for (i in 1:(hi - lo + 1)) {
      frac_tapered += (lf[i, 2] < log(0.99));
      frac_clamped += (abs(lf[i, 4] - lf[i, 3]) > eps_s * (smax - smin));
    }
    {
      vector[2] ends = arclength([pmin + 5, pmax - 5]', x_s_tot[n], b + x_p[n] * a, Kq, vlo, vhi,
                                 delta_p, m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      s_pmin[n] = ends[1];
      s_pmax[n] = ends[2];
    }
  }
  frac_tapered /= N;
  frac_clamped /= N;
}
