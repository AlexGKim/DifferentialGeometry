// The hierarchy of Section 5.2 with the turning mode removed: Eq. (curvlin)
// reduced from psi = -pi/2 + Psi_0 + x_psi Psi_1 to psi = -pi/2 + Psi_0, so
// every supernova shares one curve in magnitude space and differs only in
// where it is placed and how it traverses.
//
// Why this program exists.  At K_psi = 24 the Fisher information on c has
// four directions above ten times the prior out of 23, falling as k^-2.79,
// because c reaches the observable only through an integral of Psi_1 along
// the Frenet normal; the per-supernova x_psi has posterior sd 1.52 against a
// unit prior at that order and 0.91 at K_psi = 12, with the spread of object
// means only 0.76.  So the mode may never have earned its place, and the test
// is whether removing it costs likelihood.
//
// Three things go with it.  The gauge on the turning pair is gone: with no c
// to normalize there is no cTGc = 1, and d = d_til rather than d orthogonalized
// against c, so Psi_0 is simply a function the data determine.  The sign
// freedom (Psi_1, x_psi) -> (-Psi_1, -x_psi) is gone, leaving one discrete
// symmetry instead of two.  And mu_psi and u_psi go, with 23 + 1 + 30
// coordinates removed in all at (24, 16) on thirty supernovae.
//
// tau_psi and alpha_psi_m stay: they are the prior scale and decay of d_til,
// which is Psi_0's own expansion, and have nothing to do with the mode.
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
  // As hier.stan, but over the shorter block: z_a, b_til, d_til, with no z_c.
  matrix[2 * (Kq - 1) + (Kpsi - 1), 2 * (Kq - 1) + (Kpsi - 1)] R;
}
transformed data {
#include dgsn_grid_tdata.stan
  vector[Kq - 1] kq;
  vector[Kpsi - 1] kpsi;
  real kbar_q;
  real kbar_psi;
  array[N_sn] int sn_index;
  // Psi_1 and its amplitudes, held at zero as constants rather than as
  // transformed parameters, so they cost nothing per draw and no dead column
  // reaches the output.
  vector[Kpsi - 1] c_zero = rep_vector(0, Kpsi - 1);
  vector[N_sn] x_psi_zero = rep_vector(0, N_sn);
  for (k in 1:(Kq - 1)) kq[k] = k;
  for (k in 1:(Kpsi - 1)) kpsi[k] = k;
  kbar_q = exp(mean(log(kq)));
  kbar_psi = exp(mean(log(kpsi)));
  for (n in 1:N_sn) sn_index[n] = n;
}
parameters {
  vector[2 * (Kq - 1) + (Kpsi - 1)] w_coef;
  // Eq. (ccoef)'s scales, sampled under a half-normal that acts as a tail
  // bound and nothing more.  Three treatments have now been tried.  Flat on
  // (0, inf) leaves nothing to bound tau's excursions and gives the ridge
  // this posterior's worst direction: corr(log tau_psi, log||d||) = +0.82,
  // ESS 18 for ||d|| over 1500 draws.  Holding tau as data
  // (runs/stage3-common-fixtau) raised that to 71 and tightened the stepsize
  // spread to 1.7x, but stiffened the scale -- Rhat on ||d|| 1.004 -> 1.273,
  // kappa0 7.8 -- because departures from a held scale cost prior density
  // quadratically where before tau could move with them.  The half-normal is
  // the interpolation: tau still breathes with ||d||, but its excursions are
  // curved rather than free.
  //
  // The scale is 10 for both, which is a tail bound and not a fit to the
  // observed amplitudes.  tau_q sits at 0.81 and tau_psi at 2.24, so the
  // prior's gradient is of order tau/100 against a likelihood term of order
  // (K-1)/tau ~ 10 -- a few parts in a thousand, so it does not shape either
  // in the bulk.  Setting it anywhere near the observed values would, and
  // would also put either a template or the data's own answer inside the
  // prior, which Section 5 removed for x_s and m_peak respectively.
  real<lower=0> tau_q;
  real<lower=0> tau_psi;
  real<lower=0> alpha_q_m;         // alpha - 1/2, Eq. (coefprior)
  real<lower=0> alpha_psi_m;
  real mu_p;
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
  vector[N_sn] u_p;                // x_p = mu_p + u_p
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
  // Delta as mu_delta + deviations summing to zero; see hier_common.stan.  The
  // mean is named so a diagonal metric can scale the common mode on its own,
  // measured as the slowest-mixing direction when it was not.  mu_delta takes
  // a declared bound and the deviations a wall, sum_to_zero_vector admitting
  // none, and both are needed: with the mean held, a hyperplane ray saturates
  // the phase clamp and the density goes flat.
  real<lower=-3, upper=3> mu_delta;
  sum_to_zero_vector[N_sn] d_dev;
}
transformed parameters {
  vector[N_sn] delta_t = mu_delta + d_dev;
  // what reaches the exponent: the deviation plus q_0's degree-zero term
  vector[N_sn] x_s_tot = x_s + b_0 * inv_sqrt(vhi - vlo);
  vector[N_sn] x_p = mu_p + u_p;
  vector[N_sn] colr = c_bar + c_dev;
  vector[Kq - 1] z_a;
  vector[Kq - 1] b_til;
  vector[Kpsi - 1] d_til;
  // Eq. (gauge) on the traversal pair only; the turning pair has no gauge
  // left to impose, Psi_0 being alone.
  vector[Kq - 1] a;
  vector[Kq - 1] b;
  vector[Kpsi - 1] d;
  {
    vector[2 * (Kq - 1) + (Kpsi - 1)] raw = R * w_coef;
    int n1 = Kq - 1;
    z_a = raw[1:n1];
    b_til = raw[(n1 + 1):(2 * n1)];
    d_til = raw[(2 * n1 + 1):(2 * n1 + Kpsi - 1)];
  }
  {
    vector[Kq - 1] at = pow(kq, -(alpha_q_m + 0.5)) .* z_a;
    a = at / sqrt(dot_self(at));
    b = b_til - dot_product(b_til, a) * a;
    d = d_til;
  }
}
model {
  z_a ~ std_normal();
  b_til ~ normal(0, tau_q * pow(kq / kbar_q, -(alpha_q_m + 0.5)));
  d_til ~ normal(0, tau_psi * pow(kpsi / kbar_psi, -(alpha_psi_m + 0.5)));
  tau_q ~ normal(0, 10);           // half-normal: tau is <lower=0>
  tau_psi ~ normal(0, 10);
  alpha_q_m ~ exponential(1);      // Section 5: alpha - 1/2 ~ Exp(1)
  alpha_psi_m ~ exponential(1);
  u_p ~ std_normal();              // x_p ~ N(mu_p, 1)
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
  d_dev ~ uniform(-3, 3);   // see the declaration; mu_delta is flat on its box
  // c_dev has no prior and needs none: colour is the only placement direction
  // the forward map does not clamp, so a runaway inside the sum-to-zero
  // hyperplane drives a predicted r flux to infinity where the Gaussian-in-
  // flux likelihood kills it, and the flat prior is proper on its own.  See
  // hier_common.stan for the measurement.  m_peak, c_bar and delta_t are flat
  // on their boxes.
  target += reduce_sum(partial_sum, sn_index, grainsize,
                       dt, f, sig_f, band, start, stop, z,
                       m_peak, colr, x_s_tot, x_p, x_psi_zero, delta_t,
                       a, b, c_zero, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                       Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                       edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
}
generated quantities {
  real alpha_q = alpha_q_m + 0.5;
  real alpha_psi = alpha_psi_m + 0.5;
  // One sign convention is left, on the traversal pair; the turning pair's
  // went with Psi_1.
  int flip_q = (dot_product(w_v[(m * J + 1):(2 * m * J)],
                            (Xv_nodes * a)[(m * J + 1):(2 * m * J)]) < 0);
  vector[Kq - 1] a_rel = flip_q ? -a : a;
  vector[N_sn] x_p_rel = flip_q ? -x_p : x_p;
  real mu_p_rel = flip_q ? -mu_p : mu_p;
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
        m_peak[n], colr[n], x_s_tot[n], x_p[n], 0, delta_t[n],
        a, b, c_zero, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
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
