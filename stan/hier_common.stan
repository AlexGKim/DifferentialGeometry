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
  // Eq. (xsmean)'s device, applied to every placement parameter whose
  // population mean was previously a fitted Gaussian's.  None of the three
  // below is a population model, and nothing is assumed about how the
  // objects are distributed; each is a constraint on the realized values.
  //
  // m_peak flat on a box.  Its prior was N(m_bright, 1) with m_bright the
  // brightest good g epoch, which is a centre read out of the same
  // photometry that enters the likelihood --- the data twice.  The box
  // replaces it: zp and m_peak are exactly degenerate (Section 3), so the
  // magnitude origin is a convention already fixed at zp = 30, and the
  // bounds are that convention's own statement rather than a judgement
  // about supernovae.  m_peak = 30 is a predicted flux of one archive unit
  // and m_peak = 5 is 1e10, against a sample spanning 1.5e4 to 7.8e5, so
  // the box is four decades clear on both sides and flat everywhere the
  // data are.  It is there only to make the posterior proper: with a flat
  // prior on the whole line the likelihood tends to exp(-sum f^2/2 sig^2)
  // as m_peak -> +infinity, a positive constant, and the integral diverges.
  vector<lower=5, upper=30>[N_sn] m_peak;
  // The anchor colour's mean named, its deviations summing to zero.  The
  // old prior N(0, 1) asserted a mean anchor colour of zero; three runs put
  // it at +0.060 +- 0.013 (N = 30) and +0.0597 +- 0.0032 (N = 100), so the
  // centre was misplaced, in the way Section 5 records for x_s's templated
  // -3.0.  c_bar is flat on a box for the same reason m_peak is.
  real<lower=-5, upper=5> c_bar;
  sum_to_zero_vector[N_sn] c_dev;
  // The mean log rate as q_0's own k = 0 coefficient, and the per-supernova
  // rate as deviations summing to zero.
  //
  // q_0 is expanded in Legendre degrees 1..Kq-1 and so contains no constant;
  // b_0 restores the degree-zero term, P_0 = 1/sqrt(vhi - vlo).  It is held
  // here rather than inside the rotated coefficient block for one reason: the
  // same traversal_row builds q_1, and Eq. (gauge) requires qbar_1 = 0, so
  // widening the slice to [1:Kq] would hand the departure mode a constant and
  // break its gauge.  Since P_0 is constant, b_0 * P_0 is a constant in the
  // exponent and the two forms are identical -- b_0 is the k = 0 coefficient,
  // stored outside the vector R mixes.  It cannot take Eq. (coefprior) in any
  // case: (k/kbar)^-alpha diverges at k = 0 and kbar = exp(mean log k) is
  // undefined there.
  //
  // (b_0, x_s) carry an EXACT gauge: x_s -> x_s + eps with
  // b_0 -> b_0 - eps sqrt(vhi - vlo) leaves every prediction identical, since
  // only their sum reaches the exponent.  sum(x_s) = 0 fixes it exactly, and
  // unlike <Delta> = 0 there is nothing to pay for doing so.  What this
  // replaces, s_bar + x_s_dev, is the same model relabelled; what it buys is
  // that q_0 is now the population's log rate outright, mean included, rather
  // than a shape with its mean parked in the per-supernova block.
  //
  // The bound is on the mean log rate, not on b_0, so it survives a change of
  // phase window: b_0/sqrt(vhi - vlo) in [-20, 10] against a posterior near
  // -2.9.  It is needed.  Both tails are flat -- downward the rate vanishes
  // and the track freezes at the anchor, settling at -4.83e6; upward it
  // diverges, every epoch is driven into the arclength clamp, settling at
  // -2.93e6 -- so no width of flat prior would have been proper.
  real<lower = -20 * sqrt(vhi - vlo), upper = 10 * sqrt(vhi - vlo)> b_0;
  sum_to_zero_vector[N_sn] x_s;
  // Delta: flat on a box, per object, with no parameter for the mean.
  //
  // <Delta> is a near gauge and not an exact one.  A common shift of every
  // Delta is a shift of the phase origin, which q_0 -- a free function of v
  // -- reshapes to absorb, so the mean is pinned by the basis's finite
  // flexibility rather than by the data: the per-object posterior sd is
  // 0.359 d, thirty independent offsets would put 0.359/sqrt(30) = 0.066 d
  // on their mean, and a hierarchical fit returned 0.344 d, 5.2 times
  // larger.  An earlier version fixed <Delta> = 0 and deleted the direction
  // at a measured cost of 1.96 nats over 1568 epochs; the box confines it
  // instead, so the fit finds its own phase origin, keeps those nats, and
  // leaves the B-to-g offset in Delta rather than in q_0's minimum.
  //
  // The bound is what makes the flat prior proper, and it is needed: with
  // Delta unbounded the likelihood does not fall away.  Past |Delta| ~ 100 d
  // every epoch is outside the phase window, the softened clamp saturates,
  // the model stops depending on Delta at all, and the log density reads
  // -156709.38 at 300 d and -156709.38 at 1000 d.  Both tails are flat, so
  // no amount of width would have sufficed.
  //
  // Why +-3 d.  Over 1800 draws at N = 30 and N = 100 the largest |Delta|
  // ever drawn was 2.43 d and the largest deviation about the mean 1.085 d.
  // It is also 8.6 times the largest t0 uncertainty DR2 reports over the
  // thirty (0.350 d) and 3.7 times the largest over the training sample
  // (0.816 d), so like m_peak's box it is a statement about the archive and
  // not about supernovae.
  //
  // What this costs.  The common mode is the long direction of this block,
  // sd 0.343 and 0.437 at N = 30 against a mean deviation sd of 0.104 and
  // 0.109 -- a ratio of 3.3 and 4.0 -- and 1.64 at N = 100.  A diagonal
  // metric assigns one scale per declared coordinate, and every coordinate
  // here mixes the common mode with the deviations, so that ratio is
  // conditioning the metric has to straddle.  Naming the mean would have
  // split them; it is deliberately not named, and chain pace is therefore
  // the thing to watch.  Delta is the only one of the four per-supernova
  // parameters where the ratio exceeds one: for m_peak, colr and x_s the
  // common mode is the shorter direction.
  vector<lower=-3, upper=3>[N_sn] delta_t;
}
transformed parameters {
  // what reaches the exponent: the per-supernova deviation plus q_0's own
  // degree-zero term, P_0 = 1/sqrt(vhi - vlo)
  vector[N_sn] x_s_tot = x_s + b_0 * inv_sqrt(vhi - vlo);
  vector[N_sn] colr = c_bar + c_dev;
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
  // Propriety, as Section 5 requires: with a flat prior the likelihood tends
  // to a positive constant as x_s -> -infinity, where the supernova freezes
  // at peak.  Neither of these is a population model.
  // b_0 and x_s are flat: the box on b_0 is the whole prior, and x_s needs
  // a wall because both ends of a hyperplane ray saturate -- one object's
  // rate clamped, another's frozen, each a constant light curve, the
  // density flat at -2.47e5 from t = 50 outward.  The wall sits at 270
  // times the realized deviation sd of 0.073 and is never reached; a
  // declared bound is impossible on a sum_to_zero_vector, so it is a
  // sampling statement, which contributes exactly 0 inside.
  x_s ~ uniform(-20, 20);
  // c_dev has no prior at all, and needs none.  It is the one block here
  // whose flat prior is already proper: inside the sum-to-zero hyperplane
  // some component must run to +infinity, which drives that object's
  // predicted r flux up exponentially, and a Gaussian-in-flux likelihood
  // punishes that doubly exponentially -- walked out, the log density falls
  // through -4.67e6 at t = 3 and -2.05e12 at t = 10 to -inf, and on the least
  // favourable ray available (the object with five r epochs rising against
  // the one with a hundred and thirteen falling) it diverges faster still.
  // Colour is the only placement direction the forward map does not clamp,
  // which is why this works here and nowhere else: arclength and phase
  // saturate into constant light curves, which are finite and therefore flat.
  // The N(0,1) this replaces was shrinking by 8.3e-4, so no number moves.
  //
  // m_peak, c_bar and delta_t are flat on their boxes.
  target += reduce_sum(partial_sum, sn_index, grainsize,
                       dt, f, sig_f, band, start, stop, z,
                       m_peak, colr, x_s_tot, x_p_zero, x_psi_zero, delta_t,
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
        m_peak[n], colr[n], x_s_tot[n], 0, 0, delta_t[n],
        a_zero, b, c_zero, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
        Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
        edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
    lp_sn[n] = -0.5 * neg2_loglike(f[lo:hi], sig_f[lo:hi], exp(lf[:, 1]));
    for (i in 1:(hi - lo + 1)) {
      frac_tapered += (lf[i, 2] < log(0.99));
      frac_clamped += (abs(lf[i, 4] - lf[i, 3]) > eps_s * (smax - smin));
    }
    {
      vector[2] ends = arclength([pmin + 5, pmax - 5]', x_s_tot[n], b, Kq, vlo, vhi,
                                 delta_p, m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      s_pmin[n] = ends[1];
      s_pmax[n] = ends[2];
    }
  }
  frac_tapered /= N;
  frac_clamped /= N;
}
