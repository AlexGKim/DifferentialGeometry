// hier_base.stan -- hier_full.stan with both departure modes removed, and
// nothing else changed.  One function describes the traversal and one the
// turning: log ds/dp = x_s + b_0 P_0 + q_0(p) and psi = -pi/2 + Psi_0(s).
//
// The pair is nested, and it is worth being exact about in what sense.  The
// FORWARD MAP is nested identically: log_flux reads b + x_p a and d + x_psi c,
// so at x_p = x_psi = 0 it reads b and d, which is what this program passes.
// Every b this program can reach, hier_full can reach too -- its b is required
// only to be orthogonal to a, and a ranges over the whole sphere.  So the two
// likelihoods are the same function of the same q_0 and Psi_0.
//
// The POSTERIOR DENSITY is not the same, and cannot be: hier_full carries
// 2(Ka-1) + 2(Kc-1) more coordinates, a unit normal on each of the two
// components its gauge discards, and the oblique Jacobian log||G_s c||.  None
// of those reach the forward map, which is why they do not disturb the
// comparison, but a log density from one program is not comparable to a log
// density from the other without them.  ll_sn and chi2_sn are.
//
// So a comparison between the two is a test of the departure modes and of
// nothing else -- not of the coefficient prior, which neither carries, and
// not of how the amplitudes are parameterised, which only hier_full decides.
//
// Eq. (18)'s fitted per-band model error, the generated-quantity triple that
// separates what the target receives from the archive-error chi2, and the
// per-epoch residuals are all hier_full's and are kept.
//
// Not to be confused with hier_common.stan, which also has one function each
// but KEEPS Eq. (coefprior) and its two scales and two decay exponents.  That
// program is the model-v1.2 reference runs/stage3-abox-k88-n60; this one has
// no prior on b and d at all.
//
// Two orders per pair, not one.  Kq and Kpsi are the orders of q_0 and
// Psi_0, the mean functions; Ka and Kc are the orders of q_1 and Psi_1, the
// departure modes, and may be lower.  The reason is not the coordinate count
// -- halving them saves 14 coordinates of 420 at (12, 16) on sixty objects,
// about 3 percent -- but detectability.  A mode's shape is identified only
// through the spread of its per-object amplitudes, and that spread is read
// against the per-object posterior sd, which grows with the mode's own
// order: measured over the runs on disk at thirty objects, x_psi has sd
// 0.691 at Kpsi = 12 against a spread of 0.987, and sd 1.519 at Kpsi = 24
// against a spread of 0.628.  At the higher order the amplitudes drown and
// the effective number of objects informing the shape falls to zero.
// Keeping Ka and Kc low keeps the amplitudes measurable while Kq and Kpsi
// give the mean functions the resolution they need.
//
// free_mean is deliberately absent, so <Delta> is always sampled.  It is not
// a gauge: s(p = 0) = 0 pins the anchor, so a common shift of Delta is a
// real shift of the phase origin rather than a reparameterisation.
// runs/stage3-k88-dzero prices the constraint at 8.9 nats of sum ln L at the
// best draw and 6.6 at the mean, and shows it buys no mixing -- both arms
// realise treedepth 8, at sqrt(cond) 91 against 119.  No fit has put it at
// zero: -0.3951 +- 0.0959 at thirty objects, +0.4118 +- 0.0428 and
// +0.4345 +- 0.1979 at sixty.  The thirty- and sixty-object samples disagree
// on its sign, which is unexplained, and is a reason to leave the coordinate
// free rather than to fix it at either value.
//
// The likelihood is a reduce_sum over supernovae.  Per-epoch residuals are
// saved here, which hier.stan's header rules out for stage 4 and which is
// affordable at 3682 epochs: they are the only thing that says *where* a
// departure mode helps, a per-band sigma_b absorbing phase-structured error
// uniformly and so reporting a good fit either way.
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
  // The departure modes' own orders, Section 2.2's q_1 and Psi_1.  Ka = Kq
  // and Kc = Kpsi is the full model; lower is a truncation of the mode
  // alone, which by Cauchy interlacing cannot raise the condition number of
  // the retained block and which keeps the per-object amplitudes measurable.
  // Bounded below at 2 because a mode with no coefficient cannot be
  // normalised: the gauge needs ||a|| = 1, and a/||a|| is undefined at zero.

  // A fixed rotation of the coefficient coordinates, identity unless a
  // preconditioner has been built.  All of the posterior's stiffness lives
  // in these: of the 226 coordinates of a thirty-object fit the seven
  // stiffest directions put 96 to 100 percent of their weight here, and
  // whitening them takes the condition number from 2.4e5 to 7.0e3 -- the
  // steps a diagonal metric needs from 485 to 84 -- while whitening the
  // per-object blocks changes nothing.  R is constant, so its Jacobian is
  // constant and the posterior is untouched: this is a reparameterization.
  matrix[(Kq - 1) + (Kpsi - 1), (Kq - 1) + (Kpsi - 1)] R;
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
  // The four slice widths of the rotated block, in the order the program
  // reads them: z_a, b_til, z_c, d_til.
  int n_b = Kq - 1;
  // log_flux's signature is shared with the mode-bearing programs, so the
  // absent modes are passed as exact zeros rather than by a second signature.
  vector[Kq - 1] a_zero = rep_vector(0, Kq - 1);
  vector[Kpsi - 1] c_zero = rep_vector(0, Kpsi - 1);
  vector[N_sn] x_p_zero = rep_vector(0, N_sn);
  vector[N_sn] x_psi_zero = rep_vector(0, N_sn);
  int n_d = Kpsi - 1;
  // d(f)/f = (ln10 / 2.5) dm, so a model error of sigma_b magnitudes is a
  // fractional flux error of 0.4 ln10 sigma_b.  Eq. (modelerr).
  real MAG2FRAC = 0.4 * log(10);
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
  vector[n_b + n_d] w_coef;
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
  // No tau, no alpha, and no mu.  Eq. (coefprior)'s hierarchy is gone: on
  // every direction the likelihood can see it supplied between 1/4600 and
  // 1/284000 of the curvature, measured at two draws of
  // runs/stage3-modes-k88-n60 with the four gauge-redundant directions
  // projected out, so it was not earning its place.  What remains is the
  // gauge itself, below in the model block: unit width on the two discarded
  // components and the Jacobian of the oblique turning decomposition.  Those
  // four directions are unobservable and are the only ones a prior is
  // required for.
  //
  // With the k^-alpha weighting gone, a and c are uniform on their spheres,
  // which is what the fit chose anyway -- their effective decay exponents
  // measured -0.070 +- 0.441 and +0.119 +- 0.116 against a shared alpha of
  // 0.698 that was pressed against its own floor, P(alpha - 1/2 < 0.10) =
  // 0.359 against the Exp(1) prior's 0.095.
  //
  // mu_p and mu_psi go with the amplitudes' hierarchy; see x_p below.
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
  // The departure amplitudes, flat on a box, like every other per-supernova
  // parameter here: x_s and d_dev are flat on theirs, c_dev needs no prior,
  // and m_peak, c_bar and delta_t take boxes.  These two were the only ones
  // carrying an informative population prior, x ~ N(mu, 1) with the scale
  // fixed at 1 -- and the gauge already uses up the (q_1, x_p) -> (lam q_1,
  // x_p/lam) redundancy through ||a|| = 1, so that 1 was an assertion about
  // the mode's strength rather than a normalization.  The data disputed it in
  // every chain: sd(x_p) came out 0.849 to 0.864 and sd(x_psi) 0.879 to
  // 0.896, with no parameter able to absorb the difference.  Flat leaves the
  // whole distribution of amplitudes to the data.
  //
  // The box is needed for x_p and is prudence for x_psi.  Driving all sixty
  // objects out together, x_p's chi2 PLATEAUS -- two objects return
  // 79246.3 and 21851.4 at x_p = -500, -100 and +500 alike, the supernova
  // frozen by the clamps -- so a flat prior over that tail would not
  // integrate, exactly as Section 5 records for x_s.  x_psi does not
  // saturate, its chi2 still climbing 81 percent over the last decade of the
  // same sweep, but the sweep is numerically chaotic out there and x_p
  // overflows to NaN by 2000, so both take the wall.  Fitted ranges are
  // -0.126 to +3.783 and +2.268 to +7.303, so 20 is three times outside and
  // matches the box x_s already carries.
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
  // Eq. (modelerr): the fractional model error, in magnitudes, one per band
  // with g = 1 and r = 2 as `band` indexes them.  It is not measurement
  // error -- sig_f already carries the sky and the source -- but the gap
  // between the model curve and the supernova's own, and it is fitted
  // because its size is not known.  Here it carries a second job: a
  // departure mode that is real signal must take structure out of that gap,
  // so sigma_b falling below the no-mode fit's 0.0469 and 0.0679 mag is the
  // calibrated test of whether the mode earns its place, in magnitudes
  // rather than in nats.
  vector<lower=0>[2] sigma_b;
}
transformed parameters {
  vector[N_sn] delta_t = mu_delta + d_dev;
  // what reaches the exponent: the deviation plus q_0's degree-zero term
  vector[N_sn] x_s_tot = x_s + b_0 * inv_sqrt(vhi - vlo);
  vector[N_sn] colr = c_bar + c_dev;
  // the coefficients the priors are written on, read out of the rotation
  vector[n_b] b_til;
  vector[n_d] d_til;
  // Eq. (gauge) is VACUOUS here.  It fixes the split between a mean function
  // and a departure mode by making b orthogonal to a and d orthogonal to c
  // under G_s; with no a and no c there is nothing to be orthogonal to, so
  // b = b_til and d = d_til identically, and the four gauge-redundant
  // coordinates the mode-bearing programs carry do not exist.
  vector[Kq - 1] b;
  vector[Kpsi - 1] d;
  {
    vector[n_b + n_d] raw = R * w_coef;
    b_til = raw[1:n_b];
    d_til = raw[(n_b + 1):(n_b + n_d)];
  }
  b = b_til;
  d = d_til;
}
model {
  // b and d carry NO prior.  In hier_full the only density on either was a
  // unit normal on the component the gauge discards -- b_til along a, d_til
  // along c -- and those components do not exist without the modes.  Nor is
  // anything owed in their place: there is no gauge redundancy left, so every
  // direction of b and d reaches the forward map and is pinned by the
  // likelihood alone.  Eq. (coefprior) is deliberately absent, as in
  // hier_full.stan; hier_common.stan is the program that keeps it.
  // Half-normal, since sigma_b is <lower=0>.  0.1 mag puts 95% of the mass
  // below 0.196, so it covers the 0.0469 and 0.0679 the no-mode fit at sixty
  // objects measured, and room above them, without asserting either.
  sigma_b ~ normal(0, 0.1);
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
  target += reduce_sum(partial_sum_floor, sn_index, grainsize,
                       dt, f, sig_f, MAG2FRAC * sigma_b[band], band,
                       start, stop, z,
                       m_peak, colr, x_s_tot, x_p_zero, x_psi_zero, delta_t,
                       a_zero, b, c_zero, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                       Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                       edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
}
generated quantities {
  // No sign flags and no relabeled copies: the two discrete symmetries
  // (q_1, x_p) -> (-q_1, -x_p) and (Psi_1, x_psi) -> (-Psi_1, -x_psi) go with
  // the modes, so this program has no sign freedom to record.
  real kappa0 = turning_deriv_row(0, Kpsi, ulo, uhi, T) * d / delta;
  // Per draw, per supernova: log-likelihood, the interval check, and the
  // softening diagnostics of Section 5.2 item 3.
  // lp_sn keeps its historical meaning exactly -- -chi2/2 on the archive
  // errors, no floor and no log-determinant -- so the figures recorded for
  // every earlier run stay comparable with this one.  ll_sn is what the
  // target receives under Eq. (18), and chi2_sn is the residuals over V,
  // which approaches one per epoch when the floor is doing its job.
  vector[N_sn] lp_sn;
  vector[N_sn] ll_sn;
  vector[N_sn] chi2_sn;
  vector[N_sn] s_pmin;
  vector[N_sn] s_pmax;
  // Per epoch, the standardised residual and the arclength it sits at.  A
  // per-band sigma_b absorbs phase-structured error uniformly, so a fitted
  // floor says the error budget is honest and nothing about whether the
  // curve is right; these two are what localise the misfit, and arclength
  // rather than phase because that is the variable Psi_1 acts in and it
  // cannot be reconstructed outside Stan.
  vector[N] resid_z;
  vector[N] s_epoch;
  real frac_tapered = 0;
  real frac_clamped = 0;
  for (n in 1:N_sn) {
    int lo = start[n];
    int hi = stop[n];
    matrix[hi - lo + 1, 4] lf = log_flux(dt[lo:hi], band[lo:hi], z[n],
        m_peak[n], colr[n], x_s_tot[n], 0.0, 0.0, delta_t[n],
        a_zero, b, c_zero, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
        Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
        edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
    lp_sn[n] = -0.5 * neg2_loglike(f[lo:hi], sig_f[lo:hi], exp(lf[:, 1]), exp(lf[:, 2]));
    {
      vector[hi - lo + 1] sfr = MAG2FRAC * sigma_b[band[lo:hi]];
      vector[hi - lo + 1] V = square(sig_f[lo:hi]) + square(sfr .* f[lo:hi]);
      ll_sn[n] = -0.5 * neg2_loglike_floor(f[lo:hi], sig_f[lo:hi],
                                           exp(lf[:, 1]), exp(lf[:, 2]), sfr);
      chi2_sn[n] = sum(exp(lf[:, 2]) .* square(f[lo:hi] - exp(lf[:, 1])) ./ V);
      resid_z[lo:hi] = (f[lo:hi] - exp(lf[:, 1])) ./ sqrt(V);
      s_epoch[lo:hi] = lf[:, 3];
    }
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
