// The full six-parameter model of Section 2.2, with Eq. (18)'s model error.
// hier.stan's hierarchy -- four coefficient vectors, both departure modes,
// both gauge conditions -- plus the three things hier_common.stan gained
// afterwards and hier.stan never did: the fitted per-band model error, the
// generated-quantity triple that separates what the target receives from the
// archive-error chi2, and per-epoch residuals.
//
// Why this program exists rather than an edit to hier.stan.  hier.stan and
// hier_common.stan are what every run on record was sampled by, and
// hier_common.stan is the converged reference runs/stage3-k88-n60-floor.
// Both stay untouched so both stay reproducible; this file is where the two
// departure modes are explored.
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
  int<lower=2, upper=Kq> Ka;
  int<lower=2, upper=Kpsi> Kc;
  // A fixed rotation of the coefficient coordinates, identity unless a
  // preconditioner has been built.  All of the posterior's stiffness lives
  // in these: of the 226 coordinates of a thirty-object fit the seven
  // stiffest directions put 96 to 100 percent of their weight here, and
  // whitening them takes the condition number from 2.4e5 to 7.0e3 -- the
  // steps a diagonal metric needs from 485 to 84 -- while whitening the
  // per-object blocks changes nothing.  R is constant, so its Jacobian is
  // constant and the posterior is untouched: this is a reparameterization.
  matrix[(Ka - 1) + (Kq - 1) + (Kc - 1) + (Kpsi - 1),
         (Ka - 1) + (Kq - 1) + (Kc - 1) + (Kpsi - 1)] R;
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
  int n_a = Ka - 1;
  int n_b = Kq - 1;
  int n_c = Kc - 1;
  int n_d = Kpsi - 1;
  // d(f)/f = (ln10 / 2.5) dm, so a model error of sigma_b magnitudes is a
  // fractional flux error of 0.4 ln10 sigma_b.  Eq. (modelerr).
  real MAG2FRAC = 0.4 * log(10);
  for (k in 1:(Kq - 1)) kq[k] = k;
  for (k in 1:(Kpsi - 1)) kpsi[k] = k;
  kbar_q = exp(mean(log(kq)));
  kbar_psi = exp(mean(log(kpsi)));
  for (n in 1:N_sn) sn_index[n] = n;

  // ---- a rest-frame phase grid for the curve figures (this program only) ----
  // One synthetic mean object: x_s = 0 (the sum-to-zero population mean),
  // z = 0 so observer time IS rest-frame phase, delta_t = 0.  Both bands.
  int ng = 161;
  vector[2 * ng] dt_fig;
  array[2 * ng] int band_fig;
  for (i in 1:ng) {
    real ph = pmin + (pmax - pmin) * (i - 1.0) / (ng - 1.0);
    dt_fig[i]      = ph;  band_fig[i]      = 1;
    dt_fig[ng + i] = ph;  band_fig[ng + i] = 2;
  }

}
parameters {
  // The departure coefficients non-centred, a_tilde = tau_q k^-alpha_q z_a:
  // the gauge normalizes the functions they carry, so their scale is the
  // prior's alone.  The mean coefficients centred: the data pin them, and
  // the first hierarchical fit showed the non-centred form trading tau
  // against ||z_b|| at correlation -0.6, with 88 percent of transitions at
  // maximum treedepth.
  vector[n_a + n_b + n_c + n_d] w_coef;
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
  // to zero.  The traversal basis drops the constant and anchors each mode at v = 0, so
  // q_0(0) = 0 identically and q_0 cannot carry the mean: it has to live in
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
  // q_1, whose gauge requires q_1(0) = 0.  (b_0, x_s) carry an exact gauge --
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
  vector<lower=-20, upper=20>[N_sn] x_p;
  vector<lower=-20, upper=20>[N_sn] x_psi;
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
  vector[n_a] z_a;
  vector[n_b] b_til;
  vector[n_c] z_c;
  vector[n_d] d_til;
  // Eq. (gauge): a, b under the flat measure in v; c, d under ds through G_s.
  vector[Kq - 1] a;
  vector[Kq - 1] b;
  vector[Kpsi - 1] c;
  vector[Kpsi - 1] d;
  {
    vector[n_a + n_b + n_c + n_d] raw = R * w_coef;
    z_a = raw[1:n_a];
    b_til = raw[(n_a + 1):(n_a + n_b)];
    z_c = raw[(n_a + n_b + 1):(n_a + n_b + n_c)];
    d_til = raw[(n_a + n_b + n_c + 1):(n_a + n_b + n_c + n_d)];
  }
  {
    // The modes are carried at the mean functions' length, zero past their
    // own order, so log_flux's signature and every shared function in
    // dgsn_functions.stan are untouched by the truncation.  Padding with
    // zeros is exact: the Legendre rows nest, and turning_transform's
    // Cholesky factor nests with them.  The gauge sees only the nonzero
    // entries -- ||a|| = 1 and c' Gs c = 1 hold on the padded vectors
    // because the padding contributes nothing to either form.
    vector[Kq - 1] at = rep_vector(0, Kq - 1);
    vector[Kpsi - 1] ct = rep_vector(0, Kpsi - 1);
    at[1:n_a] = z_a;
    ct[1:n_c] = z_c;
    a = at / sqrt(quad_form(Gv, at));
    b = b_til - (b_til' * Gv * a) * a;
    c = ct / sqrt(quad_form(Gs, ct));
    d = d_til - (d_til' * Gs * c) * c;
  }
}
model {
}
generated quantities {
  // Magnitudes of one mean object on a rest-frame phase grid, with the two
  // departure modes switched on and off.  Everything is log_flux's, so the
  // basis, the quadrature and both clamps are the fitted program's own.
  real x_p_bar   = mean(x_p);
  real x_psi_bar = mean(x_psi);
  real x_p_sd    = sd(x_p);
  real x_psi_sd  = sd(x_psi);
  real kappa0 = turning_deriv_row(0, Kpsi, ulo, uhi, T) * d / delta;

  // The reference is the population mean object: BOTH modes at their mean
  // amplitude.  x_p = 0 is not "no departure" -- the gauge puts the whole
  // amplitude of q_1 into x_p, whose population mean is about +3, so the
  // mode-free point sits some two object-sds away from any real supernova.
  // mag_zero is kept only to show that distance.
  vector[2 * ng] mag_ref;     // (x_p_bar,          x_psi_bar)
  vector[2 * ng] mag_qlo;     // (x_p_bar - x_p_sd, x_psi_bar)
  vector[2 * ng] mag_qhi;     // (x_p_bar + x_p_sd, x_psi_bar)
  vector[2 * ng] mag_plo;     // (x_p_bar,          x_psi_bar - x_psi_sd)
  vector[2 * ng] mag_phi;     // (x_p_bar,          x_psi_bar + x_psi_sd)
  vector[2 * ng] mag_zero;    // (0, 0), the mode-free point
  {
    real mp = mean(m_peak);
    real cl = c_bar;                 // colr = c_bar + c_dev, and c_dev sums to zero
    // log_flux's sixth argument is x_s_tot, NOT x_s: the mean log rate
    // b_0 / sqrt(vhi - vlo) is carried outside the basis and has to be added
    // back, or exp(x_s_tot + q_0) is wrong by a factor of order e^8 and the
    // arclength pins against its clamp at both ends of the window.
    real xst = mean(x_s_tot);
    real dtm = mean(delta_t);
    // mag = zp - log_flux_col1 / (0.4 ln10); col 1 is -0.4 ln10 (mag - zp).
    real inv = 1.0 / (0.4 * log(10));
    #include full_curve_calls.stan
  }
}
