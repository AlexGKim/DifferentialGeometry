// hier_full with Eq. (basispsi)'s constant modes: Psi_1(0) = c_0 and
// Psi_0(0) = d_0 are fitted rather than held at zero.
//
// Why.  The note carried Psi_0(0) = Psi_1(0) = 0 inside Eq. (curvlin) and
// attributed the turning pair's missing additive constants to them.  They fix
// nothing: a constant in Psi_1 puts x_psi*mu into psi, a per-object rotation
// of the track, and m_peak and colr span only translations, so the
// transformation is not a symmetry and there is nothing to gauge away.  The
// conditions were restrictions on the family, and Psi_1(0) = 0 in particular
// asserted that the g maximum sits at the same arclength origin for every
// supernova.  Measured, a per-object anchor rotation is the best-identified
// direction in the per-object sector -- irreducible fraction 0.47 to 0.57
// against 0.17 for m_peak and colr and 0.04 to 0.16 for Delta -- so the
// capacity declined is capacity the data can use.
//
// What changes.  z_c and d_til each gain one coefficient, index 1 being the
// constant mode; Eq. (gauge)'s two turning conditions are imposed on the
// extended vectors under turning_gram_ext, which is turning_gram_s bordered
// by the constant.  Through Eq. (curvlin) the constants reach the likelihood
// only as the rotation alpha_n = d_0 + x_psi_n c_0 of that object's track
// about the anchor, applied in log_flux_rot.  c_0 = d_0 = 0 is hier_full.
//
// Both conditions stay inner products and neither becomes a point condition
// at the anchor: Psi_0(0) = 0 would fix nu only while c_0 != 0 and so
// degenerates at the nested model, and Psi_1'(0) = 1 is unusable because
// Psi_1'(0) sits 0.2 sigma from zero and changes sign in a third of
// consecutive draws of ladder-full-k812-n120.  The sign of c is therefore
// still unfixed, as in hier_full, and is relabelled in generated quantities.
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
  // Prior sd of Delta's deviations in days, Section 5.  Data rather than a
  // literal so that a run is identified by its data.json, and fixed rather
  // than fitted so that no scale can run away.
  real<lower=0> d_sd;
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
  // Kc and Kpsi rather than Kc - 1 and Kpsi - 1: Eq. (basispsi)'s constant
  // modes put one more coefficient in each turning block.
  matrix[(Ka - 1) + (Kq - 1) + Kc + Kpsi,
         (Ka - 1) + (Kq - 1) + Kc + Kpsi] R;
  // Eq. (gauge)'s traversal conditions: 0 anchors the pair at v = 0, as every
  // run before 2026-11-24 did, and 1 anchors its mean under the measure the
  // epochs impose.  Kept as a switch rather than a rewrite so that the
  // earlier parameterisation stays samplable and the two can be compared on
  // one data dict.  See traversal_qbar in dgsn_functions.stan.
  int<lower=0, upper=1> qbar_on;
  int<lower=1> qbar_npan;          // panels on [0, v_i]; converged at 2
}
transformed data {
#include dgsn_grid_tdata.stan
  // Eq. (gauge)'s traversal pair, anchored at the measure's centre when
  // qbar_on.  The phases are built at delta_t = 0 and clamped exactly as
  // log_flux_rot clamps them, so Qbar is a function of the data dict alone.
  vector[N] v_ep0;
  row_vector[Kq - 1] Qbar;
  row_vector[Kq - 1] anchor_v;
  matrix[Kq - 1, Kq - 1] Gvq;
  // Eq. (gauge) on the extended turning pair; Gs is its lower-right block.
  matrix[Kpsi, Kpsi] Gse = turning_gram_ext(Kpsi, ulo, uhi, gulo, guhi, delta,
                                            P0, T, 200, xg, wg);
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
  // Eq. (basispsi): index 1 is the constant mode, 2..  the Ptilde_k.
  int n_c = Kc;
  int n_d = Kpsi;
  // d(f)/f = (ln10 / 2.5) dm, so a model error of sigma_b magnitudes is a
  // fractional flux error of 0.4 ln10 sigma_b.  Eq. (modelerr).
  real MAG2FRAC = 0.4 * log(10);
  for (k in 1:(Kq - 1)) kq[k] = k;
  for (k in 1:(Kpsi - 1)) kpsi[k] = k;
  kbar_q = exp(mean(log(kq)));
  kbar_psi = exp(mean(log(kpsi)));
  for (n in 1:N_sn) sn_index[n] = n;
  for (n in 1:N_sn)
    for (i in start[n]:stop[n]) {
      real wph = (dt[i] / (1 + z[n]) - pmin) / (pmax - pmin);
      v_ep0[i] = v_of_p(pmin + (pmax - pmin) * clamp01(wph, eps), delta_p);
    }
  Qbar = qbar_on == 1
         ? traversal_qbar(v_ep0, Kq, vlo, vhi, qbar_npan, xg, wg, P0v)
         : rep_row_vector(0, Kq - 1);
  anchor_v = P0v + Qbar;
  // Gv of dgsn_grid_tdata.stan is this at shift = P0v, which is the test
  // that the gauge still means <q_1,q_1> = 1 and <q_0,q_1> = 0 under the
  // same integral.  Checked numerically before this program was run.
  Gvq = traversal_gram_shift(Kq, vlo, vhi, anchor_v, 4, xg, wg);
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
  //
  // Eq. (ampmean): each amplitude as a named mean plus deviations summing to
  // zero, both flat, exactly as colr and x_s are already written and for the
  // same reason.  Written without the split, the mean of a block of N
  // amplitudes is one direction inside it and a diagonal metric is obliged to
  // give it the scale it gives the deviations.  Measured as
  // var(<x>) * N / mean(var(x_n)), which is 1 when that scale is right: the
  // unsplit amplitudes read 6.3 to 8.2 over sixty objects and 2.4 to 46.2 over
  // 255, against exactly 0 for every block already carrying the split.  At 255
  // the unsplit <x_p> was the slowest direction in the problem -- integrated
  // autocorrelation time 28.8 draws against 2.0 for c_bar -- with split R-hat
  // 2.48 where its own deviations sat at 1.07, 34 to 68 percent of transitions
  // divergent and the median leapfrog pinned at the treedepth-12 cap.
  //
  // The box is now on each part and not on the sum, so the reachable set is
  // not quite the single box's: |x_p| reached 4.4 and |x_psi| 9.6 against 20,
  // and neither wall has been approached.  The deviations take a wall in the
  // model block, sum_to_zero_vector admitting no declared bound.
  //
  // This buys conditioning, not identification.  <q_0, q_1> = 0 fixes mu_x_p,
  // so its value states the gauge and not the supernovae, and whether q_1
  // earns its place is read from the spread of the deviations and from
  // prediction -- never from the mean.
  real<lower=-20, upper=20> mu_x_p;
  sum_to_zero_vector[N_sn] x_p_dev;
  real<lower=-20, upper=20> mu_x_psi;
  sum_to_zero_vector[N_sn] x_psi_dev;
  // Delta as mu_delta + deviations summing to zero; see hier_common.stan for
  // the measurements.  <Delta> is a near gauge worth 1.96 nats over 1568
  // epochs, and the mean is named so a diagonal metric can scale the common
  // mode on its own -- measured as the slowest-mixing direction when it was
  // not.
  //
  // The deviations take a broad normal, N(0, d_sd^2), where through
  // 2026-11-16 they took a flat wall at +-3 d.  A prior they do need: with
  // the mean held, a hyperplane ray saturates the phase clamp and the density
  // goes flat, so flat and unbounded is improper.  What the wall did badly was
  // reject where it should have restored.  Measured on stage3-full-k812-n255
  // and stage3-ampmean-k88-n60, d_dev was the only block pressed against its
  // wall -- 99.7 percent of it, against 59 for the two runs that converged --
  // and the only one whose approach tracked divergence, corr +0.52.  Per
  // object the wall sat 1.7 to 4.3 sigma away for the sparsely sampled and 53
  // to 65 for the well covered, and a crossing mid-trajectory returns -inf,
  // which the sampler records as a divergence.  Delta is not unidentified
  // there -- sigma(Delta) is 0.15 to 0.71 d against 0.045 to 0.056 for the
  // well covered -- the posterior simply wants 1.8 to 2.5 d, which is where
  // the wall was.
  //
  // d_sd = 3 d is 3.2 times the realized spread of 0.929 d, so it shapes
  // nothing: 0.25 to 5.6 percent of the curvature for the sparse objects and
  // 0.02 to 0.04 for the well covered.  What it adds is a restoring gradient
  // where the wall gave none -- chi2 of 1.0 at 3 d and 4.0 at 6 d, the
  // pushback holding tmax near the data that Section 4 says the taper gives
  // up, and which a WIDER box would have weakened instead.  Scaling it by the
  // DR2 t0 uncertainty was measured and rejected: t0_err runs to 0.816 d over
  // the 255 with a median of 0.129, roughly ten times too small, because
  // Delta is the B-to-g maximum offset and not a t0 measurement error.
  //
  // The scale is DATA and not fitted.  Fitting it is the (s_delta, u)
  // degeneracy of Section 5, whose supremum is at s_delta -> infinity.  On a
  // sum_to_zero_vector the marginal sd is d_sd sqrt((N-1)/N), 0.2 percent at
  // N = 255.
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
  // Eq. (ampmean), assembled the same way delta_t is
  vector[N_sn] x_p = mu_x_p + x_p_dev;
  vector[N_sn] x_psi = mu_x_psi + x_psi_dev;
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
  // Eq. (basispsi): Psi_1(0) and Psi_0(0), fitted.
  real c_0;
  real d_0;
  // Eq. (curvlin): the constants reach the likelihood only here, as the
  // rotation of each object's track about the anchor.
  vector[N_sn] alpha;
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
    vector[Kpsi] ct = rep_vector(0, Kpsi);
    at[1:n_a] = z_a;
    ct[1:n_c] = z_c;
    a = at / sqrt(quad_form(Gvq, at));
    b = b_til - (b_til' * Gvq * a) * a;
    {
      // Eq. (gauge) on the extended vectors, under turning_gram_ext.
      vector[Kpsi] ce = ct / sqrt(quad_form(Gse, ct));
      vector[Kpsi] de = d_til - (d_til' * Gse * ce) * ce;
      c_0 = ce[1];
      d_0 = de[1];
      c = ce[2:Kpsi];
      d = de[2:Kpsi];
    }
  }
  // Eq. (curvlin): the rotation of each track about the anchor.  At
  // c_0 = d_0 = 0 this is zero and log_flux_rot reproduces log_flux.
  alpha = d_0 + x_psi * c_0;
}
model {
  z_a ~ std_normal();
  z_c ~ std_normal();
  // Eq. (coefprior) on b and d -- the coefficients that satisfy Eq. (gauge)
  // -- and NOT on b_til and d_til.  Section 5 has the argument; in brief,
  // b_til's component along a and d_til's along c are discarded by the
  // orthogonalisation and reach neither the forward map nor the likelihood,
  // so under the old form their only curvature was Eq. (coefprior)'s own,
  // whose width grows with tau.  That put sqrt(cond) proportional to tau
  // along directions carrying no information: measured, removing the four
  // such directions took sqrt(cond) from 5.62 to 1.00 and from 9.84 to 1.38
  // at median draws of the two mode-bearing runs.
  //
  // Each pair therefore contributes three terms, and only the first survives
  // from the old form:
  //   * Eq. (coefprior) evaluated on the gauge-satisfying coefficient;
  //   * Eq. (hypernorm)'s normalisation, because that coefficient lives on a
  //     hyperplane whose orientation is itself a parameter, so the constant
  //     cannot be dropped -- it is +0.5 log(sum_k n_k^2 w_k^2) with n the
  //     UNIT normal, verified against a brute-force orthonormal-basis
  //     evaluation to 1e-8 at both orders and at random w.  It looks like a
  //     prior on the mode direction and is the opposite: it is what makes
  //     p(d_til | c) integrate to one at every c, so that c's marginal prior
  //     stays the projected normal z_c carries.  Dropping it tilts c by 1.8x
  //     over six prior draws and 6 to 8 across the prior's range;
  //   * the Jacobian of b_til -> (b, b_til.a), which is unity for the
  //     traversal pair since a is Euclidean-unit, and ||Gs c|| for the
  //     turning pair since c is Gs-normalised instead, so the decomposition
  //     is oblique in the coefficient basis.  |det[B | c]| = 1/||Gs c||,
  //     verified numerically at three orders, and -log|det| = +log||Gs c||.
  // The discarded components take unit width.  Their priors are
  // unobservable, so the choice is free, and unit is what makes the
  // precision along them independent of tau.
  //
  // None of these terms reaches the likelihood, so the transfer gate of
  // Section 5 -- sum ll_sn unchanged at the carried point -- cannot see an
  // error in any of them.  They are checked by sampling with the likelihood
  // switched off; see scripts/prior_check.py.
  {
    // b and d are flat; only the gauge is declared.  Each pair contributes
    // two terms instead of three:
    //
    //   * unit width on the discarded component, b_til . a and
    //     d_til' Gs c.  These reach neither the forward map nor the
    //     likelihood, so a prior on them is required for propriety and is
    //     otherwise unobservable; unit width is what makes the precision
    //     along them independent of everything else.
    //
    //   * the Jacobian of the decomposition.  For the traversal pair the map
    //     b_til -> (b, b_til . a) is orthogonal, a being Euclidean-unit, so
    //     the Jacobian is 1 and nothing is added.  For the turning pair the
    //     decomposition is oblique in the coefficient basis,
    //     |det[B | c]| = 1/||Gs c||, so -log|det| = +log||Gs c||.
    //
    // Eq. (hypernorm) is NOT here and must not be: it normalized a Gaussian
    // restricted to a hyperplane whose orientation is a parameter, and with
    // b and d flat there is no Gaussian to normalize.  The flat measure on
    // the hyperplane is its own (K-2)-dimensional Lebesgue measure, which
    // does not depend on the orientation, so no orientation-dependent
    // constant is owed.  The Jacobian is owed either way and stays.
    target += std_normal_lpdf(dot_product(b_til, a));
    {
      // the extended turning vector of Eq. (basispsi); Gse replaces Gs
      // throughout because d_til and the normalised mode now carry c_0.
      vector[Kpsi] ce = append_row(c_0, c);
      target += std_normal_lpdf(d_til' * Gse * ce)
                + log(sqrt(dot_self(Gse * ce)));
    }
  }
  // Half-normal, since sigma_b is <lower=0>.  0.1 mag puts 95% of the mass
  // below 0.196, so it covers the 0.0469 and 0.0679 the no-mode fit at sixty
  // objects measured, and room above them, without asserting either.
  sigma_b ~ normal(0, 0.1);
  // x_p and x_psi are flat on their declared boxes and take no statement.
  // Under Eq. (ampmean) the means carry that box as a declared bound and the
  // deviations take it as a wall here, sum_to_zero_vector admitting no
  // declaration -- the same arrangement x_s and d_dev already use.
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
  // Widened from +-3 on 2026-11-16.  Measured on stage3-full-k812-n255 and
  // stage3-ampmean-k88-n60, d_dev was the only block pressed against its
  // wall -- 99.7 percent of it, against 59 percent for the converged runs --
  // and the only one whose approach tracked divergence (corr +0.52).  Per
  // object the wall sat 1.7 to 4.3 sigma away for the sparsely sampled and
  // 53 to 65 sigma for the well covered, every mid-trajectory crossing
  // returning -inf and being recorded as a divergence.  At +-6 the same four
  // objects sit 6 to 23 sigma inside.
  //
  // This weakens the second job the bound does.  The taper deweights an
  // out-of-window epoch rather than penalising it, so once Delta pushes
  // epochs out nothing pulls it back and the bound is the pushback; Section 4
  // records the window already crossed below the old bound for sixteen of
  // thirty objects.  frac_tapered is the diagnostic and must be read against
  // the +-3 run, not assumed unchanged.
  d_dev ~ normal(0, d_sd);   // see the declaration; mu_delta is flat on its box
  x_p_dev ~ uniform(-20, 20);     // Eq. (ampmean); mu_x_p flat on its box
  x_psi_dev ~ uniform(-20, 20);   // Eq. (ampmean); mu_x_psi flat on its box
  // c_dev has no prior and needs none: colour is the only placement direction
  // the forward map does not clamp, so a runaway inside the sum-to-zero
  // hyperplane drives a predicted r flux to infinity where the Gaussian-in-
  // flux likelihood kills it, and the flat prior is proper on its own.  See
  // hier_common.stan for the measurement.  m_peak, c_bar and delta_t are flat
  // on their boxes.
  target += reduce_sum(partial_sum_floor_rot, sn_index, grainsize,
                       dt, f, sig_f, MAG2FRAC * sigma_b[band], band,
                       start, stop, z,
                       m_peak, colr, x_s_tot, x_p, x_psi, delta_t, alpha,
                       a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                       Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                       edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T,
                       anchor_v);
}
generated quantities {
  // The two sign conventions of Section 2.3, as joint relabelings.
  // int_0^{v_max} q_1 dv on the fixed nodes right of zero: the panels past
  // edge m+1 in the two-sided layout, with the plain weights.
  // The branch label is the sign of the MEAN AMPLITUDE, not of a functional
  // of the mode.  The old label was the sign of int_0^vmax q_1 dv, and that
  // integral sits near zero: measured on runs/stage3-modes-k88-n60 it
  // toggled 12 to 48 times per chain while the raw mean amplitude made zero
  // sign changes in 253 to 336 draws, agreeing with the actual branch only
  // 56 to 82 percent of the time.  Every _rel column it fed was therefore
  // being emitted with a random sign.  mean(x_p) cannot have that problem
  // while it stays away from zero, which is the condition to watch: it ran
  // +2.18 to +2.48 at sd 0.35, so about six sigma clear.
  real x_p_bar = mean(x_p);
  real x_psi_bar = mean(x_psi);
  int flip_q = (x_p_bar < 0);
  int flip_psi = (x_psi_bar < 0);
  vector[Kq - 1] a_rel = flip_q ? -a : a;
  vector[Kpsi - 1] c_rel = flip_psi ? -c : c;
  vector[N_sn] x_p_rel = flip_q ? -x_p : x_p;
  vector[N_sn] x_psi_rel = flip_psi ? -x_psi : x_psi;
  real x_p_bar_rel = flip_q ? -x_p_bar : x_p_bar;
  real x_psi_bar_rel = flip_psi ? -x_psi_bar : x_psi_bar;
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
    matrix[hi - lo + 1, 4] lf = log_flux_rot(dt[lo:hi], band[lo:hi], z[n],
        m_peak[n], colr[n], x_s_tot[n], x_p[n], x_psi[n], delta_t[n], alpha[n],
        a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
        Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
        edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T,
        anchor_v);
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
      vector[2] ends = arclength([pmin + 5, pmax - 5]', x_s_tot[n], b + x_p[n] * a, Kq, vlo, vhi,
                                 delta_p, m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      s_pmin[n] = ends[1];
      s_pmax[n] = ends[2];
    }
  }
  frac_tapered /= N;
  frac_clamped /= N;
}
