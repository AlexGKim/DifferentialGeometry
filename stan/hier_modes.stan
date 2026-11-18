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
  vector[N_sn] x_p = mu_p + u_p;
  vector[N_sn] x_psi = mu_psi + u_psi;
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
    at[1:n_a] = pow(kq[1:n_a], -(alpha_q_m + 0.5)) .* z_a;
    ct[1:n_c] = pow(kpsi[1:n_c], -(alpha_psi_m + 0.5)) .* z_c;
    a = at / sqrt(quad_form(Gv, at));
    b = b_til - (b_til' * Gv * a) * a;
    c = ct / sqrt(quad_form(Gs, ct));
    d = d_til - (d_til' * Gs * c) * c;
  }
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
    vector[Kq - 1] wq = tau_q * pow(kq / kbar_q, -(alpha_q_m + 0.5));
    vector[Kpsi - 1] wp = tau_psi * pow(kpsi / kbar_psi, -(alpha_psi_m + 0.5));
    vector[Kpsi - 1] nvec = Gs * c;
    real nnorm = sqrt(dot_self(nvec));
    // the traversal pair: b on {x : a' x = 0}, a already of unit length
    target += -0.5 * dot_self(b ./ wq) - sum(log(wq))
              + 0.5 * log(dot_self(a .* wq));
    target += std_normal_lpdf(dot_product(b_til, a));
    // the turning pair: d on {x : (Gs c)' x = 0}
    target += -0.5 * dot_self(d ./ wp) - sum(log(wp))
              + 0.5 * log(dot_self((nvec / nnorm) .* wp))
              + log(nnorm);
    target += std_normal_lpdf(d_til' * Gs * c);
  }
  // Half-normal, since sigma_b is <lower=0>.  0.1 mag puts 95% of the mass
  // below 0.196, so it covers the 0.0469 and 0.0679 the no-mode fit at sixty
  // objects measured, and room above them, without asserting either.
  sigma_b ~ normal(0, 0.1);
  tau_q ~ normal(0, 10);           // half-normal: tau is <lower=0>
  tau_psi ~ normal(0, 10);
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
                            (Xv_nodes * a)[(m * J + 1):(2 * m * J)])
                - (P0v * a) * vhi < 0);
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
        m_peak[n], colr[n], x_s_tot[n], x_p[n], x_psi[n], delta_t[n],
        a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
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
      vector[2] ends = arclength([pmin + 5, pmax - 5]', x_s_tot[n], b + x_p[n] * a, Kq, vlo, vhi,
                                 delta_p, m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      s_pmin[n] = ends[1];
      s_pmax[n] = ends[2];
    }
  }
  frac_tapered /= N;
  frac_clamped /= N;
}
