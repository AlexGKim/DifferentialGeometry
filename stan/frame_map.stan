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
  // One stored draw of hier_rot, read in to be mapped into the frame
  // coordinates and mapped back.  z_a_in and z_c_in rebuild the normalised
  // partners exactly as the sampler does, so the frame is the same frame.
  vector[Ka - 1] z_a_in;
  vector[Kc] z_c_in;
  vector[Kq - 1] b_in;             // the projected traversal vector
  vector[Kpsi] de_in;              // (d_0, d), the projected turning vector
  vector[Kq - 1] b_til_in;         // the ambient vectors, before projection
  vector[Kpsi] d_til_in;
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
  // The comoving frame carries one coefficient fewer per pair: the component
  // the projection used to discard does not exist.
  int n_bf = n_b - 1;
  int n_df = n_d - 1;
  // G = L L' for each pair, as the transposed Cholesky factor and its
  // inverse, which is what frame_lift needs.  Declared here and filled at the
  // foot of this block, because Gvq is not assigned until the anchor is.
  matrix[Kq - 1, Kq - 1] LtV;
  matrix[Kq - 1, Kq - 1] LtVinv;
  matrix[Kpsi, Kpsi] LtS;
  matrix[Kpsi, Kpsi] LtSinv;
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
  LtV = cholesky_decompose(Gvq)';
  LtVinv = inverse(LtV);
  LtS = cholesky_decompose(Gse)';
  LtSinv = inverse(LtS);
}
parameters {
}
generated quantities {
  vector[Kq - 2] y_b;
  vector[Kpsi - 1] y_d;
  real res_b;                      // ||frame_lift(frame_drop(b)) - b||_inf
  real res_d;
  real gauge_b;                    // b' Gvq a, zero by construction
  real gauge_d;                    // de' Gse ce
  real norm_b;                     // ||b||_Gvq  vs  ||y_b||
  real norm_d;
  real drop_b;                     // the prior terms rot_frame no longer has
  real drop_d;
  real dlp_pred;                   // lp(rot_frame) - lp(hier_rot), predicted
  // y from the PROJECTION of the stored ambient vector, which is what
  // hier_rot itself forms, so the two programs carry the same b and de to
  // machine precision rather than to the CSV's six figures.
  vector[Kq - 2] y_b_p;
  vector[Kpsi - 1] y_d_p;
  real round_b;                    // ||b_proj - b_in||_inf, the CSV rounding
  real round_d;
  {
    vector[Kq - 1] at = rep_vector(0, Kq - 1);
    vector[Kpsi] ct = rep_vector(0, Kpsi);
    vector[Kq - 1] a_s;
    vector[Kpsi] ce_s;
    at[1:(Ka - 1)] = z_a_in;
    ct[1:Kc] = z_c_in;
    a_s = at / sqrt(quad_form(Gvq, at));
    ce_s = ct / sqrt(quad_form(Gse, ct));
    y_b = frame_drop(LtV, a_s, b_in);
    y_d = frame_drop(LtS, ce_s, de_in);
    res_b = max(abs(frame_lift(LtV, LtVinv, a_s, y_b) - b_in));
    res_d = max(abs(frame_lift(LtS, LtSinv, ce_s, y_d) - de_in));
    gauge_b = b_in' * Gvq * a_s;
    gauge_d = de_in' * Gse * ce_s;
    norm_b = sqrt(quad_form(Gvq, b_in)) - sqrt(dot_self(y_b));
    norm_d = sqrt(quad_form(Gse, de_in)) - sqrt(dot_self(y_d));
    drop_b = std_normal_lpdf(b_til_in' * Gvq * a_s)
             + log(sqrt(dot_self(Gvq * a_s)));
    drop_d = std_normal_lpdf(d_til_in' * Gse * ce_s)
             + log(sqrt(dot_self(Gse * ce_s)));
    dlp_pred = -(drop_b + drop_d);
    {
      vector[Kq - 1] bp = b_til_in - (b_til_in' * Gvq * a_s) * a_s;
      vector[Kpsi] dp = d_til_in - (d_til_in' * Gse * ce_s) * ce_s;
      y_b_p = frame_drop(LtV, a_s, bp);
      y_d_p = frame_drop(LtS, ce_s, dp);
      round_b = max(abs(bp - b_in));
      round_d = max(abs(dp - de_in));
    }
  }
}
