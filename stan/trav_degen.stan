// The traversal near-degeneracy, measured: how much of the constant function
// the anchored traversal basis already spans under the measure the epochs
// impose, and what shift removes it.
//
// Why.  Eq. (travlog) is s_n(p) = exp(x_s,n) int_0^{v(p)} exp(qtilde_n) dp',
// with qtilde_n = (b + x_p,n a).(P - P(0)) anchored at v = 0.  Shift the
// common log rate by epsilon and the coefficients by dq:
//
//   dlog s_n(p_i) = epsilon + <dq.(P - P(0))>_i ,
//
// the average being against that epoch's own traversal integrand over
// [0, v_i].  So the shift in x_s is invisible to the extent that the
// constant -epsilon lies in the span of the anchored basis under the pooled
// measure of those per-epoch averages.  Write z_i for the averaged basis
// row, m = <z> and M = <z z'> over the epochs.  The best dq is -eps M^-1 m
// and the leftover is the Schur complement of the constant,
//
//   dlog s / epsilon  ->  1 - m' M^-1 m ,
//
// which is zero when the span contains the constant exactly and one when it
// is orthogonal to it.  R^2 = m' M^-1 m is the number this program exists to
// report.  The inversion is linear algebra on matrices Stan produced, done
// outside, as in gram_empirical.stan.
//
// The reference point.  The per-epoch average is integrand-weighted, and the
// integrand carries exp(qtilde) -- parameter dependent.  Taken here at
// qtilde = 0, the flat-rate reference, so that m is data and not a sampled
// quantity; likewise the phases are built at delta_t = 0.  Neither has to be
// exact: using the shifted basis P - P(0) - m leaves R^2 = Delta' M^-1 Delta
// for an error Delta in m, quadratic and so forgiving.
//
// Two resolutions.  npan1 and npan2 panels on [0, v_i]; the pair must agree
// or the quadrature is the limit rather than the geometry.
//
// fixed_param.  Takes any data dict the hierarchy programs take.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  int<lower=1> N_sn;
  int<lower=0> N;
  vector[N] dt;
  vector[N] f;
  vector<lower=0>[N] sig_f;
  array[N] int<lower=1, upper=2> band;
  array[N_sn] int<lower=1> start;
  array[N_sn] int<lower=1> stop;
  vector[N_sn] z;
  int<lower=1> grainsize;
  real<lower=0> d_sd;
  int<lower=2, upper=Kq> Ka;
  int<lower=2, upper=Kpsi> Kc;
  matrix[(Ka - 1) + (Kq - 1) + Kc + Kpsi,
         (Ka - 1) + (Kq - 1) + Kc + Kpsi] R;
  int<lower=1> npan1;
  int<lower=1> npan2;
}
transformed data {
#include dgsn_grid_tdata.stan
  // Phases at delta_t = 0, clamped exactly as log_flux clamps them.
  vector[N] v_ep;
  for (n in 1:N_sn)
    for (i in start[n]:stop[n]) {
      real w = (dt[i] / (1 + z[n]) - pmin) / (pmax - pmin);
      v_ep[i] = v_of_p(pmin + (pmax - pmin) * clamp01(w, eps), delta_p);
    }
}
generated quantities {
  row_vector[Kq - 1] qbar_1;        // m at npan1 panels
  row_vector[Kq - 1] qbar_2;        // m at npan2 panels
  matrix[Kq - 1, Kq - 1] M_1;       // <z z'> at npan1
  matrix[Kq - 1, Kq - 1] M_2;
  // The same two under the point-value measure gram_empirical.stan uses, so
  // the two readings of "the data's measure" can be told apart.
  row_vector[Kq - 1] qbar_pt;
  matrix[Kq - 1, Kq - 1] M_pt;
  // After the shift by qbar_2: m must come back zero and R^2 with it.  Not a
  // tautology at finite quadrature -- it checks that the shift that goes into
  // the model is the shift this program measured.
  row_vector[Kq - 1] qbar_after;
  matrix[Kq - 1, Kq - 1] M_after;
  real v_min = min(v_ep);
  real v_max = max(v_ep);
  // The row the model will actually subtract, from the function the model
  // calls.  The inline computation above is a second implementation of the
  // same integral; the two agreeing is the check that neither is wrong.
  row_vector[Kq - 1] qbar_fn = traversal_qbar(v_ep, Kq, vlo, vhi, npan2,
                                              xg, wg, P0v);
  // Gv of dgsn_grid_tdata.stan must come back from traversal_gram_shift at
  // shift = P0v, and the shifted Gram is what the new gauge uses.
  matrix[Kq - 1, Kq - 1] Gv_check = traversal_gram_shift(Kq, vlo, vhi, P0v, 4, xg, wg);
  matrix[Kq - 1, Kq - 1] Gv_gauge = Gv;
  matrix[Kq - 1, Kq - 1] Gvq_new
    = traversal_gram_shift(Kq, vlo, vhi, P0v + qbar_fn, 4, xg, wg);
  {
    matrix[N, Kq - 1] Z1;
    matrix[N, Kq - 1] Z2;
    matrix[N, Kq - 1] Zp;
    for (i in 1:N) {
      Zp[i] = traversal_row(v_ep[i], Kq, vlo, vhi) - P0v;
      for (r in 1:2) {
        int npan = r == 1 ? npan1 : npan2;
        row_vector[Kq - 1] num = rep_row_vector(0, Kq - 1);
        real den = 0;
        for (q in 1:npan) {
          real a0 = v_ep[i] * (q - 1.0) / npan;
          real a1 = v_ep[i] * (q * 1.0) / npan;
          vector[J] wj = fresh_abscissae(a1, a0, xg);
          vector[J] ww = fresh_weights(a1, a0, wg);
          for (j in 1:J) {
            real cw = cosh(wj[j]);
            num += (ww[j] * cw) * (traversal_row(wj[j], Kq, vlo, vhi) - P0v);
            den += ww[j] * cw;
          }
        }
        // v_i = 0 gives 0/0; the limit is the anchored row at zero, which is
        // zero by construction.
        if (r == 1) Z1[i] = abs(den) < 1e-12 ? rep_row_vector(0, Kq - 1) : num / den;
        else        Z2[i] = abs(den) < 1e-12 ? rep_row_vector(0, Kq - 1) : num / den;
      }
    }
    qbar_1 = rep_row_vector(0, Kq - 1);
    qbar_2 = rep_row_vector(0, Kq - 1);
    qbar_pt = rep_row_vector(0, Kq - 1);
    for (i in 1:N) { qbar_1 += Z1[i]; qbar_2 += Z2[i]; qbar_pt += Zp[i]; }
    qbar_1 /= N; qbar_2 /= N; qbar_pt /= N;
    M_1 = crossprod(Z1) / N;
    M_2 = crossprod(Z2) / N;
    M_pt = crossprod(Zp) / N;
    {
      matrix[N, Kq - 1] Za = Z2 - rep_matrix(qbar_2, N);
      qbar_after = rep_row_vector(0, Kq - 1);
      for (i in 1:N) qbar_after += Za[i];
      qbar_after /= N;
      M_after = crossprod(Za) / N;
    }
  }
}
