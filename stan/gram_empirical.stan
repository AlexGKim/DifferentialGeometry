// The turning and traversal Gram matrices under the measure the DATA impose,
// beside the ones the gauge imposes.
//
// Why this exists.  Eq. (gauge) normalises the turning pair under int f g ds
// over the gauge span and the traversal pair under int f g dv over the
// modelling window: both flat measures, weighting every part of the interval
// alike.  The epochs do not.  If the two measures disagree, a basis that is
// orthonormal under the gauge is strongly oblique under the data, and the
// coefficients the sampler moves are correlated no matter how the chain is
// run.  That is a statement about the basis, testable before any refit, and
// it is tested by comparing G_gauge with G_data.
//
// The comparison that matters is not G_data's own condition number -- that
// confounds the measure with the basis -- but the condition number of
// G_gauge^{-1/2} G_data G_gauge^{-1/2}, the data's inner product written in
// coordinates the gauge has already made orthonormal.  One means the two
// measures agree and the gauge is the right one; large means they do not.
// The whitening is linear algebra on matrices Stan produced, done outside.
//
// fixed_param only.  Takes any data dict the hierarchy programs take, plus
// the arclengths and phase coordinates the run's own epochs sit at.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  int<lower=1> N_s;               // epochs contributing to the empirical measure
  vector[N_s] s_obs;              // their arclengths, as the fit reports s_epoch
  vector[N_s] p_obs;              // their rest-frame phases
}
transformed data {
#include dgsn_grid_tdata.stan
}
generated quantities {
  matrix[Kpsi - 1, Kpsi - 1] gram_s = Gs;        // the gauge's, int f g ds
  matrix[Kq - 1, Kq - 1]     gram_v = Gv;        // the gauge's, int f g dv
  matrix[Kpsi - 1, Kpsi - 1] gram_s_data;        // the epochs', one weight each
  matrix[Kq - 1, Kq - 1]     gram_v_data;
  // The mean of each basis function under the two measures.  The gauge
  // anchors the traversal pair at v = 0, which fixes q(0) and leaves the
  // MEAN of q free; a free mean trades against the overall traversal rate,
  // so these are the numbers that say by how much.
  row_vector[Kq - 1]   qbar_data;                // mean basis row at the epochs
  row_vector[Kq - 1]   qbar_flat;                // and under flat dv
  row_vector[Kpsi - 1] psibar_data;
  // The turning basis's derivative row at the anchor, so that Psi'(0) for any
  // coefficient vector is a dot product outside.  Psi_0'(0) is kappa0 * delta;
  // Psi_1'(0) is the quantity a gauge condition Psi_1'(0) = 1 would fix, and
  // whether it can reach zero decides whether that condition is usable.
  row_vector[Kpsi - 1] dturn0 = turning_deriv_row(0, Kpsi, ulo, uhi, T);
  {
    matrix[Kpsi - 1, Kpsi - 1] As = rep_matrix(0, Kpsi - 1, Kpsi - 1);
    matrix[Kq - 1, Kq - 1]     Av = rep_matrix(0, Kq - 1, Kq - 1);
    for (i in 1:N_s) {
      row_vector[Kpsi - 1] r = (legendre_row(u_of_s(s_obs[i], delta), Kpsi, ulo, uhi)[2:Kpsi]
                                - P0) * T;
      As += r' * r;
      {
        row_vector[Kq] rv = legendre_row(v_of_p(p_obs[i], delta_p), Kq, vlo, vhi);
        row_vector[Kq - 1] rq = rv[2:Kq] - P0v;
        Av += rq' * rq;
      }
    }
    // Per epoch, so the scale is a mean and not an interval length; the
    // whitened condition number below is invariant to that scale anyway.
    gram_s_data = As / N_s;
    gram_v_data = Av / N_s;
    {
      row_vector[Kq - 1] sv = rep_row_vector(0, Kq - 1);
      row_vector[Kpsi - 1] ss = rep_row_vector(0, Kpsi - 1);
      for (i in 1:N_s) {
        sv += legendre_row(v_of_p(p_obs[i], delta_p), Kq, vlo, vhi)[2:Kq] - P0v;
        ss += (legendre_row(u_of_s(s_obs[i], delta), Kpsi, ulo, uhi)[2:Kpsi] - P0) * T;
      }
      qbar_data = sv / N_s;
      psibar_data = ss / N_s;
      // int (P_k - P_k(0)) dv / (vhi - vlo) = -P_k(0) for k >= 2, the
      // unanchored P_k integrating to zero over their own interval.
      qbar_flat = -P0v;
    }
  }
}
