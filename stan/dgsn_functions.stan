// The model of docs/tex/model-definition.tex, once.  #included by every program
// under stan/.  Nothing per-supernova is hard-coded here: every constant of
// Sections 2.3, 4 and 5.1 arrives as data.
//
// Conventions.  A "grid" is the two-sided composite Gauss-Legendre rule of
// Section 5.1: m equal panels on [lo, 0] and m on [0, hi], J nodes each, so
// that zero is edge m+1 of 2m+1 and the prefix table gives an exact zero
// there.  Indices are 1-based as in Stan: panel q spans edges[q]..edges[q+1].
functions {

  // ---- Legendre --------------------------------------------------------
  // Degrees 0..n-1 on [lo, hi], each scaled by sqrt((2k+1)/(hi-lo)) so that
  // the set is orthonormal under int f g dx over [lo, hi].
  row_vector legendre_row(real x, int n, real lo, real hi) {
    row_vector[n] P;
    real t = 2 * (x - lo) / (hi - lo) - 1;
    P[1] = 1;
    if (n > 1) P[2] = t;
    for (k in 2:(n - 1))
      P[k + 1] = ((2 * k - 1) * t * P[k] - (k - 1) * P[k - 1]) / k;
    for (k in 1:n)
      P[k] *= sqrt((2 * k - 1) / (hi - lo));
    return P;
  }

  // d/dx of legendre_row, by P'_{k+1} = P'_{k-1} + (2k+1) P_k in t.
  row_vector legendre_deriv_row(real x, int n, real lo, real hi) {
    row_vector[n] P;
    row_vector[n] D;
    real t = 2 * (x - lo) / (hi - lo) - 1;
    P[1] = 1;
    D[1] = 0;
    if (n > 1) {
      P[2] = t;
      D[2] = 1;
    }
    for (k in 2:(n - 1)) {
      P[k + 1] = ((2 * k - 1) * t * P[k] - (k - 1) * P[k - 1]) / k;
      D[k + 1] = D[k - 1] + (2 * k - 1) * P[k];
    }
    for (k in 1:n)
      D[k] *= sqrt((2 * k - 1) / (hi - lo)) * 2 / (hi - lo);
    return D;
  }

  // ---- The two bases of Section 2.3 -----------------------------------
  // Traversal: Legendre in v, degrees 1..Kq-1 (no constant).
  row_vector traversal_row(real v, int Kq, real vlo, real vhi) {
    return legendre_row(v, Kq, vlo, vhi)[2:Kq];
  }

  // Turning: Legendre in u, degrees 1..Kpsi-1, each shifted by its value at
  // s = 0 (u = 0) and re-orthonormalized by T, the transpose-inverse of the
  // Cholesky factor of the shifted Gram matrix.  T is upper triangular, so
  // mode k is a combination of degrees at most k: Gram-Schmidt in ascending
  // degree.
  row_vector turning_row(real u, int Kpsi, real ulo, real uhi,
                         row_vector P0, matrix T) {
    return (legendre_row(u, Kpsi, ulo, uhi)[2:Kpsi] - P0) * T;
  }

  row_vector turning_deriv_row(real u, int Kpsi, real ulo, real uhi, matrix T) {
    return legendre_deriv_row(u, Kpsi, ulo, uhi)[2:Kpsi] * T;
  }

  // The values P_k(0) the turning basis is shifted by.
  row_vector turning_shift(int Kpsi, real ulo, real uhi) {
    return legendre_row(0, Kpsi, ulo, uhi)[2:Kpsi];
  }

  // T for turning_row: Gram = I + (uhi - ulo) P0' P0, Gram = R R', T = inv(R)'.
  matrix turning_transform(int Kpsi, real ulo, real uhi, row_vector P0) {
    int n = Kpsi - 1;
    matrix[n, n] G = add_diag((uhi - ulo) * (P0' * P0), 1.0);
    matrix[n, n] R = cholesky_decompose(G);
    return mdivide_left_tri_low(R, identity_matrix(n))';
  }

  // The Gram matrix of the turning basis under ds, the measure the light
  // curve sees: G_s[j,k] = int Ptilde_j Ptilde_k delta cosh u du over the
  // gauge span [glo, ghi] in u, which lies inside the basis interval
  // [ulo, uhi], by a two-sided composite rule fine enough (m panels a side,
  // J points) that products of the highest modes are integrated to rounding.
  // The gauge of Section 2.3 is c' G_s c = 1, d' G_s c = 0.
  matrix turning_gram_s(int Kpsi, real ulo, real uhi, real glo, real ghi, real delta,
                        row_vector P0, matrix T, int m, vector xg, vector wg) {
    int n = Kpsi - 1;
    matrix[n, n] G = rep_matrix(0, n, n);
    vector[2 * m + 1] edges;
    for (i in 0:m) edges[i + 1] = glo * (m - i) / m;
    for (i in 0:m) edges[m + 1 + i] = ghi * i / m;
    for (q in 1:(2 * m)) {
      real half = 0.5 * (edges[q + 1] - edges[q]);
      real mid = 0.5 * (edges[q + 1] + edges[q]);
      for (j in 1:rows(xg)) {
        real u = mid + half * xg[j];
        row_vector[n] r = (legendre_row(u, Kpsi, ulo, uhi)[2:Kpsi] - P0) * T;
        G += (half * wg[j] * delta * cosh(u)) * (r' * r);
      }
    }
    return G;
  }

  // ---- Coordinates -------------------------------------------------------
  real u_of_s(real s, real delta) { return asinh(s / delta); }
  real s_of_u(real u, real delta) { return delta * sinh(u); }
  real v_of_p(real p, real delta_p) { return asinh(p / delta_p); }
  real p_of_v(real v, real delta_p) { return delta_p * sinh(v); }

  // ---- The two-sided composite Gauss-Legendre grid ----------------------
  // edges[1] = lo, edges[m+1] = 0, edges[2m+1] = hi.
  vector grid_edges(real lo, real hi, int m) {
    vector[2 * m + 1] e;
    for (i in 0:m) e[i + 1] = lo * (m - i) / m;
    for (i in 0:m) e[m + 1 + i] = hi * i / m;
    return e;
  }

  // Nodes and weights of the whole grid, panel-major: node (q-1)*J + j.
  vector grid_nodes(vector edges, vector xg) {
    int Q = rows(edges) - 1;
    int J = rows(xg);
    vector[Q * J] x;
    for (q in 1:Q) {
      real half = 0.5 * (edges[q + 1] - edges[q]);
      real mid = 0.5 * (edges[q + 1] + edges[q]);
      x[((q - 1) * J + 1):(q * J)] = mid + half * xg;
    }
    return x;
  }

  vector grid_weights(vector edges, vector wg) {
    int Q = rows(edges) - 1;
    int J = rows(wg);
    vector[Q * J] w;
    for (q in 1:Q)
      w[((q - 1) * J + 1):(q * J)] = 0.5 * (edges[q + 1] - edges[q]) * wg;
    return w;
  }

  // Which panel x lies in, for x in [lo, hi]: the last edge not above x.
  // A comparison loop rather than floor(), which Stan will not cast to an
  // integer when its argument depends on parameters.
  int panel_of(real x, vector edges) {
    int Q = rows(edges) - 1;
    int q = 1;
    for (i in 2:Q)
      if (x >= edges[i]) q = i;
    return q;
  }

  // Integral from zero to every edge, from the weighted node values fw.
  // C[m+1] is exactly zero; C[q] is the integral from 0 to edges[q].
  vector prefix_from_zero(vector fw, int m, int J) {
    int Q = 2 * m;
    vector[Q + 1] C;
    C[m + 1] = 0;
    for (q in (m + 1):Q)
      C[q + 1] = C[q] + sum(fw[((q - 1) * J + 1):(q * J)]);
    for (i in 1:m) {              // Stan ranges only ascend: walk q = m..1
      int q = m + 1 - i;
      C[q] = C[q + 1] - sum(fw[((q - 1) * J + 1):(q * J)]);
    }
    return C;
  }

  // The fresh J-point rule on [e0, x] for the remainder past the last edge.
  vector fresh_abscissae(real x, real e0, vector xg) {
    return 0.5 * (x + e0) + 0.5 * (x - e0) * xg;
  }
  vector fresh_weights(real x, real e0, vector wg) {
    return 0.5 * (x - e0) * wg;
  }

  // ---- Softening, Section 4 --------------------------------------------
  // eps * log(1 + exp(x/eps)) written so that neither branch overflows.
  real softplus(real x, real eps) {
    real z = x / eps;
    return eps * (fmax(z, 0) + log1p_exp(-abs(z)));
  }

  // Eq. (soft): the C-infinity clamp of w onto [0, 1].
  real clamp01(real w, real eps) {
    return softplus(w, eps) - softplus(w - 1, eps);
  }

  // log of the flux taper of Eq. (soft): softplus at eps_t/10, turn-off eps_t.
  real log_taper(real w, real eps_t) {
    real e = softplus(-w, eps_t / 10) + softplus(w - 1, eps_t / 10);
    return -square(e / eps_t);
  }

  // ---- The traversal, Eq. (travlog) -------------------------------------
  // s(p) for phases p inside the modelling window.  qc = b + x_p a.
  // Xv_nodes holds the traversal basis at the fixed v nodes; jacw_v the
  // weights times delta_p cosh v there.
  vector arclength(vector p, real x_s, vector qc,
                   int Kq, real vlo, real vhi, real delta_p, int m, int J,
                   vector edges_v, matrix Xv_nodes, vector jacw_v,
                   vector xg, vector wg) {
    int N = rows(p);
    vector[N] s;
    // Eq. (gauge) anchors the traversal pair at the origin, q_0(0) = q_1(0) = 0.
    // The anchored functions are b.P(v) - b.P(0) and a.P(v) - a.P(0), so both
    // constants ride with qc and ONE subtraction of qc.P(0) from the log rate
    // anchors the pair at once.  Taken out here rather than inside the
    // integrand so that it scales the prefix table and the fresh panel
    // identically; x_s is then ln(ds/dp) at p = 0 exactly.
    real x_s_a = x_s - legendre_row(0, Kq, vlo, vhi)[2:Kq] * qc;
    vector[2 * m + 1] C = prefix_from_zero(exp(Xv_nodes * qc) .* jacw_v, m, J);
    for (i in 1:N) {
      real v = v_of_p(p[i], delta_p);
      int q = panel_of(v, edges_v);
      real e0 = edges_v[q];
      vector[J] vj = fresh_abscissae(v, e0, xg);
      vector[J] wj = fresh_weights(v, e0, wg);
      real part = 0;
      for (j in 1:J)
        part += wj[j] * exp(traversal_row(vj[j], Kq, vlo, vhi) * qc)
                * delta_p * cosh(vj[j]);
      s[i] = exp(x_s_a) * (C[q] + part);
    }
    return s;
  }

  // ---- The curve, Eq. (quadrature) --------------------------------------
  // psi at u.  pc = d + x_psi c.
  real psi_at(real u, vector pc, int Kpsi, real ulo, real uhi,
              row_vector P0, matrix T) {
    return -pi() / 2 + turning_row(u, Kpsi, ulo, uhi, P0, T) * pc;
  }

  // gamma(s) as an N x 2 matrix, for arclengths s inside the interval.
  // Xu_nodes holds the turning basis at the fixed u nodes; jacw_u the
  // weights times delta cosh u there.
  matrix curve(vector s, vector pc,
               int Kpsi, real ulo, real uhi, real delta, int m, int J,
               vector edges_u, matrix Xu_nodes, vector jacw_u,
               vector xg, vector wg, row_vector P0, matrix T) {
    int N = rows(s);
    matrix[N, 2] g;
    vector[rows(Xu_nodes)] psi = -pi() / 2 + Xu_nodes * pc;
    vector[2 * m + 1] Cx = prefix_from_zero(cos(psi) .* jacw_u, m, J);
    vector[2 * m + 1] Cy = prefix_from_zero(sin(psi) .* jacw_u, m, J);
    for (i in 1:N) {
      real u = u_of_s(s[i], delta);
      int q = panel_of(u, edges_u);
      real e0 = edges_u[q];
      vector[J] uj = fresh_abscissae(u, e0, xg);
      vector[J] wj = fresh_weights(u, e0, wg);
      real px = 0;
      real py = 0;
      for (j in 1:J) {
        real ps = psi_at(uj[j], pc, Kpsi, ulo, uhi, P0, T);
        real wjac = wj[j] * delta * cosh(uj[j]);
        px += wjac * cos(ps);
        py += wjac * sin(ps);
      }
      g[i, 1] = Cx[q] + px;
      g[i, 2] = Cy[q] + py;
    }
    return g;
  }

  // ---- Observables, Section 4 ------------------------------------------
  // log of the predicted flux of Eq. (predflux) for N epochs.
  // band: 1 = g, 2 = r.  Returns the clamped phase and arclength too, so
  // that the diagnostics of Section 5.2 can be read without recomputation:
  // column 1 log flux, 2 the log likelihood weight, 3 s before the
  // arclength clamp, 4 after.
  //
  // Column 1 is the model's prediction and nothing else.  The taper used to
  // multiply it, which is not a deweighting: with f_hat driven to zero an
  // epoch leaving the window contributes its whole measured flux as a
  // residual, parameter-independent and gradientless, so the construction
  // was a barrier rather than the soft window it was wanted for, and in the
  // partial-taper zone it instructed the model to predict a fraction of the
  // true flux.  Column 2 is now a weight on that epoch's term in Eq. (like)
  // and the prediction is left alone.  Only the weight changes outside the
  // window, never the expected value.
  matrix log_flux(vector dt, array[] int band, real z,
                  real m_peak, real colr, real x_s, real x_p, real x_psi,
                  real delta_t,
                  vector a, vector b, vector c, vector d,
                  real pmin, real pmax, real smin, real smax,
                  real eps, real eps_t, real eps_s, real zp,
                  int Kq, real vlo, real vhi, real delta_p,
                  int Kpsi, real ulo, real uhi, real delta,
                  int m, int J, vector edges_v, matrix Xv_nodes, vector jacw_v,
                  vector edges_u, matrix Xu_nodes, vector jacw_u,
                  vector xg, vector wg, row_vector P0, matrix T) {
    int N = rows(dt);
    matrix[N, 4] out;
    vector[N] p_hat;
    vector[N] s_hat;
    for (i in 1:N) {
      real w = ((dt[i] - delta_t) / (1 + z) - pmin) / (pmax - pmin);
      p_hat[i] = pmin + (pmax - pmin) * clamp01(w, eps);
      out[i, 2] = log_taper(w, eps_t);
    }
    {
      vector[N] s = arclength(p_hat, x_s, b + x_p * a, Kq, vlo, vhi, delta_p,
                              m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      for (i in 1:N) {
        real vs = (s[i] - smin) / (smax - smin);
        s_hat[i] = smin + (smax - smin) * clamp01(vs, eps_s);
        out[i, 3] = s[i];
        out[i, 4] = s_hat[i];
      }
    }
    {
      matrix[N, 2] g = curve(s_hat, d + x_psi * c, Kpsi, ulo, uhi, delta, m, J,
                             edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
      for (i in 1:N) {
        real mag = (band[i] == 1) ? g[i, 1] + m_peak : g[i, 2] + m_peak - colr;
        out[i, 1] = -0.4 * log(10) * (mag - zp);
      }
    }
    return out;
  }

  // Eq. (like): -2 ln L in native flux, on the archive's errors alone.
  //
  // The variance used to carry a second term, (0.4 ln10 sigma_b f_hat)^2,
  // a fractional scatter with one fitted scalar per band.  It is gone, and
  // sigma_g and sigma_r with it: ZTF's flux_err already carries both the
  // sky and the source contribution, and the published errors are
  // calibrated -- over the baseline epochs the flux scatters at 0.977 times
  // flux_err, so DR2's own e_scale is already folded in.
  //
  // The log-determinant is gone for the same reason.  It was not optional
  // while the variance depended on the parameters: without it the density
  // falls monotonically in the variance and the fit drives sigma_b to
  // infinity.  With sig_f alone the variance is data, so sum log sig_f^2 is
  // an additive constant, and what this returns is -2 ln L up to that fixed
  // number.  Anything comparing an lp_sn here with one from a program that
  // kept the term has to restore it.
  //
  // w is the per-epoch weight of Eq. (soft), column 2 of log_flux
  // exponentiated: one inside the modelling window, falling smoothly to zero
  // outside it, so that an epoch the model does not claim to describe stops
  // contributing instead of being asked to predict no light.  The same
  // caveat as above applies to it -- an lp_sn from a program that tapered the
  // flux instead is not comparable with one from here.  On the selected
  // epochs every weight is one: the window exceeds the data window by five
  // days at each end and mu_delta is boxed at three, so no draw can carry a
  // phase within two days of the edge.
  real neg2_loglike(vector f, vector sig_f, vector f_hat, vector w) {
    return sum(w .* square((f - f_hat) ./ sig_f));
  }

  // Eq. (like) again, with Eq. (modelerr)'s fractional model error folded in.
  //
  //   V = sig_f^2 + (0.4 ln10 * sigma_b * f)^2
  //
  // sig_f already carries the sky and the source noise, so this second term
  // is not measurement error and would double-count it if it were.  It is
  // *model* error: the gap between the model curve and this supernova's own
  // curve, whose natural size is fractional because a curve good to a few
  // percent is good to a few percent bright or faint.  Without it the fit is
  // driven by whichever objects happen to be best observed -- over the sixty
  // of runs/stage3-k88-n60, Spearman corr(chi2 per epoch, peak S/N) is
  // +0.862 and corr with median sig_f/f is -0.797.
  //
  // The floor rides on the OBSERVED f, not on f_hat.  Where it dominates,
  // the flux is known to one or two percent, so (a f)^2 and (a f_hat)^2
  // agree to a few percent -- negligible in a term whose own size is being
  // estimated.  And at every epoch where the two could differ materially,
  // sig_f/f > 0.1, 426 of 3681, the floor does not dominate at all.  What
  // this buys is that V is data but for one scalar per band, so log V
  // couples to no geometric coordinate and the conditioning of the
  // coefficient block is untouched.
  //
  // log V is NOT optional once sigma_b is fitted.  Without it the density
  // rises monotonically in V and the fit drives sigma_b to infinity.  This
  // is required whichever flux the floor rides on, because sigma_b is a
  // parameter either way; riding on f only makes the term cheap, not
  // unnecessary.
  //
  // sfrac is per-epoch 0.4 ln10 sigma_b[band], so that one call covers both
  // bands.  What this returns is -2 ln L with no constant dropped: unlike
  // neg2_loglike above it is NOT comparable with an lp_sn from a program
  // that omits the log-determinant.
  real neg2_loglike_floor(vector f, vector sig_f, vector f_hat, vector w,
                          vector sfrac) {
    vector[num_elements(f)] V = square(sig_f) + square(sfrac .* f);
    return sum(w .* (square(f - f_hat) ./ V + log(V)));
  }

  // ---- The hierarchy, Section 5: a slice of supernovae for reduce_sum ---
  // -1/2 of Eq. (like) summed over a slice of supernovae.
  real partial_sum(array[] int sn, int first, int last,
                   vector dt, vector f, vector sig_f, array[] int band,
                   array[] int start, array[] int stop, vector z,
                   vector m_peak, vector colr, vector x_s, vector x_p, vector x_psi,
                   vector delta_t,
                   vector a, vector b, vector c, vector d,
                   real pmin, real pmax, real smin, real smax,
                   real eps, real eps_t, real eps_s, real zp,
                   int Kq, real vlo, real vhi, real delta_p,
                   int Kpsi, real ulo, real uhi, real delta,
                   int m, int J, vector edges_v, matrix Xv_nodes, vector jacw_v,
                   vector edges_u, matrix Xu_nodes, vector jacw_u,
                   vector xg, vector wg, row_vector P0, matrix T) {
    real total = 0;
    for (n in sn) {
      int lo = start[n];
      int hi = stop[n];
      matrix[hi - lo + 1, 4] lf = log_flux(dt[lo:hi], band[lo:hi], z[n],
          m_peak[n], colr[n], x_s[n], x_p[n], x_psi[n], delta_t[n],
          a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
          Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
          edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
      total += neg2_loglike(f[lo:hi], sig_f[lo:hi], exp(lf[:, 1]), exp(lf[:, 2]));
    }
    return -0.5 * total;
  }

  // The same slice with Eq. (modelerr)'s floor, for the programs that fit
  // sigma_b.  It is a second function rather than an argument to the first
  // because reduce_sum pins the signature, and because passing sfrac = 0
  // here would not reproduce partial_sum: the log-determinant would survive
  // as sum log sig_f^2, shifting lp__ by a constant and breaking every
  // comparison with a run made before this existed.
  real partial_sum_floor(array[] int sn, int first, int last,
                   vector dt, vector f, vector sig_f, vector sfrac,
                   array[] int band,
                   array[] int start, array[] int stop, vector z,
                   vector m_peak, vector colr, vector x_s, vector x_p, vector x_psi,
                   vector delta_t,
                   vector a, vector b, vector c, vector d,
                   real pmin, real pmax, real smin, real smax,
                   real eps, real eps_t, real eps_s, real zp,
                   int Kq, real vlo, real vhi, real delta_p,
                   int Kpsi, real ulo, real uhi, real delta,
                   int m, int J, vector edges_v, matrix Xv_nodes, vector jacw_v,
                   vector edges_u, matrix Xu_nodes, vector jacw_u,
                   vector xg, vector wg, row_vector P0, matrix T) {
    real total = 0;
    for (n in sn) {
      int lo = start[n];
      int hi = stop[n];
      matrix[hi - lo + 1, 4] lf = log_flux(dt[lo:hi], band[lo:hi], z[n],
          m_peak[n], colr[n], x_s[n], x_p[n], x_psi[n], delta_t[n],
          a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
          Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
          edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
      total += neg2_loglike_floor(f[lo:hi], sig_f[lo:hi], exp(lf[:, 1]),
                                  exp(lf[:, 2]), sfrac[lo:hi]);
    }
    return -0.5 * total;
  }

  // ---- Eq. (basispsi)'s constant modes ---------------------------------
  // Added 2026-11-23.  Nothing above this line changes, so every program on
  // record generates the same code and reproduces its runs; these names are
  // used only by hier_rot.stan.

  // The turning Gram EXTENDED by the constant mode: index 1 is Ptilde_0 = 1
  // and 2..Kpsi are the Ptilde_k of turning_row.  Eq. (gauge) is imposed on
  // the extended vectors once Psi_0 and Psi_1 carry d_0 and c_0;
  // turning_gram_s is this matrix's lower-right block.
  matrix turning_gram_ext(int Kpsi, real ulo, real uhi, real glo, real ghi,
                          real delta, row_vector P0, matrix T,
                          int m, vector xg, vector wg) {
    matrix[Kpsi, Kpsi] G = rep_matrix(0, Kpsi, Kpsi);
    vector[2 * m + 1] edges;
    for (i in 0:m) edges[i + 1] = glo * (m - i) / m;
    for (i in 0:m) edges[m + 1 + i] = ghi * i / m;
    for (q in 1:(2 * m)) {
      real half = 0.5 * (edges[q + 1] - edges[q]);
      real mid = 0.5 * (edges[q + 1] + edges[q]);
      for (j in 1:rows(xg)) {
        real u = mid + half * xg[j];
        row_vector[Kpsi] re
          = append_col(1.0, (legendre_row(u, Kpsi, ulo, uhi)[2:Kpsi] - P0) * T);
        G += (half * wg[j] * delta * cosh(u)) * (re' * re);
      }
    }
    return G;
  }

  // log_flux with the per-supernova anchor rotation of Eq. (basispsi).  A
  // constant alpha added to psi rotates gamma about the anchor exactly --
  // gamma(s) = int_0^s (cos psi, sin psi) ds' and gamma(0) = 0, so the
  // rotation comes outside the integral with no translation remainder --
  // and that is cheaper and more exact than re-running the quadrature.
  // alpha = 0 reproduces log_flux to the bit.
  matrix log_flux_rot(vector dt, array[] int band, real z,
                      real m_peak, real colr, real x_s, real x_p, real x_psi,
                      real delta_t, real alpha,
                      vector a, vector b, vector c, vector d,
                      real pmin, real pmax, real smin, real smax,
                      real eps, real eps_t, real eps_s, real zp,
                      int Kq, real vlo, real vhi, real delta_p,
                      int Kpsi, real ulo, real uhi, real delta,
                      int m, int J, vector edges_v, matrix Xv_nodes, vector jacw_v,
                      vector edges_u, matrix Xu_nodes, vector jacw_u,
                      vector xg, vector wg, row_vector P0, matrix T,
                      row_vector anchor) {
    int N = rows(dt);
    matrix[N, 4] out;
    vector[N] p_hat;
    vector[N] s_hat;
    for (i in 1:N) {
      real w = ((dt[i] - delta_t) / (1 + z) - pmin) / (pmax - pmin);
      p_hat[i] = pmin + (pmax - pmin) * clamp01(w, eps);
      out[i, 2] = log_taper(w, eps_t);
    }
    {
      vector[N] s = arclength_a(p_hat, x_s, b + x_p * a, anchor,
                                Kq, vlo, vhi, delta_p,
                                m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
      for (i in 1:N) {
        real vs = (s[i] - smin) / (smax - smin);
        s_hat[i] = smin + (smax - smin) * clamp01(vs, eps_s);
        out[i, 3] = s[i];
        out[i, 4] = s_hat[i];
      }
    }
    {
      matrix[N, 2] g = curve(s_hat, d + x_psi * c, Kpsi, ulo, uhi, delta, m, J,
                             edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
      real ca = cos(alpha);
      real sa = sin(alpha);
      for (i in 1:N) {
        real g0 = ca * g[i, 1] - sa * g[i, 2];
        real g1 = sa * g[i, 1] + ca * g[i, 2];
        real mag = (band[i] == 1) ? g0 + m_peak : g1 + m_peak - colr;
        out[i, 1] = -0.4 * log(10) * (mag - zp);
      }
    }
    return out;
  }

  // partial_sum_floor with the rotation carried per supernova.
  real partial_sum_floor_rot(array[] int sn, int first, int last,
                   vector dt, vector f, vector sig_f, vector sfrac,
                   array[] int band,
                   array[] int start, array[] int stop, vector z,
                   vector m_peak, vector colr, vector x_s, vector x_p, vector x_psi,
                   vector delta_t, vector alpha,
                   vector a, vector b, vector c, vector d,
                   real pmin, real pmax, real smin, real smax,
                   real eps, real eps_t, real eps_s, real zp,
                   int Kq, real vlo, real vhi, real delta_p,
                   int Kpsi, real ulo, real uhi, real delta,
                   int m, int J, vector edges_v, matrix Xv_nodes, vector jacw_v,
                   vector edges_u, matrix Xu_nodes, vector jacw_u,
                   vector xg, vector wg, row_vector P0, matrix T,
                   row_vector anchor) {
    real total = 0;
    for (n in sn) {
      int lo = start[n];
      int hi = stop[n];
      matrix[hi - lo + 1, 4] lf = log_flux_rot(dt[lo:hi], band[lo:hi], z[n],
          m_peak[n], colr[n], x_s[n], x_p[n], x_psi[n], delta_t[n], alpha[n],
          a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
          Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
          edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u, xg, wg, P0, T,
          anchor);
      total += neg2_loglike_floor(f[lo:hi], sig_f[lo:hi], exp(lf[:, 1]),
                                  exp(lf[:, 2]), sfrac[lo:hi]);
    }
    return -0.5 * total;
  }

  // ---- Eq. (gauge) on the traversal pair under the data's measure --------
  // Added 2026-11-24.  The anchor q_0(0) = q_1(0) = 0 fixes the additive
  // constant of q at a point; this fixes it as a mean instead.  Why: shift
  // the common log rate by epsilon and the coefficients by dq and
  //
  //   dlog s_n(p_i) = epsilon + <dq.(P - P(0))>_i ,
  //
  // the average against that epoch's traversal integrand over [0, v_i].  So
  // the shift in x_s is invisible to the extent the constant lies in the
  // anchored basis's span under the pooled measure of those averages, and
  // measured on runs/rot-k812-n60's epochs it lies 95.6 percent inside it
  // (runs/trav-degen).  Subtracting the mean puts the constant exactly
  // orthogonal to the span, so the residual is one rather than 0.0438.
  //
  // The subtracted row is data, not a parameter: the integrand average is
  // taken at the flat-rate reference qtilde = 0 and the phases at
  // delta_t = 0.  An error Delta in it leaves R^2 = Delta' M^-1 Delta,
  // quadratic, so the reference costs little.
  //
  // This is a gauge change and not a change of family.  exp(x_s) already
  // spans the constant direction, so the representable traversals are the
  // same set before and after; what moves is the meaning of x_s, from the
  // log rate at the anchor to the log rate at the measure's centre, and with
  // it what the hierarchy on x_s is a prior over.  One consequence worth
  // naming: qc = b + x_p a, so the shift qbar.qc is per supernova, equal to
  // qbar.b + x_p (qbar.a).  The fix therefore also sets the x_s--x_p
  // coupling that hier_kappa.stan fits, at kappa = -qbar.a, from the basis.
  row_vector traversal_qbar(vector v_ep, int Kq, real vlo, real vhi,
                            int npan, vector xg, vector wg,
                            row_vector P0v) {
    int n = Kq - 1;
    int N = rows(v_ep);
    int J = rows(xg);
    row_vector[n] acc = rep_row_vector(0, n);
    for (i in 1:N) {
      row_vector[n] num = rep_row_vector(0, n);
      real den = 0;
      for (q in 1:npan) {
        real a0 = v_ep[i] * (q - 1.0) / npan;
        real a1 = v_ep[i] * (q * 1.0) / npan;
        vector[J] vj = fresh_abscissae(a1, a0, xg);
        vector[J] wj = fresh_weights(a1, a0, wg);
        for (j in 1:J) {
          real cw = cosh(vj[j]);
          num += (wj[j] * cw) * (traversal_row(vj[j], Kq, vlo, vhi) - P0v);
          den += wj[j] * cw;
        }
      }
      // v_i = 0 gives 0/0 and the limit is the anchored row at zero, which
      // is zero by construction.
      acc += abs(den) < 1e-12 ? rep_row_vector(0, n) : num / den;
    }
    return acc / N;
  }

  // The traversal Gram of Eq. (gauge) for the basis P_k(v) - shift_k, under
  // int f g dv over the modelling window.  shift = P0v reproduces Gv of
  // dgsn_grid_tdata.stan, which is the test that this is the same object.
  // r' r is degree 2(Kq - 1), so a J-point rule is exact for J >= Kq - 1/2
  // and npan is a margin rather than a necessity.
  matrix traversal_gram_shift(int Kq, real vlo, real vhi, row_vector shift,
                              int npan, vector xg, vector wg) {
    int n = Kq - 1;
    int J = rows(xg);
    matrix[n, n] G = rep_matrix(0, n, n);
    for (q in 1:npan) {
      real a0 = vlo + (vhi - vlo) * (q - 1.0) / npan;
      real a1 = vlo + (vhi - vlo) * (q * 1.0) / npan;
      vector[J] vj = fresh_abscissae(a1, a0, xg);
      vector[J] wj = fresh_weights(a1, a0, wg);
      for (j in 1:J) {
        row_vector[n] r = traversal_row(vj[j], Kq, vlo, vhi) - shift;
        G += wj[j] * (r' * r);
      }
    }
    return G;
  }

  // arclength with the subtracted row passed in rather than taken as P(0).
  // anchor = P0v is arclength() exactly; anchor = P0v + qbar is the gauge
  // above.  The constant rides out of the integrand into the prefactor
  // either way, so this costs nothing beyond the row it is handed.
  vector arclength_a(vector p, real x_s, vector qc, row_vector anchor,
                     int Kq, real vlo, real vhi, real delta_p, int m, int J,
                     vector edges_v, matrix Xv_nodes, vector jacw_v,
                     vector xg, vector wg) {
    int N = rows(p);
    vector[N] s;
    real x_s_a = x_s - anchor * qc;
    vector[2 * m + 1] C = prefix_from_zero(exp(Xv_nodes * qc) .* jacw_v, m, J);
    for (i in 1:N) {
      real v = v_of_p(p[i], delta_p);
      int q = panel_of(v, edges_v);
      real e0 = edges_v[q];
      vector[J] vj = fresh_abscissae(v, e0, xg);
      vector[J] wj = fresh_weights(v, e0, wg);
      real part = 0;
      for (j in 1:J)
        part += wj[j] * exp(traversal_row(vj[j], Kq, vlo, vhi) * qc)
                * delta_p * cosh(vj[j]);
      s[i] = exp(x_s_a) * (C[q] + part);
    }
    return s;
  }

  // ---- Eq. (gauge) as a comoving frame -----------------------------------
  // Added 2026-11-25.  Eq. (gauge) asks <q_0,q_1> = 0 and <Psi_0,Psi_1> = 0.
  // Through 2026-11-24 both were imposed by projection: sample b_til (d_til)
  // in the ambient space and subtract its component along the normalised
  // partner.  That leaves the sampled coordinates attached to the ambient
  // axes while the admissible space is the partner's complement, so rotating
  // a (or c) drags the whole space the other vector lives in.  Measured on
  // runs/qbar-gfix-k812-n60 draws 251+, |corr(c_k, d_k)| runs 0.69 to 0.90 at
  // every matched order k = 1..11, uniformly positive, which is that drag and
  // not a covariance accident of one draw.  It is also curvature that varies
  // with position -- the constraint surface is curved -- which is why the
  // four constant-metric and constant-linear-map fixes of 2026-11-23 could
  // not encode it.
  //
  // Here the frame moves instead.  With G = L L' and e normalised so that
  // e' G e = 1, ehat = L' e is Euclidean-unit, and the Householder
  // reflection H = I - 2 v v'/(v'v) with v = ehat + sign(ehat_1) e_1 is
  // symmetric, involutive and takes ehat to -sign(ehat_1) e_1.  So for any y
  // of length n - 1,
  //
  //     de = L^-T H [0; y]
  //
  // satisfies de' G e = [0;y]' H (L' e) = [0;y]' (-s e_1) = 0 exactly, and
  // de' G de = y' y, so y are the components of the free part in a
  // G-orthonormal frame carried by e.  The gauge holds by construction, with
  // no component to discard and no prior owed on one: the dimension drops by
  // one per pair and the four-dimensional null space of the flux Jacobian
  // loses its two projection directions.
  //
  // The likelihood is unchanged: at three matched draws all 362 gradient
  // components the two programs share are bit-identical.  The prior is not
  // quite unchanged, and the earlier wording here ("a constant") was wrong.
  // lp(rot_frame) - lp(hier_rot) is predicted to full double precision as
  //
  //     1.837877 + (t_b^2 + t_d^2)/2 - log(||Gvq a|| ||Gse ce||),
  //
  // t being the discarded component.  The t part is an independent standard
  // normal on a coordinate reaching neither the forward map nor the
  // likelihood, so it integrates out exactly and nothing is owed for it.
  // The remainder, log(||Gvq a|| ||Gse ce||), is the projection Jacobian and
  // does depend on where a and ce point, so the two priors on the admissible
  // space differ.  Measured by stan/frame_map at the three probe draws of
  // runs/frame-probe it is 0.3776, 0.3976, 0.3834 nats: a reweighting of the
  // posterior spanning 1.020 across the sampled region, against lp ~ -31400
  // and globals with posterior sd of 10 to 20 per cent.  Too small to move a
  // reported number, and not zero.
  //
  // sign(ehat_1) is taken on the stable branch, so v never cancels; ehat_1 is
  // c_0 and b's first coefficient after whitening and is not pinned near
  // zero by anything, but the branch costs nothing either way.  The sign flip
  // the branch introduces is absorbed into y, which is free.
  // ---- Eq. (gauge) as a point on a Grassmannian --------------------------
  // Added 2026-11-26.  The gauge-invariant content of a pair is the affine
  // line {f_0 + x f_1}: a 2-plane in the K-dimensional coefficient space,
  // plus the direction of f_1 within it and the perpendicular extent of f_0.
  // So sample the plane, an angle and a length, and both conditions of
  // Eq. (gauge) hold identically -- no normalisation of a sampled vector, no
  // projection, and no frame carried by a moving partner.
  //
  // Everything is done in the whitened coordinates xhat = L' x, where
  // G = L L' so that <x,y>_G = xhat . yhat: the plane's chart, its
  // orthonormalisation and the orientation are all Euclidean there, and the
  // pair is returned by L^-T at the end.  With Qhat the whitened plane's
  // Euclidean-orthonormal 2-frame and u the unit 2-vector fixing the
  // orientation, R(theta) u and R(theta + pi/2) u are orthonormal in the
  // plane, so <f_1,f_1>_G = 1 and <f_0,f_1>_G = 0 by construction and the
  // length h is the free, physical norm of f_0.
  //
  // mgs2 is modified Gram-Schmidt on two columns: smooth wherever the two
  // are independent, which the chart guarantees near its reference plane.
  // Measured on runs/qbar-gfix-k812-n60 chain 2, sixteen sampling draws, the
  // chart centred on one of them has ||X|| at most 0.162 for the turning
  // pair and 0.519 for the traversal, so the posterior sits well inside it.
  matrix mgs2(matrix M) {
    int n = rows(M);
    vector[n] q1 = M[, 1] / sqrt(dot_self(M[, 1]));
    vector[n] v2 = M[, 2] - dot_product(M[, 2], q1) * q1;
    return append_col(q1, v2 / sqrt(dot_self(v2)));
  }

  // log det(I_2 + X' X), the chart's Jacobian up to the power -n/2 that makes
  // the measure the O(n)-invariant one on the Grassmannian.  Flat in X is a
  // different and chart-dependent prior; this term is what makes "flat" mean
  // flat on the space of planes.
  real log_det_ipxtx(matrix X) {
    matrix[2, 2] S = crossprod(X);
    return log((1 + S[1, 1]) * (1 + S[2, 2]) - S[1, 2] * S[2, 1]);
  }

  vector frame_lift(matrix Lt, matrix Ltinv, vector e, vector y) {
    int n = rows(e);
    vector[n] ehat = Lt * e;
    vector[n] v = ehat;
    vector[n] full = append_row(0, y);
    v[1] += ehat[1] >= 0 ? 1 : -1;
    return Ltinv * (full - (2 * dot_product(v, full) / dot_self(v)) * v);
  }

  // The inverse, for mapping stored draws into the frame coordinates:
  // [0; y] = H (L' de), so y is that product's tail.  Exact to rounding, and
  // the test that frame_lift is the map it claims to be.
  vector frame_drop(matrix Lt, vector e, vector de) {
    int n = rows(e);
    vector[n] ehat = Lt * e;
    vector[n] v = ehat;
    vector[n] z = Lt * de;
    v[1] += ehat[1] >= 0 ? 1 : -1;
    return (z - (2 * dot_product(v, z) / dot_self(v)) * v)[2:n];
  }
}
