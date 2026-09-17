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
  // interval, by a two-sided composite rule fine enough (m panels a side,
  // J points) that products of the highest modes are integrated to rounding.
  // The gauge of Section 2.3 is c' G_s c = 1, d' G_s c = 0.
  matrix turning_gram_s(int Kpsi, real ulo, real uhi, real delta, row_vector P0, matrix T,
                        int m, vector xg, vector wg) {
    int n = Kpsi - 1;
    matrix[n, n] G = rep_matrix(0, n, n);
    vector[2 * m + 1] edges;
    for (i in 0:m) edges[i + 1] = ulo * (m - i) / m;
    for (i in 0:m) edges[m + 1 + i] = uhi * i / m;
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
      s[i] = exp(x_s) * (C[q] + part);
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
  // log of the predicted, tapered flux of Eq. (predflux) for N epochs.
  // band: 1 = g, 2 = r.  Returns the clamped phase and arclength too, so
  // that the diagnostics of Section 5.2 can be read without recomputation:
  // column 1 log flux, 2 log taper, 3 s before the arclength clamp, 4 after.
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
        out[i, 1] = -0.4 * log(10) * (mag - zp) + out[i, 2];
      }
    }
    return out;
  }

  // Eq. (like): -2 ln L with the tapered flux in residual and variance.
  real neg2_loglike(vector f, vector sig_f, array[] int band, vector f_hat,
                    real sigma_g, real sigma_r) {
    int N = rows(f);
    real total = 0;
    real k = 0.4 * log(10);
    for (i in 1:N) {
      real sb = (band[i] == 1) ? sigma_g : sigma_r;
      real var_i = square(sig_f[i]) + square(k * sb * f_hat[i]);
      total += square(f[i] - f_hat[i]) / var_i + log(var_i);
    }
    return total;
  }
}
