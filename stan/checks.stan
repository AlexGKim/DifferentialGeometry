// Tests only.  The quadrature, clamps and taper on closed-form integrands,
// through the same grid, prefix table and partial panel the model uses.
// No model program carries a test switch; this one is the switch.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  real k;                         // constant psi' for the circle
  int<lower=0> Ns;                // arclengths for the circle
  vector[Ns] s_eval;
  int<lower=0> Np;                // phases for the unit-rate traversal
  vector[Np] p_eval;
  int<lower=0> Nx;                // points for exp(0.3x)cos(3x) on [pmin, pmax] in v
  vector[Nx] x_eval;
  int<lower=0> Nw;                // arguments for the clamp and taper
  vector[Nw] w_eval;
}
transformed data {
#include dgsn_grid_tdata.stan
}
generated quantities {
  // Circle: psi = -pi/2 + k s, integrated in u with the Jacobian, so that
  // gamma(s) = ((1 - cos ks)/k, -sin(ks)/k).
  matrix[Ns, 2] circle;
  // Unit rate: ds/dp = 1, integrated in v with the Jacobian, so s = p.
  vector[Np] s_unit;
  // A smooth test integrand on the v grid: int_0^x exp(0.3t) cos(3t) dt.
  vector[Nx] I_test;
  vector[Nw] clamped;
  vector[Nw] taper;
  vector[Nw] sp;
  {
    vector[2 * m * J] ps = -pi() / 2 + k * delta * sinh(u_nodes);
    vector[2 * m + 1] Cx = prefix_from_zero(cos(ps) .* jacw_u, m, J);
    vector[2 * m + 1] Cy = prefix_from_zero(sin(ps) .* jacw_u, m, J);
    for (i in 1:Ns) {
      real u = u_of_s(s_eval[i], delta);
      int q = panel_of(u, edges_u);
      vector[J] uj = fresh_abscissae(u, edges_u[q], xg);
      vector[J] wj = fresh_weights(u, edges_u[q], wg);
      real px = 0;
      real py = 0;
      for (j in 1:J) {
        real p_ = -pi() / 2 + k * delta * sinh(uj[j]);
        px += wj[j] * delta * cosh(uj[j]) * cos(p_);
        py += wj[j] * delta * cosh(uj[j]) * sin(p_);
      }
      circle[i, 1] = Cx[q] + px;
      circle[i, 2] = Cy[q] + py;
    }
  }
  {
    vector[2 * m + 1] C = prefix_from_zero(jacw_v, m, J);
    for (i in 1:Np) {
      real v = v_of_p(p_eval[i], delta_p);
      int q = panel_of(v, edges_v);
      vector[J] vj = fresh_abscissae(v, edges_v[q], xg);
      vector[J] wj = fresh_weights(v, edges_v[q], wg);
      real part = 0;
      for (j in 1:J) part += wj[j] * delta_p * cosh(vj[j]);
      s_unit[i] = C[q] + part;
    }
  }
  {
    vector[2 * m + 1] C = prefix_from_zero(exp(0.3 * v_nodes) .* cos(3 * v_nodes) .* w_v, m, J);
    for (i in 1:Nx) {
      real x = x_eval[i];
      int q = panel_of(x, edges_v);
      vector[J] xj = fresh_abscissae(x, edges_v[q], xg);
      vector[J] wj = fresh_weights(x, edges_v[q], wg);
      real part = 0;
      for (j in 1:J) part += wj[j] * exp(0.3 * xj[j]) * cos(3 * xj[j]);
      I_test[i] = C[q] + part;
    }
  }
  for (i in 1:Nw) {
    clamped[i] = clamp01(w_eval[i], eps);
    taper[i] = exp(log_taper(w_eval[i], eps_t));
    sp[i] = softplus(w_eval[i], eps);
  }
}
