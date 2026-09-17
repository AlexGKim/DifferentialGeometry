// The forward map at fixed parameters, run with fixed_param.  Stage 1 of
// Section 5.2: s(p), gamma(s), psi(s), the fluxes of given epochs, and the
// two bases at given points, which is what lets Python fit a SALT2 track by
// least squares without evaluating a basis itself.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
  vector[Kq - 1] a;               // gauge-fixed global coefficients
  vector[Kq - 1] b;
  vector[Kpsi - 1] c;
  vector[Kpsi - 1] d;
  real m_peak;                    // the six of Table 1
  real colr;
  real x_s;
  real x_p;
  real x_psi;
  real delta_t;                   // t_max - t0
  real z;
  int<lower=0> Np;                // phases at which s(p) is wanted
  vector[Np] p_eval;
  int<lower=0> Ns;                // arclengths at which gamma, psi are wanted
  vector[Ns] s_eval;
  int<lower=0> N;                 // epochs at which the flux is wanted
  vector[N] dt;
  array[N] int<lower=1, upper=2> band;
}
transformed data {
#include dgsn_grid_tdata.stan
}
generated quantities {
  vector[Np] s_of_p = arclength(p_eval, x_s, b + x_p * a, Kq, vlo, vhi, delta_p,
                                m, J, edges_v, Xv_nodes, jacw_v, xg, wg);
  matrix[Ns, 2] gamma = curve(s_eval, d + x_psi * c, Kpsi, ulo, uhi, delta, m, J,
                              edges_u, Xu_nodes, jacw_u, xg, wg, P0, T);
  vector[Ns] psi;
  matrix[Ns, Kpsi - 1] Xu;
  matrix[Np, Kq - 1] Xv;
  matrix[N, 4] lf = log_flux(dt, band, z, m_peak, colr, x_s, x_p, x_psi, delta_t,
                             a, b, c, d, pmin, pmax, smin, smax, eps, eps_t, eps_s, zp,
                             Kq, vlo, vhi, delta_p, Kpsi, ulo, uhi, delta, m, J,
                             edges_v, Xv_nodes, jacw_v, edges_u, Xu_nodes, jacw_u,
                             xg, wg, P0, T);
  vector[N] f_hat = exp(lf[:, 1]);
  real kappa0 = turning_deriv_row(0, Kpsi, ulo, uhi, T) * (d + x_psi * c) / delta;
  // The grid itself, for tests of orthonormality and resolution.
  vector[2 * m * J] u_nodes_out = u_nodes;
  vector[2 * m * J] v_nodes_out = v_nodes;
  vector[2 * m * J] w_u_out = w_u;
  vector[2 * m * J] w_v_out = w_v;
  matrix[2 * m * J, Kpsi - 1] Xu_nodes_out = Xu_nodes;
  matrix[2 * m * J, Kq - 1] Xv_nodes_out = Xv_nodes;
  matrix[Kpsi - 1, Kpsi - 1] T_out = T;
  matrix[Kpsi - 1, Kpsi - 1] Gs_out = Gs;
  for (i in 1:Ns) {
    real u = u_of_s(s_eval[i], delta);
    psi[i] = psi_at(u, d + x_psi * c, Kpsi, ulo, uhi, P0, T);
    Xu[i] = turning_row(u, Kpsi, ulo, uhi, P0, T);
  }
  for (i in 1:Np)
    Xv[i] = traversal_row(v_of_p(p_eval[i], delta_p), Kq, vlo, vhi);
}
