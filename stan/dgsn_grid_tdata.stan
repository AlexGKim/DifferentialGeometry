  // Everything fixed for the run: coordinates of the windows, edges, nodes,
  // weights times Jacobians, and the two bases at the nodes.
  real vlo = v_of_p(pmin, delta_p);
  real vhi = v_of_p(pmax, delta_p);
  real ulo = u_of_s(smin, delta);
  real uhi = u_of_s(smax, delta);
  vector[2 * m + 1] edges_v = grid_edges(vlo, vhi, m);
  vector[2 * m + 1] edges_u = grid_edges(ulo, uhi, m);
  vector[2 * m * J] v_nodes = grid_nodes(edges_v, xg);
  vector[2 * m * J] u_nodes = grid_nodes(edges_u, xg);
  vector[2 * m * J] w_v = grid_weights(edges_v, wg);
  vector[2 * m * J] w_u = grid_weights(edges_u, wg);
  vector[2 * m * J] jacw_v = w_v .* (delta_p * cosh(v_nodes));
  vector[2 * m * J] jacw_u = w_u .* (delta * cosh(u_nodes));
  row_vector[Kpsi - 1] P0 = turning_shift(Kpsi, ulo, uhi);
  matrix[Kpsi - 1, Kpsi - 1] T = turning_transform(Kpsi, ulo, uhi, P0);
  real gulo = u_of_s(gsmin, delta);
  real guhi = u_of_s(gsmax, delta);
  matrix[Kpsi - 1, Kpsi - 1] Gs = turning_gram_s(Kpsi, ulo, uhi, gulo, guhi, delta, P0, T,
                                                 200, xg, wg);
  matrix[2 * m * J, Kq - 1] Xv_nodes;
  matrix[2 * m * J, Kpsi - 1] Xu_nodes;
  for (n in 1:(2 * m * J)) {
    Xv_nodes[n] = traversal_row(v_nodes[n], Kq, vlo, vhi);
    Xu_nodes[n] = turning_row(u_nodes[n], Kpsi, ulo, uhi, P0, T);
  }
