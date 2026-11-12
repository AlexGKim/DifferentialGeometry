// The turning pair's gauge inner product, so that a start can be placed on
// the gauge surface without Python evaluating a basis.
//
// Why this exists.  Eq. (gauge) orthogonalises Psi_0 against Psi_1 under
// arclength, d -> d_til - (d_til' Gs c) c, so a start that carries d_til from
// a run without the mode lands at d_til minus its Gs-component along c, not
// at d_til.  Choosing c orthogonal to the carried d_til removes that, and the
// projection needs Gs.  Gs is a basis integral -- the Gram matrix of the
// shifted, re-orthonormalised Legendre functions under ds over the gauge span
// -- and the standing rule is that Python evaluates no basis and no
// quadrature.  So Stan computes it, in the one place it is already computed,
// and hands it over; the projection downstream is linear algebra on a matrix
// Stan produced.
//
// fixed_param only: there are no parameters, and transformed data has already
// done the work by the time generated quantities runs.  Takes any data dict
// the hierarchy programs take, CmdStan ignoring the names it does not declare.
#include dgsn_functions.stan
data {
#include dgsn_grid_data.stan
}
transformed data {
#include dgsn_grid_tdata.stan
}
generated quantities {
  // Gs, the gauge metric of Eq. (gauge) on the turning coefficients, and the
  // identity is the traversal pair's, the flat measure in v making that one
  // plain dot products.
  matrix[Kpsi - 1, Kpsi - 1] gram_s = Gs;
}
