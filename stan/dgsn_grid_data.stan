  // The declared constants of Sections 2.3, 4 and 5.1, and the grid.
  int<lower=2> Kpsi;              // turning order: modes 1..Kpsi-1
  int<lower=2> Kq;                // traversal order: modes 1..Kq-1
  int<lower=1> m;                 // panels per side of zero
  int<lower=1> J;                 // Gauss-Legendre points per panel
  vector[J] xg;                   // the J-point rule on [-1, 1]
  vector[J] wg;
  real<lower=0> delta;            // asinh scale in arclength, Eq. (umap)
  real<lower=0> delta_p;          // asinh scale in phase, Eq. (vmap)
  real pmin;                      // modelling window, days
  real pmax;
  real smin;                      // arclength interval, mag
  real smax;
  real<upper=0> gsmin;            // gauge span, mag: where turning norms are taken
  real<lower=0> gsmax;
  real<lower=0> eps;              // phase clamp width, fraction of range
  real<lower=0> eps_t;            // taper turn-off, fraction of range
  real<lower=0> eps_s;            // arclength clamp width, fraction of range
  real zp;                        // the archive zero point
