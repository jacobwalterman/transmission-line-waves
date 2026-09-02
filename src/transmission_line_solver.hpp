#include <Eigen/Dense>

// Builds A = M0^-1 M1, the RHS matrix of dx/dt = A x. Exposed for debugging.
Eigen::MatrixXd build_transmission_line_A(int n1, int n2, int n3, int n4,
                                           double h1, double h2, double h3,
                                           double h4, double c1, double c2,
                                           double c3, double c4, double k1,
                                           double k2, double k3, double k4);

// Solves dx/dt = M0^-1 M1 x (x = [w1..w4, v1..v4] stacked per createM1's
// state vector) via eigendecomposition, and evaluates the real part of the
// resulting state at time t, starting from initialConditions at t=0.
//
// `drop` mirrors Mathematica's createStateVectorTerms: the controller*eq
// rows are identically zero on M1's side, so A has an exact `drop`-
// dimensional null space; the corresponding near-zero eigenvalues are
// numerical noise (not physical modes) and are dropped (smallest-magnitude
// first, after sorting descending by magnitude) before reconstructing the
// state -- otherwise their ill-conditioned eigenvectors dominate the
// solution with spurious growth.
Eigen::VectorXd evaluate_transmission_line_wave(
    int n1, int n2, int n3, int n4, double h1, double h2, double h3,
    double h4, double c1, double c2, double c3, double c4, double k1,
    double k2, double k3, double k4,
    const Eigen::VectorXd &initialConditions, double t, int drop = 4);
