#include "transmission_line_solver.hpp"
#include "transmission_line_system.hpp"
#include <algorithm>
#include <complex>
#include <numeric>
#include <vector>

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::VectorXd;

Eigen::MatrixXd build_transmission_line_A(int n1, int n2, int n3, int n4,
                                           double h1, double h2, double h3,
                                           double h4, double c1, double c2,
                                           double c3, double c4, double k1,
                                           double k2, double k3, double k4) {
  MatrixXd M0 = createM0(n1, n2, n3, n4, h1, h2, h3, h4, c1, c2, c3, c4, k1,
                          k2, k3, k4);
  MatrixXd M1 = createM1(n1, n2, n3, n4, h1, h2, h3, h4, c1, c2, c3, c4);

  // Mathematica: final = Inverse[M0] . (-1 * M1)
  return -M0.colPivHouseholderQr().solve(M1);
}

Eigen::VectorXd evaluate_transmission_line_wave(
    int n1, int n2, int n3, int n4, double h1, double h2, double h3,
    double h4, double c1, double c2, double c3, double c4, double k1,
    double k2, double k3, double k4,
    const Eigen::VectorXd &initialConditions, double t, int drop) {
  // Mathematica: loadedStateVector[t_] := If[t==0, initialConditions, ...].
  // The truncated eigenbasis reconstruction only approximates the initial
  // condition (it drops `drop` modes), so at t=0 return it exactly instead.
  if (t == 0.0) {
    return initialConditions;
  }

  MatrixXd A = build_transmission_line_A(n1, n2, n3, n4, h1, h2, h3, h4, c1,
                                          c2, c3, c4, k1, k2, k3, k4);

  Eigen::EigenSolver<MatrixXd> solver(A);
  VectorXcd eigenvalues = solver.eigenvalues();
  MatrixXcd eigenvectors = solver.eigenvectors();

  // Sort descending by magnitude, matching Mathematica's Eigenvalues[] order.
  int n = static_cast<int>(eigenvalues.size());
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return std::abs(eigenvalues[a]) > std::abs(eigenvalues[b]);
  });

  VectorXcd sortedEigenvalues(n);
  MatrixXcd sortedEigenvectors(eigenvectors.rows(), n);
  for (int i = 0; i < n; ++i) {
    sortedEigenvalues[i] = eigenvalues[order[i]];
    sortedEigenvectors.col(i) = eigenvectors.col(order[i]);
  }

  VectorXcd fourierCoefficients = sortedEigenvectors.colPivHouseholderQr().solve(
      initialConditions.cast<std::complex<double>>());

  // Drop[..., -drop]: discard the trailing (smallest-magnitude) modes.
  int keep = n - drop;
  VectorXcd keptEigenvalues = sortedEigenvalues.head(keep);
  MatrixXcd keptEigenvectors = sortedEigenvectors.leftCols(keep);
  VectorXcd keptCoefficients = fourierCoefficients.head(keep);

  VectorXcd expTerms = (keptEigenvalues * t).array().exp().matrix();
  VectorXcd state = keptEigenvectors * expTerms.asDiagonal() * keptCoefficients;

  return state.real();
}
