#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::VectorXd;

MatrixXd createOrderReducedMatrix(int n, double k1, double k3) {
  double h = 1.0 / (n + 1);

  MatrixXd zero = MatrixXd::Zero(n + 1, n + 1);
  MatrixXd identity = MatrixXd::Identity(n + 1, n + 1);

  MatrixXd aMatrix = MatrixXd::Zero(n + 1, n + 1);
  aMatrix.diagonal(-1).setConstant(1.0 / (std::pow(h, 2)));
  aMatrix.diagonal().setConstant(-2.0 / (std::pow(h, 2)));
  aMatrix.diagonal(1).setConstant(1.0 / (std::pow(h, 2)));
  aMatrix(n, n) = -1.0 / (std::pow(h, 2));

  MatrixXd pMatrix = MatrixXd::Zero(n + 1, n + 1);
  pMatrix.diagonal(-1).setConstant(0.25);
  pMatrix.diagonal().setConstant(0.5);
  pMatrix.diagonal(1).setConstant(0.25);
  pMatrix(n, n) = 0.25;
  pMatrix = pMatrix.inverse();

  MatrixXd bottomLeftMatrix = pMatrix * aMatrix;

  MatrixXd bMatrix = MatrixXd::Zero(n + 1, n + 1);
  bMatrix(n, n) = -k3 / h;

  MatrixXd bottomRightMatrix = pMatrix * bMatrix - (k1 * pMatrix);

  MatrixXd orderReducedMatrix = MatrixXd::Zero(2 * (n + 1), 2 * (n + 1));
  orderReducedMatrix.block(0, 0, n + 1, n + 1) = zero;
  orderReducedMatrix.block(0, n + 1, n + 1, n + 1) = identity;
  orderReducedMatrix.block(n + 1, 0, n + 1, n + 1) = bottomLeftMatrix;
  orderReducedMatrix.block(n + 1, n + 1, n + 1, n + 1) = bottomRightMatrix;

  return orderReducedMatrix;
}

VectorXcd evaluate_solution(int n, double t) {
  double k1 = 1.0;
  double k3 = 1.0;

  VectorXd initialConditions = VectorXd::Random(2 * (n + 1));
  MatrixXd orderReducedMatrix = createOrderReducedMatrix(n, k1, k3);
  Eigen::EigenSolver<MatrixXd> solver(orderReducedMatrix);
  VectorXcd eigenvalues = solver.eigenvalues();
  MatrixXcd eigenvectors = solver.eigenvectors();
  VectorXcd fourierCoefficients = eigenvectors.colPivHouseholderQr().solve(
      initialConditions.cast<std::complex<double>>());

  VectorXcd expTerms = (eigenvalues * t).array().exp().matrix();

  VectorXcd positionsAndVelocities =
      eigenvectors * expTerms.asDiagonal() * fourierCoefficients;

  return positionsAndVelocities;
}
