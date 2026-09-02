#include "transmission_line_system.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  int n1 = 2, n2 = 2, n3 = 2, n4 = 2;
  double h1 = 0.25, h2 = 0.25, h3 = 0.25, h4 = 0.25;
  double c1 = 1.0, c2 = 2.0, c3 = 3.0, c4 = 4.0;
  double k1 = 5.0, k2 = 6.0, k3 = 7.0, k4 = 8.0;

  Eigen::MatrixXd M0 =
      createM0(n1, n2, n3, n4, h1, h2, h3, h4, c1, c2, c3, c4, k1, k2, k3, k4);
  Eigen::MatrixXd M1 =
      createM1(n1, n2, n3, n4, h1, h2, h3, h4, c1, c2, c3, c4);

  Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
  std::cout << "M0 (" << M0.rows() << "x" << M0.cols() << "):\n"
            << M0.format(fmt) << "\n\n";
  std::cout << "M1 (" << M1.rows() << "x" << M1.cols() << "):\n"
            << M1.format(fmt) << "\n";

  return 0;
}
