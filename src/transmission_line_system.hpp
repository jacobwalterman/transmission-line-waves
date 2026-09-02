#pragma once

#include "system_matrix.hpp"
#include <Eigen/Dense>
#include <cmath>

// Named quantities used across the M0 (velocity/derivative) and M1
// (position/strain) system matrices. dw*/dv* are the M0-side unknowns,
// w*/v* are the M1-side unknowns; the two never mix within one matrix.
enum class FuncId {
  dw1, dw2, dw3, dw4,
  dv1, dv2, dv3, dv4,
  w1, w2, w3, w4,
  v1, v2, v3, v4,
  body1eq, body2eq, body3eq, body4eq,
  controller1eq, controller2eq, controller3eq, controller4eq,
  reduction1eq, reduction2eq, reduction3eq, reduction4eq,
  straineq,
};

using Eq = EquationSpec<FuncId>;
using Range = RangeSpec<FuncId>;
using Term = TermTemplate<FuncId>;

inline std::vector<Range> pieceSolutionVector(int n1, int n2, int n3, int n4) {
  return {
      {FuncId::body1eq, 0, n1},
      {FuncId::controller1eq, 1, 1},
      {FuncId::controller2eq, 1, 1},
      {FuncId::body2eq, 0, n2},
      {FuncId::body3eq, 0, n3},
      {FuncId::controller3eq, 1, 1},
      {FuncId::controller4eq, 1, 1},
      {FuncId::body4eq, 0, n4},
      {FuncId::reduction1eq, 0, n1},
      {FuncId::reduction2eq, 0, n2},
      {FuncId::reduction3eq, 0, n3},
      {FuncId::reduction4eq, 0, n4},
  };
}

inline Eigen::MatrixXd createM0(int n1, int n2, int n3, int n4, double h1,
                                 double h2, double h3, double h4, double c1,
                                 double c2, double c3, double c4, double k1,
                                 double k2, double k3, double k4) {
  std::vector<Range> stateVector = {
      {FuncId::dw1, 1, n1 + 1}, {FuncId::dw2, 0, n2},
      {FuncId::dw3, 0, n3 + 1}, {FuncId::dw4, 0, n4},
      {FuncId::dv1, 0, n1 + 1}, {FuncId::dv2, 0, n2},
      {FuncId::dv3, 0, n3 + 1}, {FuncId::dv4, 0, n4 + 1},
  };
  std::vector<Range> solutionVector = pieceSolutionVector(n1, n2, n3, n4);

  std::vector<Eq> equations = {
      {FuncId::dw1, 0, 0, {{FuncId::dw1, 1, 0.0}}},
      {FuncId::dw4, n4 + 1, n4 + 1, {{FuncId::dw4, -1, 0.0}}},
      {FuncId::controller1eq,
       1,
       1,
       {
           {FuncId::dw1, n1, k1},
           {FuncId::dw2, -1, k1},
           {FuncId::dv1, n1, c1 * c1},
           {FuncId::dv2, -1, -(c2 * c2)},
       }},
      {FuncId::controller2eq,
       1,
       1,
       {
           {FuncId::dw1, n1, k2},
           {FuncId::dw2, -1, -k2},
           {FuncId::dv1, n1, c1 * c1},
           {FuncId::dv2, -1, c2 * c2},
       }},
      {FuncId::controller3eq,
       1,
       1,
       {
           {FuncId::dw3, n3, k3},
           {FuncId::dw4, -1, k3},
           {FuncId::dv3, n3, c3 * c3},
           {FuncId::dv4, -1, -(c4 * c4)},
       }},
      {FuncId::controller4eq,
       1,
       1,
       {
           {FuncId::dw3, n3, k4},
           {FuncId::dw4, -1, -k4},
           {FuncId::dv3, n3, c3 * c3},
           {FuncId::dv4, -1, c4 * c4},
       }},
      {FuncId::dv2,
       n2 + 1,
       n2 + 1,
       {{FuncId::dv3, -(n2 + 1), (c3 * c3) / (c2 * c2)}}},
      {FuncId::dw2, n2 + 1, n2 + 1, {{FuncId::dw3, -(n2 + 1), 1.0}}},

      // Body conditions
      {FuncId::body1eq,
       0,
       n1,
       {{FuncId::dw1, 1, 0.5}, {FuncId::dw1, 0, 0.5}}},
      {FuncId::body2eq,
       0,
       n2,
       {{FuncId::dw2, 1, 0.5}, {FuncId::dw2, 0, 0.5}}},
      {FuncId::body3eq,
       0,
       n3,
       {{FuncId::dw3, 0, 0.5}, {FuncId::dw3, 1, 0.5}}},
      {FuncId::body4eq,
       0,
       n4,
       {{FuncId::dw4, 0, 0.5}, {FuncId::dw4, 1, 0.5}}},
      {FuncId::reduction1eq,
       0,
       n1,
       {{FuncId::dv1, 0, 0.5}, {FuncId::dv1, 1, 0.5}}},
      {FuncId::reduction2eq,
       0,
       n2,
       {{FuncId::dv2, 0, 0.5}, {FuncId::dv2, 1, 0.5}}},
      {FuncId::reduction3eq,
       0,
       n3,
       {{FuncId::dv3, 0, 0.5}, {FuncId::dv3, 1, 0.5}}},
      {FuncId::reduction4eq,
       0,
       n4,
       {{FuncId::dv4, 0, 0.5}, {FuncId::dv4, 1, 0.5}}},
  };

  return createSystemMatrix(solutionVector, equations, stateVector);
}

inline Eigen::MatrixXd createM1(int n1, int n2, int n3, int n4, double h1,
                                 double h2, double h3, double h4, double c1,
                                 double c2, double c3, double c4) {
  std::vector<Range> stateVector = {
      {FuncId::w1, 1, n1 + 1}, {FuncId::w2, 0, n2},
      {FuncId::w3, 0, n3 + 1}, {FuncId::w4, 0, n4},
      {FuncId::v1, 0, n1 + 1}, {FuncId::v2, 0, n2},
      {FuncId::v3, 0, n3 + 1}, {FuncId::v4, 0, n4 + 1},
  };
  std::vector<Range> solutionVector = pieceSolutionVector(n1, n2, n3, n4);

  std::vector<Eq> equations = {
      {FuncId::w1, 0, 0, {{FuncId::w1, 1, 0.0}}},
      {FuncId::w4, n4 + 1, n4 + 1, {{FuncId::w4, -1, 0.0}}},
      {FuncId::controller1eq, 1, 1, {{FuncId::w1, n1, 0.0}}},
      {FuncId::controller2eq, 1, 1, {{FuncId::w1, n1, 0.0}}},
      {FuncId::controller3eq, 1, 1, {{FuncId::w3, n3, 0.0}}},
      {FuncId::controller4eq, 1, 1, {{FuncId::w3, n3, 0.0}}},
      {FuncId::straineq, 1, 1, {{FuncId::v2, n2, 0.0}}},
      {FuncId::v2,
       n2 + 1,
       n2 + 1,
       {{FuncId::v3, -(n2 + 1), (c3 * c3) / (c2 * c2)}}},
      {FuncId::w2, n2 + 1, n2 + 1, {{FuncId::w3, -(n2 + 1), 1.0}}},

      // Body conditions
      {FuncId::body1eq,
       0,
       n1,
       {{FuncId::v1, 1, -(c1 * c1) / h1}, {FuncId::v1, 0, (c1 * c1) / h1}}},
      {FuncId::body2eq,
       0,
       n2,
       {{FuncId::v2, 1, -(c2 * c2) / h2}, {FuncId::v2, 0, (c2 * c2) / h2}}},
      {FuncId::body3eq,
       0,
       n3,
       {{FuncId::v3, 0, (c3 * c3) / h3}, {FuncId::v3, 1, -(c3 * c3) / h3}}},
      {FuncId::body4eq,
       0,
       n4,
       {{FuncId::v4, 0, (c4 * c4) / h4}, {FuncId::v4, 1, -(c4 * c4) / h4}}},
      {FuncId::reduction1eq,
       0,
       n1,
       {{FuncId::w1, 0, 1.0 / h1}, {FuncId::w1, 1, -1.0 / h1}}},
      {FuncId::reduction2eq,
       0,
       n2,
       {{FuncId::w2, 0, 1.0 / h2}, {FuncId::w2, 1, -1.0 / h2}}},
      {FuncId::reduction3eq,
       0,
       n3,
       {{FuncId::w3, 0, 1.0 / h3}, {FuncId::w3, 1, -1.0 / h3}}},
      {FuncId::reduction4eq,
       0,
       n4,
       {{FuncId::w4, 0, 1.0 / h4}, {FuncId::w4, 1, -1.0 / h4}}},
  };

  return createSystemMatrix(solutionVector, equations, stateVector);
}
