#pragma once

#include <Eigen/Dense>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Templated on FuncId (an enum identifying named quantities) so this stays
// independent of any particular PDE's variable set.

template <typename FuncId> struct TermTemplate {
  FuncId func;
  int offset; // referenced index is (equation index i) + offset
  double coefficient;
};

template <typename FuncId> struct EquationSpec {
  FuncId lhs;
  int start, end; // inclusive; equation is instantiated for i in [start, end]
  std::vector<TermTemplate<FuncId>> terms;
};

// Declares a contiguous block of named, indexed quantities — used both for
// the raw unknowns (state vector) and for the desired output equations
// (solution vector).
template <typename FuncId> struct RangeSpec {
  FuncId id;
  int start, end; // inclusive
};

template <typename FuncId> struct FuncRef {
  FuncId id;
  int index;
  bool operator==(const FuncRef &other) const {
    return id == other.id && index == other.index;
  }
};

template <typename FuncId> struct FuncRefHash {
  size_t operator()(const FuncRef<FuncId> &ref) const {
    size_t h1 = std::hash<int>()(static_cast<int>(ref.id));
    size_t h2 = std::hash<int>()(ref.index);
    return h1 * 1000003u ^ h2;
  }
};

template <typename FuncId>
using IndexMap = std::unordered_map<FuncRef<FuncId>, int, FuncRefHash<FuncId>>;

// Assigns each (id, i) pair a sequential column/row position, in the order
// the specs and their index ranges are given.
template <typename FuncId>
IndexMap<FuncId> buildIndexMap(const std::vector<RangeSpec<FuncId>> &specs) {
  IndexMap<FuncId> indices;
  int position = 0;
  for (const auto &spec : specs) {
    for (int i = spec.start; i <= spec.end; ++i) {
      indices[FuncRef<FuncId>{spec.id, i}] = position++;
    }
  }
  return indices;
}

template <typename FuncId>
using TermMap =
    std::unordered_map<FuncRef<FuncId>,
                       std::vector<std::pair<FuncRef<FuncId>, double>>,
                       FuncRefHash<FuncId>>;

// Instantiates every EquationSpec over its index range, resolving each
// term's offset against the concrete equation index.
template <typename FuncId>
TermMap<FuncId>
buildTermMap(const std::vector<EquationSpec<FuncId>> &equations) {
  TermMap<FuncId> functionToTerms;
  for (const auto &equation : equations) {
    for (int i = equation.start; i <= equation.end; ++i) {
      std::vector<std::pair<FuncRef<FuncId>, double>> resolvedTerms;
      resolvedTerms.reserve(equation.terms.size());
      for (const auto &term : equation.terms) {
        resolvedTerms.push_back(
            {FuncRef<FuncId>{term.func, i + term.offset}, term.coefficient});
      }
      functionToTerms[FuncRef<FuncId>{equation.lhs, i}] =
          std::move(resolvedTerms);
    }
  }
  return functionToTerms;
}

// Flattens `equations` (which may reference each other) down to a pure
// linear map from raw unknowns (`stateSpecs`) to each output row
// (`solutionSpecs`), inlining any intermediate/derived quantities along
// the way.
template <typename FuncId>
Eigen::MatrixXd
createSystemMatrix(const std::vector<RangeSpec<FuncId>> &solutionSpecs,
                   const std::vector<EquationSpec<FuncId>> &equations,
                   const std::vector<RangeSpec<FuncId>> &stateSpecs) {
  auto stateToColumn = buildIndexMap(stateSpecs);
  auto solutionToRow = buildIndexMap(solutionSpecs);
  auto functionToTerms = buildTermMap(equations);

  Eigen::MatrixXd systemMatrix =
      Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(solutionToRow.size()),
                            static_cast<Eigen::Index>(stateToColumn.size()));

  std::function<void(int, FuncRef<FuncId>, double)> populateSystemMatrix =
      [&](int row, FuncRef<FuncId> f, double coefficient) {
        auto columnIt = stateToColumn.find(f);
        if (columnIt != stateToColumn.end()) {
          systemMatrix(row, columnIt->second) += coefficient;
          return;
        }
        auto termsIt = functionToTerms.find(f);
        if (termsIt != functionToTerms.end()) {
          for (const auto &[subFunc, subCoefficient] : termsIt->second) {
            populateSystemMatrix(row, subFunc, coefficient * subCoefficient);
          }
          return;
        }
        throw std::runtime_error("Definition for function id=" +
                                 std::to_string(static_cast<int>(f.id)) +
                                 " index=" + std::to_string(f.index) +
                                 " was not found.");
      };

  for (const auto &[key, row] : solutionToRow) {
    populateSystemMatrix(row, key, 1.0);
  }

  return systemMatrix;
}
