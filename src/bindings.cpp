#include "solver.hpp"
#include "transmission_line_solver.hpp"
#include <nanobind/eigen/dense.h> // enables Eigen <-> numpy conversion
#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_transmission_line_waves, m) {
  m.def("evaluate_solution", &evaluate_solution, nb::arg("n"), nb::arg("t"),
        "Evaluate the solution vector at time t");
  m.def("build_transmission_line_A", &build_transmission_line_A,
        nb::arg("n1"), nb::arg("n2"), nb::arg("n3"), nb::arg("n4"),
        nb::arg("h1"), nb::arg("h2"), nb::arg("h3"), nb::arg("h4"),
        nb::arg("c1"), nb::arg("c2"), nb::arg("c3"), nb::arg("c4"),
        nb::arg("k1"), nb::arg("k2"), nb::arg("k3"), nb::arg("k4"),
        "Build A = M0^-1 M1 (debug helper)");
  m.def("evaluate_transmission_line_wave", &evaluate_transmission_line_wave,
        nb::arg("n1"), nb::arg("n2"), nb::arg("n3"), nb::arg("n4"),
        nb::arg("h1"), nb::arg("h2"), nb::arg("h3"), nb::arg("h4"),
        nb::arg("c1"), nb::arg("c2"), nb::arg("c3"), nb::arg("c4"),
        nb::arg("k1"), nb::arg("k2"), nb::arg("k3"), nb::arg("k4"),
        nb::arg("initial_conditions"), nb::arg("t"), nb::arg("drop") = 4,
        "Evaluate the [w1..w4, v1..v4] state vector at time t");
}
