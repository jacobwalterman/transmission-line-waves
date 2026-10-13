#pragma
#include "solver.hpp"
#include "transmission_line_solver.hpp"
#include <nanobind/eigen/dense.h> // enables Eigen <-> numpy conversion
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

NB_MODULE(_transmission_line_waves, m) {
  nb::class_<Params>(m, "Params")
      .def(nb::init<>())
      .def(nb::init<int, int, int, int, double, double, double, double, double,
                    double, double, double, double, double, double, double>(),
           nb::arg("n1"), nb::arg("n2"), nb::arg("n3"), nb::arg("n4"),
           nb::arg("h1"), nb::arg("h2"), nb::arg("h3"), nb::arg("h4"),
           nb::arg("c1"), nb::arg("c2"), nb::arg("c3"), nb::arg("c4"),
           nb::arg("k1"), nb::arg("k2"), nb::arg("k3"), nb::arg("k4"))
      .def_rw("n1", &Params::n1)
      .def_rw("n2", &Params::n2)
      .def_rw("n3", &Params::n3)
      .def_rw("n4", &Params::n4)
      .def_rw("h1", &Params::h1)
      .def_rw("h2", &Params::h2)
      .def_rw("h3", &Params::h3)
      .def_rw("h4", &Params::h4)
      .def_rw("c1", &Params::c1)
      .def_rw("c2", &Params::c2)
      .def_rw("c3", &Params::c3)
      .def_rw("c4", &Params::c4)
      .def_rw("k1", &Params::k1)
      .def_rw("k2", &Params::k2)
      .def_rw("k3", &Params::k3)
      .def_rw("k4", &Params::k4);
  nb::class_<EigenSystem>(m, "EigenSystem")
      .def("sortedEigenvalues",
           [](const EigenSystem &es) { return es.sortedEigenvalues; })
      .def("sortedEigenvectors",
           [](const EigenSystem &es) { return es.sortedEigenvectors; });
  m.def("decompose", &decompose);
  m.def("coefficients", &coefficients);
  m.def("evolve", &evolve);
  m.def("evaluate_solution", &evaluate_solution, nb::arg("n"), nb::arg("t"),
        "Evaluate the solution vector at time t");
  m.def("build_transmission_line_A", &build_transmission_line_A, nb::arg("n1"),
        nb::arg("n2"), nb::arg("n3"), nb::arg("n4"), nb::arg("h1"),
        nb::arg("h2"), nb::arg("h3"), nb::arg("h4"), nb::arg("c1"),
        nb::arg("c2"), nb::arg("c3"), nb::arg("c4"), nb::arg("k1"),
        nb::arg("k2"), nb::arg("k3"), nb::arg("k4"),
        "Build A = M0^-1 M1 (debug helper)");
  m.def("evaluate_transmission_line_wave", &evaluate_transmission_line_wave,
        nb::arg("n1"), nb::arg("n2"), nb::arg("n3"), nb::arg("n4"),
        nb::arg("h1"), nb::arg("h2"), nb::arg("h3"), nb::arg("h4"),
        nb::arg("c1"), nb::arg("c2"), nb::arg("c3"), nb::arg("c4"),
        nb::arg("k1"), nb::arg("k2"), nb::arg("k3"), nb::arg("k4"),
        nb::arg("initial_conditions"), nb::arg("t"), nb::arg("drop") = 4,
        "Evaluate the [w1..w4, v1..v4] state vector at time t");
}
