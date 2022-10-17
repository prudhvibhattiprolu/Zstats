#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "FCTools.h"

namespace py = pybind11;

PYBIND11_MODULE(FC, mod)
{
    
    mod.def("UpperLimit", &UpperLimit, R"pbdoc(
    Inputs
    ------

    n: Observed number of Poisson events
    b: Poisson mean of the background
    CL: Confidence level that is set to 0.90 (90% CL) by default
    more_info: Set to *True* to print out more information about the
               computation. Set to *False* by default.
               (Also see the C++ code for more details.)

    Returns
    -------

    Computes FC upper limit given number of observed events (n),
    background mean (b) at a chosen confidence level (CL)

    References
    ----------

    G. J. Feldman, R. D. Cousins [arXiv:physics/9711021]
    )pbdoc", py::arg("n"), py::arg("b"), py::arg("CL") = 0.9,
    py::arg("more_info") = 0);
    
    mod.def("Sensitivity", &Sensitivity, R"pbdoc(
    Inputs
    ------

    b: Poisson mean of the background
    CL: Confidence level that is set to 0.90 (90% CL) by default
    more_info: Set to *True* to print out more information about the
               computation. Set to *False* by default.
               (Also see the C++ code for more details.)

    Returns
    -------

    Computes FC experimental sensitivity given a background mean (b) at a
    chosen confidence level (CL)

    References
    ----------

    G. J. Feldman, R. D. Cousins [arXiv:physics/9711021]
    )pbdoc", py::arg("b"), py::arg("CL") = 0.9, py::arg("more_info") = 0);
};
