#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "BayesianTools.h"

namespace py = pybind11;

PYBIND11_MODULE(ProtonDecayTools, mod)
{
    mod.def("DeltaPN", &DeltaPN, R"pbdoc(
    Inputs
    ------

    narray: array of number of observed events at a multi-channel proton
            decay experiment

    bhatarray: array of background means at a multi-channel
               proton decay experiment

    dbhatarray: array of background uncertainties at a multi-channel
                proton decay experiment

    epsarray: array of signal efficiencies at a multi-channel proton
              decay experiment

    depsarray: array of the uncertainties in signal efficiencies at a
               multi-channel proton decay experiment

    exposurearray: array of the exposures in kiloton-years at a
                   multi-channel proton decay experiment

    NpbyLifetime: Number of protons per kilton divided by proton partial
                  lifetime in years

    Returns
    -------

    Computes the product of Poisson probabilities in each channel that is
    marginalized over th background and signal efficiency to take the
    corresponding undertainties into account

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    )pbdoc", py::arg("narray"), py::arg("bhatarray"), py::arg("dbhatarray"),
    py::arg("epsarray"), py::arg("depsarray"), py::arg("exposurearray"),
    py::arg("NpbyLifetime"));
    
    mod.def("CLExcl", &CLExcl, R"pbdoc(
    Inputs
    ------

    narray: array of number of observed events at a multi-channel proton
            decay experiment

    bhatarray: array of background means at a multi-channel
               proton decay experiment

    dbhatarray: array of background uncertainties at a multi-channel
                proton decay experiment

    epsarray: array of signal efficiencies at a multi-channel proton
              decay experiment

    depsarray: array of the uncertainties in signal efficiencies at a
               multi-channel proton decay experiment

    exposurearray: array of the exposures in kiloton-years at a
                   multi-channel proton decay experiment

    NpbyLifetime: Number of protons per kilton divided by proton partial
             lifetime in years

    Returns
    -------

    Computes the Bayesian statistic CLExcl at a multi-channel proton decay
    experiment which can be interpreted as (1 - Confidence Level)

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    )pbdoc", py::arg("narray"), py::arg("bhatarray"), py::arg("dbhatarray"),
    py::arg("epsarray"), py::arg("depsarray"), py::arg("exposurearray"),
    py::arg("NpbyLifetime"), py::arg("Normalized") = 1);
    
    mod.def("CLDisc", &CLDisc, R"pbdoc(
    Inputs
    ------

    narray: array of number of observed events at a multi-channel proton
            decay experiment

    bhatarray: array of background means at a multi-channel
               proton decay experiment

    dbhatarray: array of background uncertainties at a multi-channel
                proton decay experiment

    epsarray: array of signal efficiencies at a multi-channel proton
              decay experiment

    depsarray: array of the uncertainties in signal efficiencies at a
               multi-channel proton decay experiment

    exposurearray: array of the exposures in kiloton-years at a
                   multi-channel proton decay experiment

    NpbyLifetime: number of protons per kilton divided by proton partial
             lifetime in years

    Returns
    -------

    Computes the Bayesian statistic CLDisc at a multi-channel proton decay
    experiment which can be used instead of the p-value for discovery

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    )pbdoc", py::arg("narray"), py::arg("bhatarray"), py::arg("dbhatarray"),
    py::arg("epsarray"), py::arg("depsarray"), py::arg("exposurearray"),
        py::arg("NpbyLifetime"));
    
    mod.def("WidthExclObs", &WidthExclObs, R"pbdoc(
    Inputs
    ------

    narray: array of number of observed events at a multi-channel proton
            decay experiment

    bhatarray: array of background means at a multi-channel
               proton decay experiment

    dbhatarray: array of background uncertainties at a multi-channel
                proton decay experiment

    epsarray: array of signal efficiencies at a multi-channel proton
              decay experiment

    depsarray: array of the uncertainties in signal efficiencies at a
               multi-channel proton decay experiment

    exposurearray: array of the exposures in kiloton-years at a
                   multi-channel proton decay experiment

    Np: number of protons per kiloton of detector material

    CL: confidence level for expected exclusion

    NpbyLifetime_guess: initial guess for number of protons per kiloton divided by
                   proton partial lifetime in years used by the root solving
                   algorithm

    rootfactor: factor used by the root solving algorithm to solve for the root
                that is ideally in between NpbyLifetime_guess/rootfactor and 
                NpbyLifetime_guess*rootfactor

    niter: maximum number of iterations while solving for the root

    Returns
    -------
    
    Computes the observed upper limit on proton partial width in year^-1 at a
    at chosen confidence level CL at a multi-channel proton decay experiment

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    )pbdoc", py::arg("narray"), py::arg("bhatarray"), py::arg("dbhatarray"), py::arg("epsarray"),
    py::arg("depsarray"), py::arg("exposurearray"), py::arg("Np"),
        py::arg("CL"), py::arg("NpbyLifetime_guess"), py::arg("rootfactor"),
        py::arg("niter"));
    
    mod.def("WidthExclExp", &WidthExclExp, R"pbdoc(
    Inputs
    ------

    bhatarray: array of background means at a multi-channel
               proton decay experiment

    dbhatarray: array of background uncertainties at a multi-channel
                proton decay experiment

    epsarray: array of signal efficiencies at a multi-channel proton
              decay experiment

    depsarray: array of the uncertainties in signal efficiencies at a
               multi-channel proton decay experiment

    exposurearray: array of the exposures in kiloton-years at a
                   multi-channel proton decay experiment

    Np: number of protons per kiloton of detector material

    CL: confidence level for expected exclusion

    NpbyLifetime_guess: initial guess for number of protons per kiloton divided by
                   proton partial lifetime in years used by the root solving
                   algorithm

    rootfactor: factor used by the root solving algorithm to solve for the root
                that is ideally in between NpbyLifetime_guess/rootfactor and 
                NpbyLifetime_guess*rootfactor

    niter: maximum number of iterations while solving for the root

    Returns
    -------
    
    Computes the expected upper limit on proton partial width in year^-1 at a
    at chosen confidence level CL at a multi-channel proton decay experiment

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    )pbdoc", py::arg("bhatarray"), py::arg("dbhatarray"), py::arg("epsarray"),
    py::arg("depsarray"), py::arg("exposurearray"), py::arg("Np"),
        py::arg("CL"), py::arg("NpbyLifetime_guess"), py::arg("rootfactor"),
        py::arg("niter"));
    
    mod.def("WidthDiscExp", &WidthDiscExp, R"pbdoc(
    Inputs
    ------

    bhatarray: array of background means at a multi-channel
               proton decay experiment

    dbhatarray: array of background uncertainties at a multi-channel
                proton decay experiment

    epsarray: array of signal efficiencies at a multi-channel proton
              decay experiment

    depsarray: array of the uncertainties in signal efficiencies at a
               multi-channel proton decay experiment

    exposurearray: array of the exposures in kiloton-years at a
                   multi-channel proton decay experiment

    Np: number of protons per kiloton of detector material

    Z: significance for expected discovery

    NpbyLifetime_guess: initial guess for number of protons per kiloton divided by
                   proton partial lifetime in years used by the root solving
                   algorithm

    rootfactor: factor used by the root solving algorithm to solve for the root
                that is ideally in between NpbyLifetime_guess/rootfactor and 
                NpbyLifetime_guess*rootfactor

    niter: maximum number of iterations while solving for the root

    Returns
    -------
    
    Computes the proton partial width in year^-1 needed for an expected
    discovery at a chosen significance Z at a multi-channel proton decay
    experiment

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    )pbdoc", py::arg("bhatarray"), py::arg("dbhatarray"), py::arg("epsarray"),
    py::arg("depsarray"), py::arg("exposurearray"), py::arg("Np"),
        py::arg("Z"), py::arg("NpbyLifetime_guess"), py::arg("rootfactor"),
        py::arg("niter"));
    
};
