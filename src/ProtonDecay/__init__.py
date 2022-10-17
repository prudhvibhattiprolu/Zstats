import sys
import itertools
import numpy as np

#Import local packages

sys.path.append('../')
import ProtonDecayTools as CppTools
import FC #C++ tools (with python bindings) to compute Feldman-Cousins upper limits on signal and experimental sensitivity as detailed in arXiv:physics/9711021
import Zstats
del sys.path[-1]

############### Single-channel proton decay experiments ###############

#Computes proton partial lifetime given the signal, exposure in kton.years for various detector materials
def ProtonLifetime(exposure, s, eps, Np):
    '''
    Inputs
    ------

    exposure: exposure in kiloton-years

    s: signal needed for discovery/exclusion

    eps: signal selection efficiency (0 <= eps <= 1)

    Np: number of protons per kiloton of detector material

    Returns
    -------

    Proton partial lifetime reach/limit in years at a single-channel proton
    decay experiment for a given exposure of the experiment in kiloton-years,
    signal needed for (expected) discovery/exclusion, and the
    signal selection efficiency that is taken to be perfectly known,
    and number of protons/kiloton that can be estimated based on the detector
    material used

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Proton partial lifetime in years
    proton_lifetime = Np*eps*exposure/s
    return proton_lifetime

#Computes the expected exclusion reach for proton lifetime in years
def LifetimeExclExp_Unc0(b, eps, exposure, Np, CL=0.9, CLExclbool=True):
    '''
    Inputs
    ------

    b: estimated background at the experiment

    eps: signal selection efficiency (0 <= eps <= 1)

    exposure: exposure in kiloton-years

    Np: number of protons per kiloton of detector material

    CL: confidence level for an expected exclusion

    CLExclbool: Set to *False* by default. Set this to *True* to use
                a more conservative Bayesian statistic CLExcl instead of
                the standard frequentist pExcl

    Returns
    -------

    Exclusion reach for proton partial lifetime in years at (planned) future
    proton decay experiments (such as DUNE, JUNO, THEIA) that treat the
    estimated background and signal efficiency to be perfectly known while
    modeling the experiment as a single-channel counting experiment

    Examples
    --------

    Details of a few planned experiments are given below:

        DUNE:

            detector material: liquid argon
            Np = 2.71 10^32 protons/kiloton
            Nkton = 40 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at DUNE are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                b/Mton-year = 0.25, 1, 2.5
                eps = 0.4 +/- 0.1

                Zero background limit: (b, eps) = (0, 0.46)

                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (15, 0.47)

        JUNO:

            detector material: liquid scintillator
            Np = 6.75 10^33/20 protons/kiloton
            Nkton = 20 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at JUNO are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = {(2.5, 0.55), (1.5, 0.26)}

        THEIA-25/100:

            detector material: water based liquid scintillator
            Np = 3.35 10^32 protons/kiloton
            Nkton = 17 (80) kiloton fiducial mass of detector material at
                    THEIA-25 (THEIA-100)

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at THEIA-25/THEIA-100 are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = (2.5, 0.55)


                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (0.3, 0.40)

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Expected upper limit on signal at chosen CL
    s = Zstats.sExclExp(bhat=b, CL=CL, CLExclbool=CLExclbool)

    return ProtonLifetime(exposure, s, eps, Np)

#Computes the expected discovery reach for proton lifetime in years for a chosen runtime in years
def LifetimeDiscExp_Unc0(b, eps, exposure, Np, Z=3, CLDiscbool=True, s_min=1):
    '''
    Inputs
    ------

    b: estimated background at the experiment

    eps: signal selection efficiency (0 <= eps <= 1)

    exposure: exposure in kiloton-years

    Np: number of protons per kiloton of detector material

    Z: significance for an expected discovery

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLDisc

    s_min: minimum signal needed to claim a discovery. Set 1 by default.

    Returns
    -------

    Discovery reach for proton partial lifetime in years at (planned) future
    proton decay experiments (such as DUNE, JUNO, THEIA) that treat the
    estimated background and signal efficiency to be perfectly known while
    modeling the experiment as a single-channel counting experiment

    Examples
    --------

    Details of a few planned experiments are given below:

        DUNE:

            detector material: liquid argon
            Np = 2.71 10^32 protons/kiloton
            Nkton = 40 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at DUNE are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                b/Mton-year = 0.25, 1, 2.5
                eps = 0.4 +/- 0.1

                Zero background limit: (b, eps) = (0, 0.46)

                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (15, 0.47)

        JUNO:

            detector material: liquid scintillator
            Np = 6.75 10^33/20 protons/kiloton
            Nkton = 20 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at JUNO are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = {(2.5, 0.55), (1.5, 0.26)}

        THEIA-25/100:

            detector material: water based liquid scintillator
            Np = 3.35 10^32 protons/kiloton
            Nkton = 17 (80) kiloton fiducial mass of detector material at
                    THEIA-25 (THEIA-100)

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at THEIA-25/THEIA-100 are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = (2.5, 0.55)


                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (0.3, 0.40)

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #signal needed for an expected discovery at chosen Z
    s = Zstats.sDiscExp(bhat=b, Z=Z, CLDiscbool=CLDiscbool) if b != 0 else s_min

    if s < s_min:
        s = s_min

    return ProtonLifetime(exposure, s, eps, Np)

#Computes experiment runtime given upper limit on signal (for exclusion) or expected signal (for discovery), proton partial lifetime in years, etc.
def Runtime(proton_lifetime, s, eps, Nkton, Np):
    '''
    Inputs
    ------

    proton_lifetime: proton partial lifetime in years

    s: signal needed for expected discovery/exclusion

    eps: signal selection efficiency (0 <= eps <= 1)

    Nkton: number of kiltons of detector material

    Np: number of protons per kiloton of detector material

    Returns
    -------

    Computes the runtime of an experiment in years needed for an expected
    discovery/exclusion at a single-channel proton decay experiment for a given
    proton partial lifetime in years, signal needed for expected
    discovery/exclusion, signal selection efficiency that is taken to be
    perfectly known, number of kiltons of detector material, and number of
    protons/kiloton that can be estimated based on the detector material used

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Expt run time in years
    time_excl = s*proton_lifetime/(Np*Nkton*eps)

    return time_excl

#Computes the expected runtime in years required to exclude proton decay for a chosen partial lifetime in years
def RuntimeExclExp_Unc0(proton_lifetime, b, eps, Nkton, Np, CL=0.9, CLExclbool=True):
    '''
    Inputs
    ------

    proton_lifetime: proton partial lifetime in years

    b: estimated background at the experiment

    eps: signal selection efficiency (0 <= eps <= 1)

    Nkton: number of kiltons of detector material

    Np: number of protons per kiloton of detector material

    CL: confidence level for an expected exclusion

    CLExclbool: Set to *False* by default. Set this to *True* to use
                a more conservative Bayesian statistic CLExcl instead of
                the standard frequentist pExcl

    Returns
    -------

    Runtime in years needed for an expected exclusion of proton partial
    lifetime in years at a chosen confidence level at (planned) future
    proton decay experiments (such as DUNE, JUNO, THEIA) that treat the
    estimated background and signal efficiency to be perfectly known while
    modeling the experiment as a single-channel counting experiment

    Examples
    --------

    Details of a few planned experiments are given below:

        DUNE:

            detector material: liquid argon
            Np = 2.71 10^32 protons/kiloton
            Nkton = 40 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at DUNE are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                b/Mton-year = 0.25, 1, 2.5
                eps = 0.4 +/- 0.1

                Zero background limit: (b, eps) = (0, 0.46)

                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (15, 0.47)

        JUNO:

            detector material: liquid scintillator
            Np = 6.75 10^33/20 protons/kiloton
            Nkton = 20 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at JUNO are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = {(2.5, 0.55), (1.5, 0.26)}

        THEIA-25/100:

            detector material: water based liquid scintillator
            Np = 3.35 10^32 protons/kiloton
            Nkton = 17 (80) kiloton fiducial mass of detector material at
                    THEIA-25 (THEIA-100)

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at THEIA-25/THEIA-100 are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = (2.5, 0.55)


                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (0.3, 0.40)

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Expected upper limit on signal at chosen CL
    s = Zstats.sExclExp(bhat=b, CL=CL, CLExclbool=CLExclbool)

    return Runtime(proton_lifetime, s, eps, Nkton, Np)

#Computes the expected runtime in years required to discover proton decay for a chosen partial lifetime in years
def RuntimeDiscExp_Unc0(proton_lifetime, b, eps, Nkton, Np, Z=3, CLDiscbool=True, s_min=1):
    '''
    Inputs
    ------

    proton_lifetime: proton partial lifetime in years

    b: estimated background at the experiment

    eps: signal selection efficiency (0 <= eps <= 1)

    Nkton: number of kiltons of detector material

    Np: number of protons per kiloton of detector material

    Z: significance for an expected discovery

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLDisc

    s_min: minimum signal needed to claim a discovery. Set 1 by default.

    Returns
    -------

    Runtime in years needed for an expected discovery of proton partial
    lifetime in years at a chosen significance at (planned) future
    proton decay experiments (such as DUNE, JUNO, THEIA) that treat the
    estimated background and signal efficiency to be perfectly known while
    modeling the experiment as a single-channel counting experiment

    Examples
    --------

    Details of a few planned experiments are given below:

        DUNE:

            detector material: liquid argon
            Np = 2.71 10^32 protons/kiloton
            Nkton = 40 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at DUNE are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                b/Mton-year = 0.25, 1, 2.5
                eps = 0.4 +/- 0.1

                Zero background limit: (b, eps) = (0, 0.46)

                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (15, 0.47)

        JUNO:

            detector material: liquid scintillator
            Np = 6.75 10^33/20 protons/kiloton
            Nkton = 20 kiloton fiducial mass of detector material

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at JUNO are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = {(2.5, 0.55), (1.5, 0.26)}

        THEIA-25/100:

            detector material: water based liquid scintillator
            Np = 3.35 10^32 protons/kiloton
            Nkton = 17 (80) kiloton fiducial mass of detector material at
                    THEIA-25 (THEIA-100)

            Plausible background rates and signal efficiencies for various
            proton decay modes that can probed at THEIA-25/THEIA-100 are

                p -> \overline nu +  K^+ decay mode:
                ------------------------------------

                (b/Mton-year, eps) = (2.5, 0.55)


                p -> e^+ + \pi^0 decay mode:
                ----------------------------

                (b/Mton-year, eps) = (0.3, 0.40)

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #signal needed for an expected discovery at chosen Z
    s = Zstats.sDiscExp(bhat=b, Z=Z, CLDiscbool=CLDiscbool) if b != 0 else s_min

    if s < s_min:
        s = s_min

    return Runtime(proton_lifetime, s, eps, Nkton, Np)

def _DiscSignal(b, Z, CLDiscbool=True):
    '''
    Inputs
    ------

    b: estimated background at the experiment

    Z: significance for an expected discovery

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLDisc

    Returns
    -------

    Signal needed for discovery before imposing a minimum signal needed for
    discovery

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #signal needed for an expected discovery at chosen Z
    s = Zstats.sDiscExp(bhat=b, Z=Z, CLDiscbool=CLDiscbool)

    return s

############### Multi-channel proton decay experiments ###############

def LifetimeExclExpN(bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, CL=0.9):
    '''
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

    Returns
    -------

    Exclusion reach for proton partial lifetime in years at a chosen confidence
    level CL at (planned) future proton decay experiments (such as
    Hyper-Kamiokande) with backgrounds and signal efficiencies along with their
    corresponding uncertainties estimated while modeling the experiment as a
    multi-channel counting experiment

    Examples
    --------

    Details of Hyper-Kamiokande are given below:

        detector material: water
        Np = 3.34 10^32 protons/kiloton
        Nkton = 186 kiloton fiducial mass of detector material

        Estimated background rates and signal efficiencies in various
        search channels along with their corresponding uncertainties for
        various proton decay modes that can probed at Hyper-Kamiokande are

            p -> \overline nu +  K^+ decay mode:
            ------------------------------------

            b/Mton-year = [1916, 0.9, 0.7]

            db/Mton-year = [0.0, 0.2, 0.2]

            eps = [0.31, 0.127, 0.108]

            deps = [0, 0.024, 0.011]

            p -> e^+ + \pi^0 decay mode:
            ----------------------------

            b/Mton-year = [0.06, 0.62]

            db/Mton-year = [0.02, 0.2]

            eps = [0.187, 0.194]

            deps = [0.012, 0.029]

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    Width_Unc0 = CppTools.WidthExclExp(bhatarray, [0.0]*len(bhatarray), epsarray, [0.0]*len(bhatarray), exposurearray, Np, CL, Np/1e34, 15, 50)

    Width = CppTools.WidthExclExp(bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, CL, Np*Width_Unc0, 2, 25)

    return 1/Width

def LifetimeDiscExpN(bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, Z=3):
    '''
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

    Returns
    -------

    Discovery reach for proton partial lifetime in years at a chosen
    significance at (planned) future proton decay experiments (such as
    Hyper-Kamiokande) with backgrounds and signal efficiencies along with their
    corresponding uncertainties estimated while modeling the experiment as a
    multi-channel counting experiment

    Examples
    --------

    Details of Hyper-Kamiokande are given below:

        detector material: water
        Np = 3.34 10^32 protons/kiloton
        Nkton = 186 kiloton fiducial mass of detector material

        Estimated background rates and signal efficiencies in various
        search channels along with their corresponding uncertainties for
        various proton decay modes that can probed at Hyper-Kamiokande are

            p -> \overline nu +  K^+ decay mode:
            ------------------------------------

            b/Mton-year = [1916, 0.9, 0.7]

            db/Mton-year = [0.0, 0.2, 0.2]

            eps = [0.31, 0.127, 0.108]

            deps = [0, 0.024, 0.011]

            p -> e^+ + \pi^0 decay mode:
            ----------------------------

            b/Mton-year = [0.06, 0.62]

            db/Mton-year = [0.02, 0.2]

            eps = [0.187, 0.194]

            deps = [0.012, 0.029]

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    Width_Unc0 = CppTools.WidthDiscExp(bhatarray, [0.0]*len(bhatarray), epsarray, [0.0]*len(bhatarray), exposurearray, Np, Z, Np/1e34, 15, 50)

    Width = CppTools.WidthDiscExp(bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, Z, Np*Width_Unc0, 2, 25)

    return 1/Width

def ConfidenceLevelExclObsN(proton_lifetime, narray, bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np):
    '''
    Inputs
    ------

    proton_lifetime: proton partial lifetime in years

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

    Returns
    -------

    The confidence level for (observed) exclusion of proton partial
    lifetime in years based on the data collected at a multi-channel proton
    decay experiment (such as Super-Kamiokande)

    Examples
    --------

    (Published) data collected so far at Super-Kamiokande is given below:

        detector material: water
        Np = 3.34 10^32 protons/kiloton

        Data along with the estimated background rates and signal
        efficiencies in various search channels along with their
        corresponding uncertainties for various proton decay modes that
        can probed at Super-Kamiokande are

            p -> \overline nu +  K^+ decay mode:
            ------------------------------------
            //Observed number of events
            n = [0, 0, 0, 0,
                 0, 0, 0, 0,
                 177, 78, 85, 226]

            //Expected backgrounds
            b = [0.07, 0.14, 0.03, 0.13,
                 0.18, 0.17, 0.07, 0.17,
                 193.21, 94.27, 69.0, 223.14]

            //Expected background uncertainties
            db = [0.02, 0.03, 0.01, 0.03,
                  0.04, 0.03, 0.01, 0.03,
                  3.58, 1.72, 1.28, 4.1]

            //Signal selection efficiencies
            eps = [0.079, 0.063, 0.077, 0.091,
                   0.078, 0.067, 0.079, 0.1,
                   0.339, 0.306, 0.326, 0.376]

            //Uncertainties on the signal selection efficiencies
            deps = [0.001, 0.001, 0.001, 0.001,
                    0.001, 0.001, 0.001, 0.001,
                    0.003, 0.003, 0.003, 0.003]

            //Corresponding exposures in kiloton-years
            exposure = [91.7, 49.2, 31.9, 87.3,
                        91.7, 49.2, 31.9, 87.3,
                        91.7, 49.2, 31.9, 87.3]

            p -> e^+ + \pi^0 decay mode:
            ----------------------------
            //Observed number of events
            n = [0, 0, 0, 0,
                 0, 0, 0, 0]

            //Expected backgrounds
            b = [0.01, 0.01, 0.01, 0.01,
                 0.15, 0.11, 0.07, 0.25]

            //Expected background uncertainties
            db = [0.01, 0.01, 0, 0,
                  0.06, 0.04, 0.03, 0.11]

            //Signal selection efficiencies
            eps = [0.183, 0.166, 0.187, 0.182,
                   0.2, 0.194, 0.203, 0.192]

            //Uncertainties on the signal selection efficiencies
            deps = [0.017, 0.017, 0.017, 0.015,
                    0.033, 0.03, 0.033, 0.031]

            //Corresponding exposures in kiloton-years
            exposure = [111.4, 59.4, 38.6, 241.3,
                        111.4, 59.4, 38.6, 241.3]

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    return 1 - CppTools.CLExcl(narray, bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np/proton_lifetime)

def LifetimeExclObsN(narray, bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, CL=0.9):
    '''
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

    Returns
    -------

    Exclusion reach for proton partial lifetime in years at a chosen confidence
    level CL at (planned) future proton decay experiments (such as
    Hyper-Kamiokande) with backgrounds and signal efficiencies along with their
    corresponding uncertainties estimated while modeling the experiment as a
    multi-channel counting experiment

    Examples
    --------

    Details of Hyper-Kamiokande are given below:

        detector material: water
        Np = 3.34 10^32 protons/kiloton
        Nkton = 186 kiloton fiducial mass of detector material

        Estimated background rates and signal efficiencies in various
        search channels along with their corresponding uncertainties for
        various proton decay modes that can probed at Hyper-Kamiokande are

            p -> \overline nu +  K^+ decay mode:
            ------------------------------------

            b/Mton-year = [1916, 0.9, 0.7]

            db/Mton-year = [0.0, 0.2, 0.2]

            eps = [0.31, 0.127, 0.108]

            deps = [0, 0.024, 0.011]

            p -> e^+ + \pi^0 decay mode:
            ----------------------------

            b/Mton-year = [0.06, 0.62]

            db/Mton-year = [0.02, 0.2]

            eps = [0.187, 0.194]

            deps = [0.012, 0.029]

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    Width_Unc0 = CppTools.WidthExclObs(narray, bhatarray, [0.0]*len(bhatarray), epsarray, [0.0]*len(bhatarray), exposurearray, Np, CL, Np/1e34, 10, 60)

    Width = CppTools.WidthExclObs(narray, bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, CL, Np*Width_Unc0, 2, 25)

    return 1/Width

## Expected signal based on Hyper-K's methods
#def sExclExpHyperK(b, CL=0.9):
#    maxcount = Zstats._MaxPoissonCount(0, b)
#    result = 0.0
#    for n in range(maxcount + 1):
#        stemp = Zstats.sExclObs(n=n, bhat=b, CL=CL, CLExclbool=True)
#        result = result + Zstats.DeltaP(n, np.inf, np.inf, 0, b)*(1/stemp)
#    return 1/result
#
#def LifetimeExclExpNHyperK(bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, CL=0.9):
#
#    probn = lambda narraytemp: CppTools.DeltaPN(narraytemp, bhatarray, dbhatarray, [0.0]*len(bhatarray), [0.0]*len(bhatarray), [0.0]*len(bhatarray), 0.0)
#    lifetimen = lambda narraytemp: LifetimeExclObsN(narraytemp, bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, CL)
#
#    probnavg = probn(bhatarray)
#    probmin = probnavg/(1e5*len(bhatarray))
#
#    result = 0.0
#
#    if len(bhatarray) == 1:
#
#        n1max = Zstats._MaxPoissonCount(0.0, bhatarray[0])
#
#        for n in range(n1max + 1):
#            probstep = probn([n])
#            if probstep >= probmin:
#                lifetimestep = lifetimen([n])
#                result = result + probstep*lifetimestep
#
#    if len(bhatarray) == 2:
#
#        n1max = Zstats._MaxPoissonCount(0.0, bhatarray[0])
#        n2max = Zstats._MaxPoissonCount(0.0, bhatarray[1])
#
#        narraylist = itertools.product(range(n1max + 1), range(n2max + 1))
#
#        for narray in narraylist:
#            probstep = probn(list(narray))
#            if probstep >= probmin:
#                lifetimestep = lifetimen(list(narray))
#                result = result + probstep*lifetimestep
#
#    if len(bhatarray) == 3:
#
#        n1max = Zstats._MaxPoissonCount(0.0, bhatarray[0])
#        n2max = Zstats._MaxPoissonCount(0.0, bhatarray[1])
#        n3max = Zstats._MaxPoissonCount(0.0, bhatarray[2])
#
#        narraylist = itertools.product(range(n1max + 1), range(n2max + 1), range(n3max + 1))
#
#        for narray in narraylist:
#            probstep = probn(list(narray))
#            if probstep >= probmin:
#                lifetimestep = lifetimen(list(narray))
#                result = result + probstep*lifetimestep
#
#    return result
