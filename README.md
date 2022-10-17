
# Zstats v2.0

### Statistical measures for discovery and exclusion of new physics signals at multi-channel counting experiments

### Prudhvi N. Bhattiprolu, Stephen P. Martin, James D. Wells

This package provides the tools to compute the observed and projected discovery and
exclusion significances at counting experiments, possibly with multiple independent
search channels and uncertain backgrounds, based on References:
[arXiv:2009.07249 [physics.data-an]](https://arxiv.org/abs/2009.07249) (Ref. [1]) and
[arXiv:2210.07735 [hep-ph]](https://arxiv.org/abs/2210.07735) (Ref. [2]).
Also included are tools to compute the (expected) upper limits on the signal
and the signal needed for an expected discovery using various methods.
As an application, these tools are used to study the statistical
significances for proton decay experiments, including the experiments with various
independent search channels with possibly uncertain backgrounds and signal efficiencies.

Ref. [1] advocated for the "exact Asimov significance" as the standard significance
measure for projected exclusions and discovery sensitivities in counting experiments.
Ref. [2] generalized the methods in Ref. [1] to include Bayesian and modified frequentist
statistics, and multi-channel counting experiments, with application to proton decay experiments
at current and future neutrino detectors. Ref. [2] argued in favor of conservative
Bayesian-motivated statistical measures which do not suffer from various counterintuitive flaws
associated with (modified) frequentist statistical measures, especially for multi-channel
counting experiments.


## Installation

The installation requires Python3 with the following packages (that are not built-in):

* numpy
* scipy
* mpmath

and a C++(11) compiler (e.g. gcc). (The installation also requires some C++ headers from [boost](https://www.boost.org/) and [pybind11](https://pybind11.readthedocs.io/en/stable/) which are already built into this package, located at ` src/ProtonDecay/BayesianCpp/boost` and `extern/pybind11`, and need not be installed seperately.)

### Automatic (**recommended**)

To automatically install Zstats as a python package using pip, run the following from command line

```bash
python -m pip install "git+https://github.com/prudhvibhattiprolu/Zstats.git#egg=Zstats"
```

and the package can be readily imported from python.

### Manual

Instead, if you want to edit the source code and/or build it into your own code/package (or if the automatic installation did not successfully go through) the package can also be manually cloned and installed, by doing the following

```bash
# Git clone
git clone https://github.com/prudhvibhattiprolu/Zstats.git

# cd into the package
cd Zstats
```
At this stage, if needed, you can make edits to the Python code at `src/Zstats/__init__.py` and `src/ProtonDecay/__init__.py`, and C++ headers/code at `src/ProtonDecay/BayesianCpp/`
and `src/FeldmanCousins/`, and then build the package:

```bash
# build the package manually by doing
python setup.py --quiet build_ext --inplace clean --all
```

 If the installation is successful, two new files, namely `src/FC.*.so` and `src/ProtonDecayTools.*.so` files should appear. To start using the package in python, `cd src/` and launch python to import Zstats
 as a package.
 
(Note that the C++(11) code is only used in the modules `Zstats.ProtonDecay` and `Zstats.FC` but not
in functions in the main module. So if those modules are not required, one can also edit the setup.py file to remove those modules from the build as done in `src/_setup_nocpp.py`.)

## Usage

To start using the package

```python
import Zstats
```

All the functions (along with the documentation) in the main package can be listed by using the Python help function 

```python
help(Zstats)
```

The documentation for each function can be accessed by doing

```python
# To access the documentation of the function Zstats.ZDiscExp, say, we can either do
help(Zstats.ZDiscExp)
# or
print(Zstats.ZDiscExp.__doc__)
```

The package also contains two additional modules `Zstats.ProtonDecay` and `Zstats.FC` and the documentation of each of the function within these modules can also be accessed in the same way as shown above.

The main functions in the package `Zstats` and its modules: `Zstats.ProtonDecay` and `Zstats.FC` are
listed below with a brief description. (The full functionality can easily be looked up using the
Python help function.)

## Main functionality

### `Zstats` package:

This package primarily computes the observed and expected significances for
discovery and exclusion possibly with multiple independent search channels
and background uncertainties. Additionally, for single-channel experiments,
tools are also provided to compute the observed and expected upper limits on the signal,
along with the number of observed events needed for discovery, and the signal needed
for an expected discovery.

Definitions of the input parameters used in the functions below:

* `n`: number of observed events in the signal region
* `bhat`: Poisson mean of the background in the signal region
* `dbhat`: background uncertainty in the signal region
* `s`: Poisson mean of the signal
* `Z`: significance needed for discovery (e.g. 5 sigma)
* `CL`: confidence level for exclusion (e.g. 95% CL)

For the functions that compute the significances for multi-channel experiments,
the inputs should be *python lists* or *numpy arrays* which will be
referred to as (`narray`, `bhatarray`, `dbhatarray`, `sarray`) here. (Also in the case of multi-channel
counting experiments with uncertain backgrounds, we take the probability distribution
of the true background in the multi-channel experiment to be the product of
probability distributions of the true background in each of the search channel.)

All the functions use conservative Bayesian-motivated statistical measures
`CLExcl` and `CLDisc`, defined and argued for in Ref. [2],
unless explicitly mentioned in the function's docstring.
If required, most of functions in this package, especially the ones for single-channel experiments,
can also use the standard frequentist p-values:
`pExcl` and `pDisc` instead of
`CLExcl` and `CLDisc`
by simply setting the boolean input parameters`CLExclbool` (for exclusion)
and `CLDiscbool` (for discovery) to `False`.

#### Observed

##### Single-channel

`ZDiscObs(n, bhat, dbhat, s)`:
computes the observed significance for discovery for a single-channel counting experiment

`ZExclObs(n, bhat, dbhat, s)`:
computes the observed significance for exclusion for a single-channel counting experiment

`nDiscObs(bhat, dbhat, Z)`:
computes the number of observed events needed for `Z` sigma discovery using the standard frequentist method by default. Note that the standard p-value for discovery does not depend on `s`.

`sExclObs(n, bhat, dbhat, CL)`:
computes the observed upper limit on the signal at confidence level `CL`

##### Multi-channel

`ZDiscObsN(narray, bhatarray, dbhatarray, sarray)`:
computes the observed significance for discovery for a counting experiment with multiple independent search channels

`ZExclObsN(narray, bhatarray, dbhatarray, sarray)`:
computes the observed significance for exclusion for a counting experiment with multiple independent search channels

#### Expected

##### Single-channel

`ZDiscExp(s, bhat, dbhat)`:
computes the exact Asimov expected significance for discovery for a single-channel counting experiment

`ZExclExp(s, bhat, dbhat)`:
computes the exact Asimov expected significance for exclusion for a single-channel counting experiment

`sDiscExp(bhat, dbhat, Z)`:
computes the signal needed for an expected `Z` sigma discovery

`sExclExp(bhat, dbhat, CL)`:
computes the expected upper limit on the signal at confidence level `CL`

(Note: the functions `ZDiscExp` and `ZExclExp` also contain a
boolean input parameter `asimov_only`that is set to `True` by default such that only the
*exact Asimov significance* is returned. Instead if `asimov_only=False`, these functions return
an ordered list of  
{

* exact Asimov significance,
* mean expected significance,
* mean expected significance with positive significances only,
* median expected significance ( = 50% quantile of `n`) for the default setting of the input parameter `quantile=0.5`,
* significance of mean of p-values,
* probaility of obtaining a significance greater than some input value of significance `Zcriteria`

}  
 
Ref. [1] advocated for the *exact Asimov significance* as the standard significance measure
for projected exclusions and discovery sensitivities in counting experiments.)

##### Multi-channel

`ZDiscExpN(sarray, bhatarray, dbhatarray)`:
computes the exact Asimov expected significance for discovery for a counting experiment with multiple independent search channels 

`ZExclExpN(sarray, bhatarray, dbhatarray)`:
computes the exact Asimov expected significance for exclusion for a counting experiment with multiple independent search channels

###### Other functions

`ZDiscExpCCGV(s, b, db)`:
computes the Asimov approximation to the median expected significance for the discovery case,
obtained in Refs. [3, 4], based on a profile likelihood ratio method

`ZExclExpKM(s, b, db)`:
computes the Asimov approximation to the median expected significance for the exclusion case,
obtained in Ref. [5], based on a profile likelihood ratio method

(Note: The true background mean `b` that appears in the above two functions is in principle unknown if background uncertainty `db` > 0.)

### `Zstats.ProtonDecay` module:

In this module, the tools in the Zstats package are extended to apply to current
and future proton decay experiments. In particular, tools are provided to compute the confidence
level of excluding a proton partial lifetime based on data already taken (e.g. Super-Kamiokande),
along with the discovery and exclusion reach estimates of proton partial lifetimes at (planned)
future proton decay experiments (such as DUNE, JUNO, Hyper-Kamiokande, and THEIA), possibly with
multiple independent channels and with uncertain backgrounds and signal efficiencies.

In addition to input parameters defined above, we also need to define the following for this module:

* `proton_tau`: proton partial lifetime in years
* `Np`: number of protons per kiloton of detector material
* `Nkton`: number of kiltons of detector material
* `exposure`: exposure (= `runtime * Nkton`) of the experiment in kiloton-years
* `eps`: signal selection efficiency (0 <= `eps` <= 1)
* `deps`: uncertainty in the signal selecttion efficiency

Once again, for multi-channel experiments, the some of the above inputs
(`n`, `bhat`, `dbhat`, `eps`, `deps`, `exposure`)
should be
*python lists* or *numpy arrays* which will be referred to as
(`narray`, `bhatarray`, `dbhatarray`, `epsarray`, `depsarray`, `exposurearray`)
here.

The estimates for the backgrounds and the signal selection efficiencies in a specific proton decay
mode have been obtained by the DUNE, JUNO, and THEIA (detector concept) collaborations by modeling
the experiments as single-channel counting experiments with background and signal efficiency
perfectly known,
whereas Hyper-Kamiokande searches for proton decay are modeled as multi-channel counting experiments with uncertain backgrounds and signal efficiencies
based on the signal regions and search strategies used at Super-Kamiokande.

(The Python help function used on the functions for single-channel proton decay experiments will also print out some estimates of background and signal efficiency estimates for various proton decay modes at
DUNE, JUNO, and THEIA.  
And the help function on
`Lifetime*ExpN`(`ConfidenceLevelExclObsN` ) will print out the estimates
for `bhatarray`, `dbhatarray`, `epsarray`, `depsarray` (along with `narray`)
for various proton decay modes at Hyper-Kamiokande (Super-Kamiokande).)

#### Single-channel proton decay experiments (with `dbhat` = `deps` = 0)

`LifetimeDiscExp_Unc0(b, eps, exposure, Np, Z)`:
computes the discovery reach for proton partial lifetime in years at a chosen significance `Z`
at (planned) future proton decay experiments that treat the estimated background and
signal efficiency to be perfectly known while modeling the experiment as a
single-channel counting experiment

`LifetimeExclExp_Unc0(b, eps, exposure, Np, CL)`:
computes the exclusion reach for proton partial lifetime in years at a chosen confidence
level `CL` at (planned) future proton decay experiments that treat the estimated
background and signal efficiency to be perfectly known while modeling the experiment as
a single-channel counting experiment

`RuntimeDiscExp_Unc0(proton_tau, b, eps, Nkton, Np, Z)`:
computes the runtime in years needed for an expected discovery of proton partial lifetime in years
at a chosen significance `Z` at (planned) future proton decay experiments that treat the estimated
background and signal efficiency to be perfectly known while modeling the experiment as a
single-channel counting experiment

`RuntimeExclExp_Unc0(proton_tau, b, eps, Nkton, Np, CL)`:
computes the runtime in years needed for an expected exclusion of proton partial lifetime in years
at a chosen confidence level `CL` at (planned) future proton decay experiments that treat the estimated
background and signal efficiency to be perfectly known while modeling the experiment as a
single-channel counting experiment

#### Multi-channel proton decay experiments

`LifetimeDiscExpN(bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, Z)`:
computes the discovery reach for proton partial lifetime in years at a chosen
significance `Z` at (planned) future proton decay experiments with backgrounds
and signal efficiencies along with their corresponding uncertainties estimated
while modeling the experiment as a multi-channel counting experiment

`LifetimeExclExpN(bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np, CL)`:
computes the exclusion reach for proton partial lifetime in years at a chosen
significance `CL` at (planned) future proton decay experiments with backgrounds
and signal efficiencies along with their corresponding uncertainties estimated
while modeling the experiment as a multi-channel counting experiment

`ConfidenceLevelExclObsN(proton_tau, bhatarray, dbhatarray, epsarray, depsarray, exposurearray, Np)`:
computes the confidence level for (observed) exclusion of proton partial lifetime in years
based on the data collected at a multi-channel proton decay experiment (such as Super-Kamiokande)

### `Zstats.FC` module:

This module contains tools to compute the Feldman-Cousins upper limits on signal along the
experimental sensitivity, both detailed and advocated in Ref. [6]. This method is only applicable
for experiments with single-channel and if the background is known perfectly.

#### Observed

`UpperLimit(n, b, CL)`:
computes the Feldman-Cousins upper limit at confidence level `CL`

#### Expected

`Sensitivity(b, CL)`:
computes the Feldman-Cousins experimental sensitivity at confidence level `CL`

## Examples

To illustrate the usage of the code, this repository also has:

* short programs at `examples/2009.07249/makefig*.py` that produce the data in each of the 10 figures in Reference [1]
* some code snippets in a Python notebook at `examples/2210.07735/makefigs.ipynb` that generate the data in each of the figures in Reference [2] 

<!---
**Note: Each computation with functions `Zdisc`, `Zexcl`, particularly when `asimov_only` is set to `False`, takes a lot more time for the uncertain background case.**
-->

## References

[1]: P. N. Bhattiprolu, S. P. Martin, J. D. Wells, “Criteria for projected discovery and exclusion sensitivities of counting experiments,” arXiv: 2009.07249 [physics.data-an].

[2]: P. N. Bhattiprolu, S. P. Martin, J. D. Wells, “Statistical significances and projections for proton decay experiments,” arXiv: 2210.07735 [hep-ph].

[3]: G. Cowan, K. Cranmer, E. Gross and O. Vitells, “Asymptotic formulae for likelihood-based tests of new physics,” Eur. Phys. J. C 71, 1554 (2011) [arXiv:1007.1727 [physics.data-an]].

[4]: G. Cowan, “Two developments in tests for discovery: use of weighted Monte Carlo events and an improved measure”, Progress on Statistical Issues in Searches,” SLAC, June 4 - 6, 2012.

[5]: N. Kumar and S. P. Martin, “Vectorlike Leptons at the Large Hadron Collider,” Phys. Rev. D 92, no.11, 115018 (2015) [arXiv:1510.03456 [hep-ph]].

[6] G. J. Feldman and R. D. Cousins, “A Unified approach to the classical statistical analysis of
small signals,” Phys. Rev. D 57, 3873-3889 (1998) [arXiv:physics/9711021 [physics.data-an]].
