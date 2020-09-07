# Zstats package
### Criteria for projected discovery and exclusion sensitivities of counting experiments (https://arxiv.org/abs/2009.*****)
### P. N. Bhattiprolu, S. P. Martin, J. D. Wells

Python tools to compute the **exact Asimov significance** <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{A}"> advocated in Reference [1] as the standard significance measure for projected exclusions and discovery sensitivities in counting experiments, both when the background is known and when it is subject to some uncertainty.

**(The documentation of all functions in this package can also be accessed using the python help funtion)**

## Installation

### From GitHub
From command line, run

    $ python -m pip install "git+https://github.com/prudhvibhattiprolu/Zstats.git#egg=Zstats"

## Significance <img src="https://latex.codecogs.com/gif.latex?(Z)"> measures in Reference [1]

<img src="https://latex.codecogs.com/gif.latex?Z^\textrm{A}"> : Exact Asimov significance

<img src="https://latex.codecogs.com/gif.latex?Z^\textrm{mean}"> : Arithmetic mean of <img src="https://latex.codecogs.com/gif.latex?Z"> values

<img src="https://latex.codecogs.com/gif.latex?Z^\textrm{mean}(Z>0)"> : Arithmetic mean of positive <img src="https://latex.codecogs.com/gif.latex?Z"> values only 

<img src="https://latex.codecogs.com/gif.latex?Z^\textrm{med}"> : Median expected significance

<img src="https://latex.codecogs.com/gif.latex?Z^{p\,\textrm{mean}}"> : Significance obtained from the arithmetic mean of the <img src="https://latex.codecogs.com/gif.latex?p">-values

<img src="https://latex.codecogs.com/gif.latex?P(Z_\textrm{disc}>Z)"> or <img src="https://latex.codecogs.com/gif.latex?P(Z_\textrm{excl}>Z)"> : Probability of obtaining <img src="https://latex.codecogs.com/gif.latex?Z_\textrm{disc}"> or <img src="https://latex.codecogs.com/gif.latex?Z_\textrm{excl}"> greater than a certain <img src="https://latex.codecogs.com/gif.latex?Z"> of choice, in a large number of pseudo-experiments simulated for discovery or exclusion case
 
<img src="https://latex.codecogs.com/gif.latex?Z^\textrm{Cowan}_\textrm{disc}"> : Asimov approximation to <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{med}"> for the discovery case, obtained in Ref. [2], based on a likelihood ratio method

<img src="https://latex.codecogs.com/gif.latex?Z^\textrm{KM}_\textrm{excl}"> : Asimov approximation to <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{med}"> for the exclusion case, obtained in Ref. [3], based on a likelihood ratio method

## Expected significances with Zstats package

### Parameters

`s` : Poisson mean of the signal.

`bhat` : Poisson mean of the background in the signal region.

`dbhat` : Poisson uncertainty of the background in the signal region. Set to `0` by default.

`asimov_only` :  This boolean parameter is set to `True` by default and the functions `Zdisc`, `Zexcl` only return the **exact Asimov significance** <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{A}">. Instead, if this option is set to `False`, the functions `Zdisc`, `Zexcl` return an ordered list of <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{A}">, <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{mean}">, <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{mean}(Z>0)">, <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{med}">, <img src="https://latex.codecogs.com/gif.latex?Z^{p\,\textrm{mean}}">, <img src="https://latex.codecogs.com/gif.latex?P(Z_\textrm{disc}>\textrm{Zcriteria})"> or <img src="https://latex.codecogs.com/gif.latex?P(Z_\textrm{excl}>\textrm{Zcriteria})">.

`Zcriteria` : This parameter is used only when `asimov_only` is set to `False` to compute <img src="https://latex.codecogs.com/gif.latex?P(Z_\textrm{disc}>\textrm{Zcriteria})"> or <img src="https://latex.codecogs.com/gif.latex?P(Z_\textrm{excl}>\textrm{Zcriteria})">. By default `Zcriteria=5` in the function `Zdisc` , and `Zcriteria=1.645` in the function `Zexcl`.

`more_info` : This boolean parameter is set to `False` by default. Set this to `True` to show more information about the computations.


### Functions

`Zdisc(s, bhat, dbhat=0, asimov_only=True, Zcriteria=5.0, more_info=False)` : Computes the expected significance for discovery given `s`, `b`, and `dbhat`. 

`Zexcl(s, bhat, dbhat=0, asimov_only=True, Zcriteria=1.645, more_info=False)` : Computes the expected significance for exclusion given `s`, `b`, and `dbhat`.

**More information about each of the input parameter is given above**

### Other functions

`ZdiscAsimovCowan(s, b, db=0)` : Computes <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{Cowan}_\textrm{disc}"> given `s`, true background mean `b`, and background uncertainty `db`.

`ZexclAsimovKM(s, b, db=0)` :  Computes <img src="https://latex.codecogs.com/gif.latex?Z^\textrm{KM}_\textrm{excl}"> given `s`, true background mean `b`, and background uncertainty `db`.

**Note: The true background mean `b` is, in principle, unknown if background uncertainty `db > 0`.**

## References

[1]: P. N. Bhattiprolu, S. P. Martin, J. D. Wells, "Criteria for projected discovery and exclusion sensitivities of counting experiments," arXiv: 2009.***** [physics.data-an].

[2]: G. Cowan, “Two developments in tests for discovery: use of weighted Monte Carlo events and an improved measure”, Progress on Statistical Issues in Searches,” SLAC, June 4 - 6, 2012.

[3]: N. Kumar and S. P. Martin, “Vectorlike Leptons at the Large Hadron Collider,” Phys. Rev. D 92, no.11, 115018 (2015) [arXiv:1510.03456 [hep-ph]].
