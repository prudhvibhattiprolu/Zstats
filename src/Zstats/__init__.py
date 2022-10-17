#Tools to compute the exclusion and discovery significances at current and
#future single- and multi-channel counting experiments based on
#arXiv:2009.07249 [physics.data-an] and arXiv:2210.07735 [physics.data-an].
#Also included are tools to compute the upper limit on signal for exclusion
#and the signal needed for discovery using various methods.
#As an application, these tools are applied to study the statistical
#significance of proton decay experiments.

#arXiv:2009.07249 [physics.data-an] advocated for the 
#*exact Asimov significance* as the standard significance measure for
#projected exclusions and discovery sensitivities in counting experiments.
#arXiv:2210.07735 [physics.data-an] generalized
#arXiv:2009.07249 [physics.data-an] to include Bayesian and
#modified frequentist statistics, and multi-channel counting experiments,
#with application to proton decay at current and future neutrino detectors.

#Import as a python package (if Zstats package is in path):

## >>> import Zstats as Z
### Now, e.g., to compute the expected discovery significance:
## >>> Z.ZDiscExp(5, 1)
### for more information about each function/input parameters
### use the Python help function:
## >>> help(Z.ZDiscExp) #or
## >>> print(Z.ZDiscExp.__doc__)

#Python packages used
import numpy as np
import scipy
import scipy.special as sc
from scipy.integrate import quad
import mpmath as mp #mpmath can be made faster by installing *gmpy2*
import itertools as it

import sys

sys.path.append('../')
#Import local packages
import FC #C++ tools (with python bindings) to compute Feldman-Cousins upper limits on signal and experimental sensitivity as detailed in arXiv:physics/9711021
import ProtonDecay #Python + C++ tools (with python bindings) to compute current limits and future exclusion and discovery reaches for proton partial lifetimes at (single/multi)-channel proton decay experiments 
del sys.path[-1]

#Convert a one-sided p-value to a significance Z
def Zfromp(p):
    '''
    Converts one-sided p-value to significance Z using:

        Z = Sqrt[2] InverseErfc[2 p]
    '''
    return np.sqrt(2)*sc.erfcinv(2*p)

#Convert a significance to a one-sided p-value
def pfromZ(Z):
    '''
    Converts significance Z to one-sided p-value using:

        p = 0.5 Erfc[Z/Sqrt[2]]
    '''
    return 0.5*sc.erfc(Z/np.sqrt(2))

############### Single-channel counting experiments ###############

#Poisson probability with Bayesian marginalization for the uncertain background case
def DeltaP(n, m, tau, s, b):
    '''
    Inputs
    ------

    n: number of observed events in the signal (on) region

    m: number of observed events in the supposed
       background-only (off) region

    tau: ratio of background means in off and on regions

    s: Poisson mean of the signal

    b: Poisson mean of the background in the on region

    Note: In the case of uncertain background,
    *m* and *tau* should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    and
    dbhat: background uncertainty in the on region

    Returns
    -------

    Known background case (*m* and *tau* are both infinite):
    returns the Poisson probability of observing n-events given a mean = s + b

    Uncertain background case (*m* and *tau* are finite):
    returns the Poisson probability (with Bayesian marginalization) of
    observing *n* events in the on region, given *m* and *tau* in the
    on-off problem, and number of signal events *s*
    (here, *b* is a dummy variable and should not affect the result)

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2009.07249 [physics.data-an]]
    '''
    n = float(n)
    m = float(m)
    tau=float(tau)
    b=float(b)
    s=float(s)
    if m == scipy.inf and tau == scipy.inf:
        result = float(mp.exp(-s-b)*mp.power(s+b,n)/mp.factorial(n))
    else:
        if s == 0 :
            result = float((mp.gamma(n+m+1)/(mp.gamma(n+1)*mp.gamma(m+1)))*(mp.power(tau,m+1)/mp.power(1+tau,n+m+1)))
        else:
            if n.is_integer():
                result = float((mp.exp(-s)*mp.power(tau,m+1)/mp.gamma(m+1))*sum([(mp.power(s,k)/mp.factorial(n-k))*(mp.gamma(n+m+1-k)/(mp.factorial(k)*mp.power(1+tau,n+m+1-k))) for k in mp.arange(n+1)]))
            else:
                result = float((mp.exp(-s)*mp.power(tau,m+1)/(mp.gamma(m+1)*mp.gamma(n+1)))*mp.quad(lambda bb: mp.power(bb,m)*mp.exp(-bb*(1+tau))*mp.power(s+bb,n),[0,mp.inf]))
    return result

# Maximum Poisson count with a significant probability
def _MaxPoissonCount(s, b):
    '''
    Inputs
    ------

    s: Poisson mean of the signal

    b: Poisson mean of the background in the on region

    Returns
    -------

    Maximum Poisson count with a significant probability
    '''
    totalprob_max = 0.9999999999
    totalprob = 0.0

    n=0
    while totalprob <= totalprob_max:
        prob = DeltaP(n, scipy.inf, scipy.inf, s, b)
        totalprob = totalprob + prob
        n = n + 1
    return int(n)

# Discovery - Frequentist
def pDisc(n, m, tau, b, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal (on) region

    m: number of observed events in the supposed
       background-only (off) region

    tau: ratio of background means in off and on regions

    b(hat): Poisson mean of the background (in the on region)

    more_info: Set to False (or 0) by default. Set this to True to show
               more information about the computation of p-values

    Note: In the case of uncertain background,
    m and tau should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    and
    dbhat: background uncertainty in the on region

    Returns
    -------

    p-value for discovery for:

    Known background case - when *m* and *tau* are infinite

    Uncertain background case - when *m* and *tau* are finite
    (here *b* is a dummy variable and should not affect the result)

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2009.07249 [physics.data-an]]
    '''
    n = float(n)
    m = float(m)
    tau=float(tau)
    b=float(b)
    if more_info:
        print('Inputs:')
        print('n = %s' %(n))
        print('m = %s' %(m))
        print('tau = %s' %(tau))
        print('b = %s\n' %(b))
        print('s = 0')
    if m == scipy.inf and tau == scipy.inf:
        #p-value for discovery when background is known exactly.
        pBi = sc.gammainc(n,b) if n != 0 else 1.0 #has round-off errors for large n and small b, for example, when (n,b)=(>=118,0.1)
        if more_info:
            print('p-value for discovery case when background is known exaclty:\n')
            print('\tpDisc = %s\n' %(pBi))
    else:
        #p-value for discovery when background is not known exactly.
        w = 1/(1+tau)
        pBi = sc.betainc(n,m+1,w) if n != 0 else 1.0
        if more_info:
            print('p-value for discovery case when background is not known exaclty:\n')
            print('\tpDisc = %s\n' %(pBi))
    return pBi

##p-value/power for the discovery case used in modified frequentist approach instead of the standard p-value
#def CLDisc_modifiedfrequentist(n, m, tau, b, s, more_info=False):
#    n=float(n)
#    m=float(m)
#    tau=float(tau)
#    b=float(b)
#    s=float(s)
#    if b != 0:
#        if m == scipy.inf and tau == scipy.inf:
#            result=sc.gammainc(n,b)/sc.gammainc(n,s+b) if n != 0 else 1.0
#        else:
#            num=pDisc(n,m,tau,b, more_info)
#            den=1 - pExcl(n-1,m,tau,b,s, more_info)
#            result=num/den if den != 0 else np.nan
#    else:
#        if (n != 0) and (s != 0):
#            result=0.0
#        else:
#            result=1.0
#    return result

#Discovery - Bayesian (using Bayes Factors)
def CLDisc(n, m, tau, b, s):
    '''
    Inputs
    ------

    n: number of observed events in the signal (on) region

    m: number of observed events in the supposed
       background-only (off) region

    tau: ratio of background means in off and on regions

    b(hat): Poisson mean of the background (in the on region)

    more_info: Set to False (or 0) by default. Set this to True to show
               more information about the computation of p-values

    Note: In the case of uncertain background,
    m and tau should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    and
    dbhat: background uncertainty in the on region

    Returns
    -------

    Bayesian statistic CLExcl for discovery obtained using Bayes factors
    that can be used in place of p-value for:

    Known background case - when *m* and *tau* are infinite

    Uncertain background case - when *m* and *tau* are finite
    (here *b* is a dummy variable and should not affect the result)

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    n=float(n)
    m=float(m)
    tau=float(tau)
    b=float(b)
    s=float(s)
    return DeltaP(n, m, tau, 0, b)/DeltaP(n, m, tau, s, b)

#Observed discovery significances using the standard frequentist pDisc or
#a Bayesian statistic CLDisc
def ZDiscObs(n, bhat, dbhat, s, CLDiscbool=True, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal region

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region

    s: Poisson mean of the signal. This is used only when
       *CLDiscbool* is set to *True*

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                discovery instead of a Bayesian statistic CLDisc

    more_info: Set to *False* by default. Set this to *True* to show
               more information about the computation of significance

    Returns
    -------

    This function returns the observed significance for discovery
    given *n*, *b*, *dbhat*, and *s*

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2009.07249 [physics.data-an]]
    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Trading (bhat, dbhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    if more_info:
        print('Inputs:')
        print('n = %s' %(n))
        print('b hat = %s' %(bhat))
        print('dbhat = %s\n' %(dbhat))
        print('s = %s' %(s))
        print('Trading (bhat, dbhat) for (m, tau):')
        print('m = %s' %(m))
        print('tau = %s\n' %(tau))
        if CLDiscbool: print('\nTaking a conservative approach of using a Bayesian statistic CLDisc instead of the standard p-value:\n')
    #Computing the observed significance
    ZObs = Zfromp(pDisc(n,m,tau,bhat)) if not CLDiscbool else Zfromp(CLDisc(n,m,tau,bhat,s))
    if more_info:
        print('\tObserved Z = %s\n' %(ZObs))
    return ZObs

#Expected discovery significances using the standard frequentist pDisc or
#a Bayesian statistic CLDisc
def ZDiscExp(s, bhat, dbhat=0, asimov_only=True, CLDiscbool=True, quantile=0.5, Zcriteria=5.0, more_info=False):
    '''
    Inputs
    ------

    s: Poisson mean of the signal

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region.
           Set to *0* by default

    asimov_only: By default this is set to *True* and the function
                 only returns the exact Asimov discovery significance.
                 Instead if this option is set to *False*, the function
                 returns a list of
                 {
                  Exact Asimov significance,
                  Mean Z with Z != -Infinity,
                  Mean Z with Z > 0,
                  Median Z = Z(100% times *quantile* of *n*),
                  Z obtained from Mean p-value,
                  Probability of obtaining ZDisc > Zcriteria
                 }

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                discovery instead of a Bayesian statistic CLDisc

    quantile: This option is used only when *asimov_only* is set to *False*
              where the fourth element of the list returned by the function is
              Z(50% quantile of *n*) = Median Z for the default setting of
              *quantile*=0.5. For any other setting of *quantile*,
              Z(100% times *quantile* of *n*) is returned

    Zcriteria: Set to 5.0 sigma by default. This is used only when
               *asimov_only* is set to *False* to compute
               the probability of obtaining ZDisc > Zcriteria

    more_info: Set to *False* by default. Set this to *True* to show
               more information about the computation of significance(s)

    Returns
    -------

    This function returns the expected significance(s) for discovery
    given *s*, *b* and *dbhat*

    If the option *asimov_only* is set to *True* (default case):
    only returns the exact Asimov discovery significance

    If the option *asimov_only* is set to *False*:
    returns a list of
    {
        Exact Asimov significance,

        Mean Z with Z != -Infinity,

        Mean Z with Z > 0,

        Z(100% times *quantile* of number of pseudo-experiments *n*),

        Z obtained from Mean p-value,

        Probability of obtaining ZDisc > Zcriteria
    }

    To use a more conservative Bayesian statistic CLDisc instead of
    the standard frequentist pDisc, set *CLDiscbool* to *True*

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2009.07249 [physics.data-an]]
    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Trading (bhat, dbhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    #Mean n
    if bhat != 0:
        nmean = s + bhat + dbhat**2/bhat
    else:
        nmean = s
    if more_info:
        print('Inputs:')
        print('s = %s' %(s))
        print('b hat = %s' %(bhat))
        print('dbhat = %s\n' %(dbhat))
        print('Trading (bhat, dbhat) for (m, tau):')
        print('m = %s' %(m))
        print('tau = %s\n' %(tau))
        if CLDiscbool: print('\nTaking a conservative approach of using a Bayesian statistic CLDisc instead of the standard p-value:\n')
    #Computing exact Asimov significance
    Zasimov = Zfromp(pDisc(nmean,m,tau,bhat)) if not CLDiscbool else Zfromp(CLDisc(nmean,m,tau,bhat,s))
    if asimov_only:
        if more_info:
            print('Mean n = %s\n' %(nmean))
            print('Recommended: (The exact Asimov discovery significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
        return Zasimov
    else:
        nobs = -1
        probtotal = 0
        foundQuantile = False
        Zmeanlist = []
        pmeanlist = []
        pZlist = []
        if more_info:
            print('Listing all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability:\n')
            print('%-10s %-25s %-25s %-25s %-10s' %('n', 'Probability', 'Cumulative Probability', 'pvalue', 'Significance (ZDisc)'))
        #Considering all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability
        while probtotal < 0.99999999:
            nobs = nobs + 1
            probval = DeltaP(nobs,m,tau,s,bhat)
            pval = pDisc(nobs,m,tau,bhat) if not CLDiscbool else CLDisc(nobs,m,tau,bhat,s)
            Zval = Zfromp(pval)
            Zmeanlist.append(probval*Zval)
            pmeanlist.append(probval*pval)
            if Zval > Zcriteria: pZlist.append(probval)
            probtotal = probtotal + probval
            if more_info: print('%-10s %-25s %-25s %-25s %-10s' %(nobs, probval, probtotal, pval, Zval))
            if probtotal > quantile and not foundQuantile:
                #Finding Median Z
                nquantile = nobs
                Zquantile = Zval
                foundQuantile = True
        #Computing Mean Z if Z != -Infinity
        Zmean_not_neginf = sum([i for i in Zmeanlist if not np.isneginf(i)])
        #Computing Mean Z if Z > 0
        Zmean_not_neg = sum([i for i in Zmeanlist if i > 0])
        #Computing Mean pvalue
        pmean = sum(pmeanlist)
        #Converting pmean to significance
        Zpmean = Zfromp(pmean)
        #Computing probability of obtaining a significance  Z > Zcriteria
        pZcriteria = sum(pZlist)
        if more_info:
            print('\nMean n = %s' %(nmean))
            if quantile == 0.5:
                print('Median n = 50%% quantile of the number of pseudo-experiments n = %s\n' %(nquantile))
            else:
                print('%s%% quantile of the number of pseudo-experiments n = %s\n' %(quantile*100, nquantile))
            print('Discovery significances:\n')
            print('Recommended: (The exact Asimov significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
            print('Reasonable:\n')
            print('\tMean Z (Z != -Infinity) = %s\n' %(Zmean_not_neginf))
            print('\tMean Z (Z > 0) = %s\n' %(Zmean_not_neg))
            print('Not Recommended:\n')
            if quantile == 0.5:
                print('\tMedian Z = Z(50%% quantile of the number of pseudo-experiments n) = %s\n' %(Zquantile))
            else:
                print('\tZ(%s%% quantile of the number of pseudo-experiments n) = %s\n' %(quantile*100, Zquantile))
            print('\tZ of Mean p-value = %s\n' %(Zpmean))
            print('Probability of obtaining a significance (ZDisc) greater than Zcriteria = %s is %s\n' %(Zcriteria, pZcriteria))
        return [Zasimov, Zmean_not_neginf, Zmean_not_neg, Zquantile, Zpmean, pZcriteria]

#Number of observed events needed for discovery obtained using
#standard frequentist pDisc or a Bayesian statistic CLDisc
def nDiscObs(bhat, dbhat=0, Z=3, CLDiscbool=False, s=1, integer_only=True, more_info=False):
    '''
    Inputs
    ------

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region set to 0 by default

    Z: Significance needed for discovery. Set to *3* sigma by default

    CLDiscbool: Set to *True* to use a Bayesian statistic CLDisc instead of
                the standard frequentist p-value for discovery

    s: Poisson mean of the signal set to 1 by default.
       This is used only when *CLDiscbool* is set to *True*

    integer_only: Set to true by default to return the result
                  rounded-up to the next integer

    more_info: Set to *False* by default. Set this to *True* to show
               more information about the computation of significance(s)

    Returns
    -------

    The number of observed events needed for *Z* sigma discovery using the
    standard frequentist method by default
    '''
    #Trading (bhat, dbhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    #Equation to solve for n
    if not CLDiscbool:
        #p-value (standard frequentist)
        solve = lambda n: pDisc(n, m, tau, bhat) - pfromZ(Z)
    else:
        #CLDisc (Bayesian)
        solve = lambda n: CLDisc(n, m, tau, bhat, s) - pfromZ(Z)
    try:
        #Initial guess for mpmath to find root
        if dbhat == 0:
            n_guess = bhat
        else:
            #To speed up root finding for the uncertain background case, use the result for the known background case as the initial guess
            n_guess = nDiscObs(bhat, dbhat=0, Z=Z, CLDiscbool=CLDiscbool, s=s, integer_only=False)
        #Solve for root
        result = mp.findroot(solve, n_guess, verbose=more_info)
        #Sometimes the result could be NaN, so raise an error if that is the case
        if mp.isnan(result): raise TypeError
    except:
        #In case of errors while solving for root in the try block, the except block is executed
        try:
            #n_guess = 1
            n_guess = nDiscObs(bhat, dbhat=0, Z=Z, CLDiscbool=0, s=s, integer_only=False)
            result = mp.findroot(solve, n_guess, verbose=more_info)
            if mp.isnan(result): raise TypeError
        except:
            if more_info: print("Solution doesn't exist or has numerical instabilities.")
            result=mp.nan
            pass
    if (integer_only) and (not mp.isnan(result)): result = np.ceil(result)
    return result

# Signal needed for discovery obtained using the standard frequentist pDisc or
#a Bayesian statistic CLDisc
def _sDisc(bhat, dbhat=0, expected_limit=True, CLDiscbool=True, n=0, Z=3, more_info=False):
    #Trading (bhat, dbhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    if expected_limit:
        #Setting n equal to its expected mean
        if bhat != 0:
            nmean = lambda s: s + bhat + dbhat**2/bhat
        else:
            nmean = lambda s: s
    #Significance Z expressed as p-value
    pcriteria = pfromZ(Z)
    #Compute the number of observed events needed for discovery if expected_limit is set to *False*
    if not expected_limit: n_limit = nDiscObs(bhat=bhat, dbhat=dbhat, Z=Z, CLDiscbool=False, s=1, integer_only=True, more_info=False)
    #Equation to solve for the signal needed for discovery given bhat, dbhat, Z
    if CLDiscbool:
        #CLDisc (Bayesian)
        if expected_limit:
            solve = lambda s: CLDisc(nmean(s), m, tau, bhat, s) - pcriteria
        else:
            solve = lambda s: CLDisc(n, m, tau, bhat, s) - pcriteria
            if (n < n_limit) and (more_info):
                print("No solution unless the number of observed events are at least %s for %s sigma discovery." %(n_limit, Z))
    else:
        #p-value (standard frequentist)
        if expected_limit:
            solve = lambda s: pDisc(nmean(s), m, tau, bhat) - pcriteria
        else:
            solve = lambda s: pDisc(n, m, tau, bhat) - pcriteria
            if n >= n_limit:
                sysexit_text = "Also, according to pure frequentist approach, for the given expected background, the number of observed events should be atleast %s and so this could be a %s sigma discovery." %(n_limit, Z)
            else:
                sysexit_text = "Also, according to pure frequentist approach, for the given expected background, the number of observed events should be atleast %s to claim a %s sigma discovery." %(n_limit, Z)
            raise SystemExit("For a fixed *n*, *pDisc* does not depend on *s*, so the solution for *s* doesn't exist. \nInstead, the expected signal for discovery is %s, obtained by replacing *n* by its expected mean.\n" %(_sDisc(bhat, dbhat, expected_limit=True, CLDiscbool=CLDiscbool, n=n, Z=Z, more_info=False)) + sysexit_text)
    try:
        #Initial guess for mpmath to find root
        if dbhat == 0:
            #s_guess = nmean(0) if expected_limit else n-bhat
            s_guess = np.sqrt(-bhat*np.log(pcriteria)) if expected_limit else n-bhat
        else:
            #To speed up root finding for the uncertain background case, use the result for the known background case as the initial guess
            s_guess = _sDisc(bhat, 0, expected_limit, CLDiscbool, n, Z, more_info=False)
        #Solve for root
        result = mp.findroot(solve, s_guess, verbose=more_info)
        #Sometimes the result could be NaN, so raise an error if that is the case
        if mp.isnan(result): raise TypeError
    except:
        #In case of errors while solving for root in the try block, the except block is executed
        try:
            if expected_limit:
                s_guess = 1
                result = mp.findroot(solve, s_guess, verbose=more_info)
            else:
                s_guess = (0, n-bhat) if dbhat==0 else (_sDisc(bhat, 0, expected_limit, CLDiscbool, n, Z, more_info=False), n-(m+1)/tau)
                rootsolvetol = None if dbhat==0 else 1e-9
                result = mp.findroot(solve, s_guess, verbose=more_info, solver='ridder', maxsteps=200, tol=rootsolvetol)
            if mp.isnan(result): raise TypeError
        except:
            if more_info: print("Solution doesn't exist or has numerical instabilities.")
            result=mp.nan
    return float(result)

def sDiscExp(bhat, dbhat=0, Z=3, CLDiscbool=True, more_info=False):
    '''
    Inputs
    ------

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region.
           Set to *0* by default

    Z: Required significance for discovery. Set to *3* sigma by default

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                discovery instead of a Bayesian statistic CLDisc

    more_info: Set to *False* by default. Set this to *True* to show more
               information about the computation

    Returns
    -------

    The signal needed for an expected *Z* sigma discovery using a more
    conservative Bayesian statistic CLDisc by default

    If *CLDiscbool* is set to *False*, the standard p-value is used instead of
    a more conservative Bayesian statistic CLDisc

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    return _sDisc(bhat, dbhat=0, expected_limit=True, CLDiscbool=CLDiscbool, n=0, Z=Z, more_info=more_info)

def _sDiscObs(n, bhat, dbhat=0, Z=3, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal region

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region.
           Set to *0* by default

    Z: Required significance for discovery. Set to *3* sigma by default

    more_info: Set to *False* by default. Set this to *True* to show more
               information about the computation

    Returns
    -------

    The signal needed for a *Z* sigma discovery for a fixed *n* using a more
    conservative Bayesian statistic CLDisc

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    return _sDisc(bhat, dbhat, expected_limit=False, CLDiscbool=True, n=n, Z=Z, more_info=more_info)

# Discovery - likelihood ratio method (Asimov)
def ZDiscExpCCGV(s, b, db=0):
    '''
    Inputs
    ------

    s: Poisson mean of signal

    b: Poisson mean of background

    db: uncertainty in background

    Returns
    -------

    The expected significance for discovery based on a likelihood ratio method

    References
    ----------

    G. Cowan, K. Cranmer, E. Gross and O. Vitells [arXiv:1007.1727 [physics.data-an]]

    G. Cowan [https://www.pp.rhul.ac.uk/~cowan/stat/cowan_slac_4jun12.pdf]
    '''
    if db == 0:
        ZCCGV = np.sqrt(2 * ((s + b) * np.log(1 + s/b) - s))
    else:
        ZCCGV = np.sqrt(2 * ((s + b) * np.log(((s + b) * (b + db**2))/(b**2 + (s + b) * db**2)) - (b**2/db**2) * np.log(1 + (s * db**2/(b * (b + db**2))))))
    return ZCCGV


# Exclusion - Frequentist
def pExcl(n, m, tau, b, s, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal (on) region

    m: number of observed events in the supposed
       background-only (off) region

    tau: ratio of background means in off and on regions

    b(hat): Poisson mean of the background (in the on region)

    more_info: Set to False (or 0) by default. Set this to True to show
               more information about the computation of p-values

    Note: In the case of uncertain background,
    m and tau should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    and
    dbhat: background uncertainty in the on region

    Returns
    -------

    p-value for exclusion for:

    Known background case - when *m* and *tau* are infinite

    Uncertain background case - when *m* and *tau* are finite
    (here *b* is a dummy variable and should not affect the result)

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2009.07249 [physics.data-an]]
    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    n = float(n)
    m = float(m)
    tau=float(tau)
    b=float(b)
    s=float(s)
    if more_info:
        print('Inputs:')
        print('n = %s' %(n))
        print('m = %s' %(m))
        print('tau = %s' %(tau))
        print('b = %s\n' %(b))
        print('s = %s' %(s))
    if m == scipy.inf and tau == scipy.inf:
        #p-value for exclusion when background is known exactly.
        pBi = sc.gammaincc(n+1,s+b)
        if more_info:
            print('p-value for exclusion case when background is known exactly:\n')
            print('\tpExcl = %s\n' %(pBi))
    else:
        if more_info: print('p-value for exclusion case when background is not known exactly:\n')
        #p-value for exclusion when background is not known exactly.
        if s == 0:
            #pExcl for uncertain background case when s=0 (needed for computing CLs) can be more conveniently computed using incomplete beta functions
            pBi = sc.betainc(m+1,n+1,tau/(1+tau))
            if more_info:
                print('Here, s = 0, so pExcl can be computed using incomplete beta functions:\n')
                print('\tpExcl = %s\n' %(pBi))
        else:
            if n.is_integer():
                pBi = float((mp.power(tau,m+1)/mp.gamma(m+1))*sum([(mp.gamma(k+m+1)/(mp.factorial(k)*mp.power(1+tau,k+m+1)))*sc.gammaincc(n-float(k)+1,s) for k in mp.arange(n+1)]))
                if more_info:
                    print('Here, n is an integer. pExcl evaluated using the sum formula is:\n')
                    print('\tpExcl = %s\n' %(pBi))
            else:
                #only used for Asimov approximations in the exclusion case (if n is not an integer)
                pBilist = []

                #first computing pBisum = pExcl(l,m,tau,s) using sum formula with l = closest integer to n
                pBisum = float((mp.power(tau,m+1)/mp.gamma(m+1))*sum([(mp.gamma(k+m+1)/(mp.factorial(k)*mp.power(1+tau,k+m+1)))*sc.gammaincc(n-float(k)+1,s) for k in mp.arange(n+1)]))#always reliable
                pBilist.append(pBisum)


                #then computing pBiint = pExcl(n,m,tau,s) using the first integral formula (Second line of Eq. 15 in arXiv: 2009.07249)
                #Computing with scipy's quad and gammaincc functions
                pBiint = float((mp.power(tau,m+1)/mp.gamma(m+1))*quad(lambda bb: mp.power(bb,m)*mp.exp(-bb*tau)*sc.gammaincc(n+1,s+bb),0,scipy.inf)[0])#in principle same as pBibyparts, but may have numerical instabilities

                #To compute the integral to arbitrary precision, use mpmath's quad (for numerical integration) and inc. gamma functions
                #To do so, uncomment the below two lines and comment out the line somewhere above, starting with *pBiint*, which uses scipy to do the integral
                #mp.mp.dps = 15 #This sets the global *decimal precision* of *mpmath*. In case of errors, try increasing this number (which in turn increases the computation time).
                #pBiint = float((mp.power(tau,m+1)/mp.gamma(m+1))*mp.quad(lambda bb: mp.power(bb,m)*mp.exp(-bb*tau)*mp.gammainc(n+1,s+bb,regularized=True),[0,mp.inf]))

                if pBiint >= 0 and not (np.isnan(pBiint) or np.isinf(pBiint)):
                    percentdiffint = 100*2*abs((pBisum-pBiint)/(pBisum+pBiint))
                    if percentdiffint < 40:
                        pBilist.append(pBiint)

                #finally computing pBibyparts = pExcl(n,m,tau,s) using the second integral formula (Third/Last line of Eq. 15 in arXiv: 2009.07249)
                #Computing with scipy's quad and gammaincc functions
                pBibyparts = float(sc.gammaincc(n+1,s) - (quad(lambda bb: mp.power(s+bb,n)*mp.exp(-(s+bb))*sc.gammaincc(m+1,tau*bb),0,scipy.inf)[0]/mp.gamma(n+1)))#can be inaccurate when n >> b

                #To compute the integral to arbitrary precision, use mpmath's quad (for numerical integration) and inc. gamma functions
                #To do so, uncomment the below two lines and comment out the line somewhere above, starting with *pBibyparts*, which uses scipy to do the integral
                #mp.mp.dps = 15 #This sets the global *decimal precision* of *mpmath*. In case of errors, try increasing this number (which in turn increases the computation time).
                #pBibyparts = float(mp.gammainc(n+1,s,regularized=True) - (mp.quad(lambda bb: mp.power(s+bb,n)*mp.exp(-(s+bb))*mp.gammainc(m+1,tau*bb,regularized=True),[0,mp.inf])/mp.gamma(n+1)))

                if pBibyparts >= 0 and not (np.isnan(pBibyparts) or np.isinf(pBibyparts)):
                    percentdiffbyparts = 100*2*abs((pBisum-pBibyparts)/(pBisum+pBibyparts))
                    if percentdiffbyparts < 40:
                        pBilist.append(pBibyparts)

                #Find the best solution among {pBisum (always reliable but approximate), pBibyparts, pBiint}
                if len(pBilist) == 3:
                    if mp.fabs(mp.fsub(pBilist[0],pBilist[1],exact=True)) <= mp.fabs(mp.fsub(pBilist[0],pBilist[2],exact=True)):
                        pBi = pBilist[1]
                    else:
                        pBi = pBilist[2]
                elif len(pBilist) == 2:
                    pBi = pBilist[1]
                elif len(pBilist) == 1:
                    pBi = pBilist[0]
                if more_info:
                    print('Here, n is not an integer:\n')
                    print('First, computing psum = pExcl(l,m,tau,s) using a sum formula with l = closest integer to n: psum = %s' %(pBisum))
                    print('(Note: psum is an approximate answer when n != integer, but can be computed more reliably)\n')
                    print('Then, computing pint1 = pExcl(n,m,tau,s) using the first integral formula (second line of Eq. 15 in arXiv: 2009.07249): pint1 = %s\n' %(pBiint))
                    print('Finally, computing pint2 = pExcl(n,m,tau,s) using the second integral formula (third/last line of Eq. 15 in arXiv: 2009.07249): pint2 = %s\n' %(pBibyparts))
                    print('Now, finding the closest of {pint1, pint2} to psum:')
                    print('\tpExcl = %s\n' %(pBi))
                    print('Note: In some cases only one of the integral formulas gives a sensible result. If so, this function returns that.')
                    print('And, in some other rare cases, both pint1 and pint2 are not sensible, where, this function returns psum instead, which is an approximation to the true result.\n')
    return pBi

#Exclusion - modified frequentist

#p-value/p-value with s=0 (= CLs) for the exclusion case used in modified frequentist approach instead of the standard p-value
def CLs(n, m, tau, b, s, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal (on) region

    m: number of observed events in the supposed
       background-only (off) region

    tau: ratio of background means in off and on regions

    b(hat): Poisson mean of the background (in the on region)

    more_info: Set to False (or 0) by default. Set this to True to show
               more information about the computation of p-values

    Note: In the case of uncertain background,
    m and tau should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    and
    dbhat: background uncertainty in the on region

    Returns
    -------

    The CLs (= pExcl/pExcl with s=0) for the exclusion case that is used in
    modified frequetist approach instead of the standard p-value for:

    Known background case - when *m* and *tau* are infinite

    Uncertain background case - when *m* and *tau* are finite
    (here *b* is a dummy variable and should not affect the result)

    Note: CLs is numerically equal to a Bayesian statistic CLExcl for
          for a single-channel counting experiment

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    n=float(n)
    m=float(m)
    tau=float(tau)
    b=float(b)
    s=float(s)
    result=pExcl(n,m,tau,b,s,more_info)/pExcl(n,m,tau,b,0,more_info)
    if np.isnan(result) or np.isinf(result):
        if (m == scipy.inf) and (tau == scipy.inf):
            result = float(mp.gammainc(n+1, s+b)/mp.gammainc(n+1, b))
    return result

#Exclusion - Bayesian
#Bayesian statistic CLExcl that can be used in place of the standard pExcl
def CLExcl(n, m, tau, b, s, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal (on) region

    m: number of observed events in the supposed
       background-only (off) region

    tau: ratio of background means in off and on regions

    b(hat): Poisson mean of the background (in the on region)

    more_info: Set to False (or 0) by default. Set this to True to show
               more information about the computation of p-values

    Note: In the case of uncertain background,
    m and tau should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    and
    dbhat: background uncertainty in the on region

    Returns
    -------

    A Bayesian statistic CLExcl for the exclusion case that can be used
    in the place of the standard p-value for:

    Known background case - when *m* and *tau* are infinite

    Uncertain background case - when *m* and *tau* are finite
    (here *b* is a dummy variable and should not affect the result)

    Note: CLExcl is numerically equal to the modified frequentist CLs for
          for a single-channel counting experiment

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    return CLs(n, m, tau, b, s, more_info)

#Observed exclusion significances using the standard frequentist pExcl or
#a Bayesian statistic CLExcl (= modified frequentist CLs for single-channel)
def ZExclObs(n, bhat, dbhat, s, CLExclbool=True, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal region

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region

    s: Poisson mean of the signal.

    CLExclbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLExcl

    more_info: Set to *False* by default. Set this to *True* to show
               more information about the computation of significance

    Returns
    -------

    This function returns the observed significance for exclusion
    given *n*, *b*, *dbhat*, and *s*

    To use a more conservative Bayesian statistic CLExcl instead of
    the standard frequentist pExcl, set *CLExclbool* to *True*

    Note: CLExcl is numerically equal to the modified frequentist CLs for
          for a single-channel counting experiment

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2009.07249 [physics.data-an]]
    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Trading (bhat, dbhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    if more_info:
        print('Inputs:')
        print('n = %s' %(n))
        print('b hat = %s' %(bhat))
        print('dbhat = %s\n' %(dbhat))
        print('s = %s' %(s))
        print('Trading (bhat, dbhat) for (m, tau):')
        print('m = %s' %(m))
        print('tau = %s\n' %(tau))
        if CLExclbool: print('\nTaking a conservative approach of using a Bayesian statistic CLExcl (= modified frequentist CLs for single-channel) instead of the standard p-value:\n')
    #Computing the observed significance
    ZObs = Zfromp(pExcl(n,m,tau,bhat,s)) if not CLExclbool else Zfromp(CLExcl(n,m,tau,bhat,s))
    if more_info:
        print('\tObserved Z = %s\n' %(ZObs))
    return ZObs

#Expected exclusion significances using the standard frequentist pExcl or
#a Bayesian statistic CLExcl (= modified frequentist CLs for single-channel)
def ZExclExp(s, bhat, dbhat=0, asimov_only=True, CLExclbool=True, quantile=0.5, Zcriteria=1.645, more_info=False):
    '''
    Inputs
    ------

    s: Poisson mean of the signal

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region.
           Set to *0* by default

    asimov_only: By default this is set to *True* and the function
                 only returns the exact Asimov exclusion significance.
                 Instead if this option is set to *False*, the function
                 returns a list of
                 {
                  Exact Asimov significance,
                  Mean Z with Z != -Infinity,
                  Mean Z with Z > 0,
                  Median Z = Z(100% times *quantile* of *n*),
                  Z obtained from Mean p-value,
                  Probability of obtaining ZExcl > Zcriteria
                 }

    CLExclbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLExcl

    quantile: This option is used only when *asimov_only* is set to *False*
              where the fourth element of the list returned by the function is
              Z(50% quantile of *n*) = Median Z for the default setting of
              *quantile*=0.5. For any other setting of *quantile*,
              Z(100% times *quantile* of *n*) is returned

    Zcriteria: Set to 1.645 sigma (<-> 95% CL) by default. This is used only
               when *asimov_only* is set to *False* to compute
               the probability of obtaining ZExcl > Zcriteria

    more_info: Set to *False* by default. Set this to *True* to show
               more information about the computation of significance(s)

    Returns
    -------

    This function returns the expected significance(s) for exclusion
    given *s*, *b* and *dbhat*

    If the option *asimov_only* is set to *True* (default case):
    only returns the exact Asimov exclusion significance

    If the option *asimov_only* is set to *False*:
    returns a list of
    {
        Exact Asimov significance,

        Mean Z with Z != -Infinity,

        Mean Z with Z > 0,

        Z(100% times *quantile* of number of pseudo-experiments *n*),

        Z obtained from Mean p-value,

        Probability of obtaining ZExcl > Zcriteria
    }

    To use a more conservative Bayesian statistic CLExcl instead of
    the standard frequentist pExcl, set *CLExclbool* to *True*

    Note: CLExcl is numerically equal to the modified frequentist CLs for
          for a single-channel counting experiment

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2009.07249 [physics.data-an]]
    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    #Trading (bhat, dbhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    #Mean n
    if bhat != 0:
        nmean = bhat + dbhat**2/bhat
    else:
        nmean = 0
    if more_info:
        print('Inputs:')
        print('s = %s' %(s))
        print('bhat = %s' %(bhat))
        print('dbhat = %s\n' %(dbhat))
        print('Trading (bhat, dbhat) for (m, tau):')
        print('m = %s' %(m))
        print('tau = %s\n' %(tau))
        if CLExclbool: print('\nTaking a conservative approach of using a Bayesian statistic CLExcl (= modified frequentist CLs for single-channel) instead of the standard p-value:\n')
    #Computing exact Asimov significance
    Zasimov = Zfromp(pExcl(nmean,m,tau,bhat,s)) if not CLExclbool else Zfromp(CLs(nmean,m,tau,bhat,s))
    if asimov_only:
        if more_info:
            print('Mean n = %s\n' %(nmean))
            print('Recommended: (The exact Asimov exclusion significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
        return Zasimov
    else:
        nobs = -1
        probtotal = 0
        foundQuantile = False
        Zmeanlist = []
        pmeanlist = []
        pZlist = []
        if more_info:
            print('Listing all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability:\n')
            print('%-10s %-25s %-25s %-25s %-10s' %('n', 'Probability', 'Cumulative Probability', 'pvalue', 'Significance (ZExcl)'))
        #Considering all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability
        while probtotal < 0.99999999:
            nobs = nobs + 1
            probval = DeltaP(nobs,m,tau,0,bhat)
            pval = pExcl(nobs,m,tau,bhat,s) if not CLExclbool else CLs(nobs,m,tau,bhat,s)
            Zval = Zfromp(pval)
            Zmeanlist.append(probval*Zval)
            pmeanlist.append(probval*pval)
            if Zval > Zcriteria: pZlist.append(probval)
            probtotal = probtotal + probval
            if more_info: print('%-10s %-25s %-25s %-25s %-10s' %(nobs, probval, probtotal, pval, Zval))
            if probtotal > quantile and not foundQuantile:
                #Finding Median Z
                nquantile = nobs
                Zquantile = Zval
                foundQuantile = True
        #Computing Mean Z
        Zmean = sum(Zmeanlist)
        #Computing Mean Z if Z > 0
        Zmean_not_neg = sum([i for i in Zmeanlist if i > 0])
        #Computing Mean pvalue
        pmean = sum(pmeanlist)
        #Converting pmean to significance
        Zpmean = Zfromp(pmean)
        #Computing probability of obtaining a significance  Z > Zcriteria
        pZcriteria = sum(pZlist)
        if more_info:
            print('\nMean n = %s' %(nmean))
            if quantile == 0.5:
                print('Median n = 50%% quantile of the number of pseudo-experiments n = %s\n' %(nquantile))
            else:
                print('%s%% quantile of the number of pseudo-experiments n = %s\n' %(quantile*100, nquantile))
            print('Exclusion significances:\n')
            print('Recommended: (The exact Asimov significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
            print('Reasonable:\n')
            print('\tMean Z = %s\n' %(Zmean))
            print('\tMean Z (Z > 0) = %s\n' %(Zmean_not_neg))
            print('Not Recommended:\n')
            if quantile == 0.5:
                print('\tMedian Z = Z(50%% quantile of the number of pseudo-experiments n) = %s\n' %(Zquantile))
            else:
                print('\tZ(%s%% quantile of the number of pseudo-experiments n) = %s\n' %(quantile*100, Zquantile))
            print('\tZ of Mean p-value = %s\n' %(Zpmean))
            print('Probability of obtaining a significance (ZExcl) greater than Zcriteria = %s is %s\n' %(Zcriteria, pZcriteria))
        return [Zasimov, Zmean, Zmean_not_neg, Zquantile, Zpmean, pZcriteria]

# Upper limit on the signal obtained using the standard frequentist pExcl or
#a Bayesian statistic CLExcl (= modified frequentist CLs for single-channel)
def _sExcl(bhat, dbhat=0, expected_limit=True, CLExclbool=True, n=0, CL=0.9, more_info=False):
    #Trading (bhat, dbhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    if expected_limit:
        #Setting n equal to its expected mean
        if bhat != 0:
            n = bhat + dbhat**2/bhat
        else:
            n = 0
    #Confidence level CL expressed as p-value
    pcriteria = 1 - CL
    #Equation to solve for the upper limit on the signal given bhat, dbhat, CL
    if CLExclbool:
        #CLExcl (Bayesian) (= CLs (modified frequentist) for single-channel)
        solve = lambda s: CLs(n, m, tau, bhat, s) - pcriteria
    else:
        #p-value (standard frequentist)
        solve = lambda s: pExcl(n, m, tau, bhat, s) - pcriteria
    try:
        #Initial guess for mpmath to find root
        if dbhat == 0:
            s_guess = n
        else:
            #To speed up root finding for the uncertain background case, use the result for the known background case as the initial guess
            s_guess = _sExcl(bhat, 0, expected_limit, CLExclbool, n, CL, more_info=False)
        #Solve for root
        result = mp.findroot(solve, s_guess, verbose=more_info)
        #Sometimes the result could be NaN, so raise an error if that is the case
        if mp.isnan(result): raise TypeError
    except:
        #In case of errors while solving for root in the try block, the except block is executed
        try:
            if not expected_limit:
                s_guess = 1 if n <= bhat else n - bhat
            else:
                s_guess = 1
            result = mp.findroot(solve, s_guess, verbose=more_info)
            if mp.isnan(result): raise TypeError
        except:
            if more_info: print("Solution doesn't exist or has numerical instabilities.")
            result=mp.nan
    return float(result)

def sExclObs(n, bhat, dbhat=0, CL=0.9, CLExclbool=True, more_info=False):
    '''
    Inputs
    ------

    n: number of observed events in the signal region.
       This is used only when *expected_limit* is set to *False*

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region.
           Set to *0* by default

    CL: Required confidence level for exclusion. Set to *0.90* by default

    CLExclbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLExcl

    more_info: Set to *False* by default. Set this to *True* to show more
               information about the computation

    Returns
    -------

    The observed upper limit on the signal at *CL* confidence level using
    a more conservative Bayesian statistic CLExcl by default

    If *CLExclbool* is set to *False*, the standard p-value is used instead of
    a more conservative Bayesian statistic CLExcl

    Note: CLExcl is numerically equal to the modified frequentist CLs for
          for a single-channel counting experiment

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    return _sExcl(bhat, dbhat, expected_limit=False, CLExclbool=CLExclbool, n=n, CL=CL, more_info=more_info)

def sExclExp(bhat, dbhat=0, CL=0.9, CLExclbool=True, more_info=False):
    '''
    Inputs
    ------

    bhat: Poisson mean of the background in the signal region

    dbhat: background uncertainty in the signal region.
           Set to *0* by default

    CL: Required confidence level for exclusion. Set to *0.90* by default

    CLExclbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLExcl

    more_info: Set to *False* by default. Set this to *True* to show more
               information about the computation

    Returns
    -------

    The expected upper limit on the signal at *CL* confidence level using
    a more conservative Bayesian statistic CLExcl by default

    If *CLExclbool* is set to *False*, the standard p-value is used instead of
    a more conservative Bayesian statistic CLExcl

    Note: CLExcl is numerically equal to the modified frequentist CLs for
          for a single-channel counting experiment

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    return _sExcl(bhat, dbhat, expected_limit=True, CLExclbool=CLExclbool, n=0, CL=CL, more_info=more_info)

# Exclusion - likelihood ratio method (Asimov)

def ZExclExpKM(s, b, db=0):
    '''
    Inputs
    ------

    s: Poisson mean of signal

    b: Poisson mean of background

    db: uncertainty in background

    Returns
    -------

    The expected significance for exclusion based on a likelihood ratio method

    References
    ----------

    N. Kumar, S. P. Martin [arXiv: 1510.03456 [hep-ph]]
    '''
    if b != 0:
        if db == 0:
            ZKM = np.sqrt(2 * (s - b * np.log(1 + s/b)))
        else:
            x = np.sqrt((s + b)**2 - (4 * s * b * db**2/(b + db**2)))
            ZKM = np.sqrt(2 * (s - b * np.log((b + s + x)/(2 * b)) - (b**2/db**2) * np.log((b - s + x)/(2 * b))) - (b + s - x) * (1 + (b/db**2)))
    else:
        ZKM = np.sqrt(2 * s)
    return ZKM


############### Multi-channel counting experiments ###############

##Discovery - frequentist

def _RestrictQDiscN(karray, narray, barray, sarray, more_info = False):
    '''
    Inputs
    ------

    karray: array of number of events in a mult-channel counting experiment

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    more_info: Set to *False* by default. Set this to *True* to show
               more information about the computation of significance

    Returns
    -------

    Checks if the log-likelihood ratio lnQ(narray) <= lnQ(karray)

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    summ = 0.0
    for i in range(len(karray)):
        if barray[i] == 0.0:
            if narray[i] == karray[i]:
                deltalnQ = 0
            elif narray[i] > karray[i]:
                deltalnQ = np.inf
            elif narray[i] < karray[i]:
                deltalnQ = -np.inf
        else:
            deltalnQ = (narray[i] - karray[i])*np.log(1 + sarray[i]/barray[i])
        summ = summ + deltalnQ
        if more_info:
            print("\tki = %s, ni = %s, bi = %s, s_i = %s:" %(karray[i], narray[i], barray[i], sarray[i]))
            print("\t\tDelta lnQ = %s" %(deltalnQ))
            print("\t\trunning sum = %s\n" %(summ))
    summ = round(summ, 14)
    if more_info: print("total sum = %s\n" %(summ))
    if summ <= 0.0:
        return True
    else:
        return False

def pDiscN_Unc0(narray, barray, sarray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    Returns
    -------

    The frequentist p-value for discovery for multi-channel counting
    experiments with perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    poissonprod = lambda karray: np.prod([DeltaP(karray[i], scipy.inf, scipy.inf, 0, barray[i]) for i in range(len(barray))])

    nmaxarray = [_MaxPoissonCount(sarray[i], barray[i]) for i in range(len(barray))]

    if len(narray) == 1:
        karray_list = [[i] for i in range(nmaxarray[0]+1)]
    elif len(narray) == 2:
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1))
        karray_list = [list(i) for i in karray_list]
    elif len(narray) == 3:
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1), range(nmaxarray[2]+1))
        karray_list = [list(i) for i in karray_list]
    elif len(narray) == 4:
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1), range(nmaxarray[2]+1), range(nmaxarray[3]+1))
        karray_list = [list(i) for i in karray_list]
    else:
        flatten_list = lambda irregular_list:[element for item in irregular_list for element in flatten_list(item)] if type(irregular_list) is list else [irregular_list]
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1), range(nmaxarray[2]+1), range(nmaxarray[3]+1))
        karray_list = [list(i) for i in karray_list]
        for j in range(4, len(narray)):
            karray_list = it.product(karray_list, range(nmaxarray[j]+1))
            karray_list = [list(i) for i in karray_list]
            karray_list = [flatten_list(i) for i in karray_list]

    pval = 0.0

    for karray in karray_list:
        if _RestrictQDiscN(karray, narray, barray, sarray):
            pval = pval + poissonprod(karray)

    return pval

##Discovery - Bayesian

def CLDiscN(narray, marray, tauarray, sarray, barray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    marray: number of observed events in the supposed background-only (off)
            regions in a multi-channel counting experiment

    tauarray: ratio of background means in off and on regions in a
              multi-channel counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    b(hat)array: array of Poisson means of the background in a multi-channel
                 counting experiment

    Returns
    -------

    A Bayesian statistic CLDisc for discovery for multi-channel counting
    experiments

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    poissonprod = lambda mu: np.prod([DeltaP(narray[i], marray[i], tauarray[i], mu*sarray[i], barray[i]) for i in range(len(barray))])

    Num = poissonprod(0.0)
    Norm = poissonprod(1.0)

    BF = Num/Norm

    return BF

def CLDiscN_Unc0(narray, barray, sarray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    Returns
    -------

    A Bayesian statistic CLDisc for discovery for multi-channel counting
    experiments with perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    poissonprod = lambda mu: np.prod([DeltaP(narray[i], scipy.inf, scipy.inf, mu*sarray[i], barray[i]) for i in range(len(barray))])

    Num = poissonprod(0.0)
    Norm = poissonprod(1.0)

    BF = Num/Norm

    return BF


def ZDiscObsN_Unc0(narray, barray, sarray, CLDiscbool=True):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLDisc

    Returns
    -------

    Discovery significance for multi-channel counting experiments with
    perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    if CLDiscbool:
        pres = CLDiscN_Unc0(narray, barray, sarray)
    else:
        pres = pDiscN_Unc0(narray, barray, sarray)

    return Zfromp(pres)

def ZDiscExpN_Unc0(sarray, barray, CLDiscbool=True):
    '''
    Inputs
    ------

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    CLDiscbool: Set to *False* to use the standard frequentist p-value for
                exclusion instead of a Bayesian statistic CLDisc

    Returns
    -------

    Expected (exact Asimov) discovery significance for multi-channel counting
    experiments with perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    narray = [sarray[i] + barray[i] for i in range(len(barray))]

    if CLDiscbool:
        pres = CLDiscN_Unc0(narray, barray, sarray)
    else:
        pres = pDiscN_Unc0(narray, barray, sarray)

    return Zfromp(pres)

def ZDiscObsN(narray, bhatarray, dbhatarray, sarray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    bhatarray: array of Poisson means of the background in a multi-channel
               counting experiment

    dbhatarray: array of background uncertainties in a multi-channel
                counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    Returns
    -------

    Discovery significance for multi-channel counting experiments using a
    Bayesian statistic CLDisc

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    marray = []
    tauarray = []

    for i in range(len(bhatarray)):
        bhattemp = bhatarray[i]
        dbhattemp = dbhatarray[i]
        if dbhattemp != 0.0:
            marray.append(bhattemp**2/dbhattemp**2)
            tauarray.append(bhattemp/dbhattemp**2)
        else:
            marray.append(np.inf)
            tauarray.append(np.inf)

    pres = CLDiscN(narray, marray, tauarray, sarray, bhatarray)

    return Zfromp(pres)

def ZDiscExpN(sarray, bhatarray, dbhatarray):
    '''
    Inputs
    ------

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    bhatarray: array of Poisson means of the background in a multi-channel
               counting experiment

    dbhatarray: array of background uncertainties in a multi-channel
                counting experiment

    Returns
    -------

    Expected (exact Asimov) discovery significance for multi-channel counting
    experiments using a Bayesian statistic CLDisc

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    marray = []
    tauarray = []

    for i in range(len(bhatarray)):
        bhattemp = bhatarray[i]
        dbhattemp = dbhatarray[i]
        if dbhattemp != 0.0:
            marray.append(bhattemp**2/dbhattemp**2)
            tauarray.append(bhattemp/dbhattemp**2)
        else:
            marray.append(np.inf)
            tauarray.append(np.inf)

    narray = [sarray[i] + bhatarray[i] + dbhatarray[i]**2/bhatarray[i] if bhatarray[i] != 0 else sarray[i] for i in range(len(bhatarray))]

    pres = CLDiscN(narray, marray, tauarray, sarray, bhatarray)

    return Zfromp(pres)


##Exclusion - frequentist

def _RestrictQExclN(karray, narray, barray, sarray, more_info = False):
    '''
    Inputs
    ------

    karray: array of number of events in a mult-channel counting experiment

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    more_info: Set to *False* by default. Set this to *True* to show
               more information about the computation of significance

    Returns
    -------

    Checks if the log-likelihood ratio lnQ(karray) >= lnQ(narray)

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    summ = 0.0
    for i in range(len(karray)):
        if barray[i] == 0.0:
            if narray[i] == karray[i]:
                deltalnQ = 0
            elif narray[i] > karray[i]:
                deltalnQ = np.inf
            elif narray[i] < karray[i]:
                deltalnQ = -np.inf
        else:
            deltalnQ = (narray[i] - karray[i])*np.log(1 + sarray[i]/barray[i])
        summ = summ + deltalnQ
        if more_info:
            print("\tki = %s, ni = %s, bi = %s, s_i = %s:" %(karray[i], narray[i], barray[i], sarray[i]))
            print("\t\tDelta lnQ = %s" %(deltalnQ))
            print("\t\trunning sum = %s\n" %(summ))
    summ = round(summ, 14)
    if more_info: print("total sum = %s\n" %(summ))
    if summ >= 0.0:
        return True
    else:
        return False

def pExclN_Unc0(narray, barray, sarray, mu=1):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    mu: Set to *1* by default. Set it to *0* to compute pExcl with signal
        set to 0 in each channel

    Returns
    -------

    The frequentist p-value for exclusion for multi-channel counting
    experiments with perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    poissonprod = lambda karray: np.prod([DeltaP(karray[i], scipy.inf, scipy.inf, mu*sarray[i], barray[i]) for i in range(len(barray))])

    nmaxarray = [_MaxPoissonCount(sarray[i], barray[i]) for i in range(len(barray))]

    if len(narray) == 1:
        karray_list = [[i] for i in range(nmaxarray[0]+1)]
    elif len(narray) == 2:
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1))
        karray_list = [list(i) for i in karray_list]
    elif len(narray) == 3:
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1), range(nmaxarray[2]+1))
        karray_list = [list(i) for i in karray_list]
    elif len(narray) == 4:
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1), range(nmaxarray[2]+1), range(nmaxarray[3]+1))
        karray_list = [list(i) for i in karray_list]
    else:
        flatten_list = lambda irregular_list:[element for item in irregular_list for element in flatten_list(item)] if type(irregular_list) is list else [irregular_list]
        karray_list = it.product(range(nmaxarray[0]+1), range(nmaxarray[1]+1), range(nmaxarray[2]+1), range(nmaxarray[3]+1))
        karray_list = [list(i) for i in karray_list]
        for j in range(4, len(narray)):
            karray_list = it.product(karray_list, range(nmaxarray[j]+1))
            karray_list = [list(i) for i in karray_list]
            karray_list = [flatten_list(i) for i in karray_list]

    pval = 0.0

    for karray in karray_list:
        if _RestrictQExclN(karray, narray, barray, sarray):
            pval = pval + poissonprod(karray)

    return pval

def CLsN_Unc0(narray, barray, sarray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    Returns
    -------

    The modified frequentist CLs for exclusion for multi-channel counting
    experiments with perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    pval = pExclN_Unc0(narray, barray, sarray, mu=1)
    power = pExclN_Unc0(narray, barray, sarray, mu=0)

    return pval/power

##Exclusion - Bayesian

def CLExclN(narray, marray, tauarray, sarray, barray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    marray: number of observed events in the supposed background-only (off)
            regions in a multi-channel counting experiment

    tauarray: ratio of background means in off and on regions in a
              multi-channel counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    b(hat)array: array of Poisson means of the background in a multi-channel
                 counting experiment

    Returns
    -------

    A Bayesian statistic CLExcl for exclusion for multi-channel counting
    experiments

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    stot = sum(sarray)
    poissonprod = lambda sprime: np.prod([DeltaP(narray[i], marray[i], tauarray[i], (sarray[i]/stot)*sprime, barray[i]) for i in range(len(barray))])
    Integrate_poissonprod = lambda s: scipy.integrate.quad(lambda sprime: poissonprod(sprime), s, scipy.inf, limit=200)[0]

    Norm = Integrate_poissonprod(0)
    Num = Integrate_poissonprod(stot)

    return Num/Norm

def CLExclN_Unc0(narray, barray, sarray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    Returns
    -------

    A Bayesian statistic CLExcl for exclusion for multi-channel counting
    experiments with perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    stot = sum(sarray)
    poissonprod = lambda gamma: np.prod([DeltaP(narray[i], scipy.inf, scipy.inf, (sarray[i]/stot)*gamma, barray[i]) for i in range(len(barray))])
    Integrate_poissonprod = lambda gamma: scipy.integrate.quad(lambda gamma: poissonprod(gamma), gamma, scipy.inf, limit=200)[0]

    Norm = Integrate_poissonprod(0)
    Num = Integrate_poissonprod(stot)

    return Num/Norm

def ZExclObsN_Unc0(narray, barray, sarray, CLExclbool=True, CLsbool=False):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    CLExclbool: Set to *False* to use the (modified) frequentist method for
                exclusion instead of a Bayesian statistic CLExcl

    CLsbool: If *CLExclbool* is set to *False*, set this to *True* to use the
             modified frequentist CLs instead of the standard frequentist
             p-value

    Returns
    -------

    Exclusion significance for multi-channel counting experiments with
    perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    if CLExclbool:
        pval = CLExclN_Unc0(narray, barray, sarray)
    else:
        if CLsbool:
            pval = CLsN_Unc0(narray, barray, sarray)
        else:
            pval = pExclN_Unc0(narray, barray, sarray)

    return Zfromp(pval)

def ZExclExpN_Unc0(sarray, barray, CLExclbool=True, CLsbool=False):
    '''
    Inputs
    ------

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    barray: array of Poisson means of the background in a multi-channel
            counting experiment

    CLExclbool: Set to *False* to use the (modified) frequentist method for
                exclusion instead of a Bayesian statistic CLExcl

    CLsbool: If *CLExclbool* is set to *False*, set this to *True* to use the
             modified frequentist CLs instead of the standard frequentist
             p-value

    Returns
    -------

    Expected (exact Asimov) exclusion significance for multi-channel counting
    experiments with perfectly known backgrounds

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    narray=[b for b in barray]

    if CLExclbool:
        pval = CLExclN_Unc0(narray, barray, sarray)
    else:
        if CLsbool:
            pval = CLsN_Unc0(narray, barray, sarray)
        else:
            pval = pExclN_Unc0(narray, barray, sarray)

    return Zfromp(pval)

def ZExclObsN(narray, bhatarray, dbhatarray, sarray):
    '''
    Inputs
    ------

    narray: array of number of observed events in a mult-channel counting
            experiment

    bhatarray: array of Poisson means of the background in a multi-channel
               counting experiment

    dbhatarray: array of background uncertainties in a multi-channel
                counting experiment

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    Returns
    -------

    Exclusion significance for multi-channel counting experiments using a
    Bayesian statistic CLExcl

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    marray = []
    tauarray = []

    for i in range(len(bhatarray)):
        bhattemp = bhatarray[i]
        dbhattemp = dbhatarray[i]
        if dbhattemp != 0.0:
            marray.append(bhattemp**2/dbhattemp**2)
            tauarray.append(bhattemp/dbhattemp**2)
        else:
            marray.append(np.inf)
            tauarray.append(np.inf)

    pres = CLExclN(narray, marray, tauarray, sarray, bhatarray)

    return Zfromp(pres)

def ZExclExpN(sarray, bhatarray, dbhatarray):
    '''
    Inputs
    ------

    sarray: array of Poisson means of the signal in a multi-channel
            counting experiment

    bhatarray: array of Poisson means of the background in a multi-channel
               counting experiment

    dbhatarray: array of background uncertainties in a multi-channel
                counting experiment

    Returns
    -------

    Expected (exact Asimov) exclusion significance for multi-channel counting
    experiments using a Bayesian statistic CLExcl

    References
    ----------

    P.N.Bhattiprolu, S.P.Martin, J.D.Wells [arXiv:2210.07735 [physics.data-an]]
    '''
    marray = []
    tauarray = []

    for i in range(len(bhatarray)):
        bhattemp = bhatarray[i]
        dbhattemp = dbhatarray[i]
        if dbhattemp != 0.0:
            marray.append(bhattemp**2/dbhattemp**2)
            tauarray.append(bhattemp/dbhattemp**2)
        else:
            marray.append(np.inf)
            tauarray.append(np.inf)

    narray = [bhatarray[i] + dbhatarray[i]**2/bhatarray[i] if bhatarray[i] != 0 else 0.0 for i in range(len(bhatarray))]

    pres = CLExclN(narray, marray, tauarray, sarray, bhatarray)

    return Zfromp(pres)
