#Python tools to compute the "exact Asimov significance" advocated in arXiv: 2009.07249 [physics.data-an] as the standard significance measure for projected exclusions and discovery sensitivities in counting experiments.

#Main Functions: (More information is given as a docstring to the functions, which can be accessed using the *help* command)
#Zdisc(s,bhat,dbhat): computes the discovery significance given the Poisson means of signal *s* and background *bhat* events, and Poisson uncertainty *dbhat* of the background. *dbhat* is set to 0 by default.
#Zexcl(s,bhat,dbhat): computes the exclusion significance given the Poisson means of signal *s* and background *bhat* events, and Poisson uncertainty *dbhat* of the background. *dbhat* is set to 0 by default.

#Import as a python package:
#If Zstats package is in path
#>>> import Zstats as Z
#>>> Z.Zdisc(5, 1) #for example
#>>> help(Z.Zdisc) #or
#>>> print(Z.disc.__doc__) #for more info about the usage of *Zdisc* function

#Python packages used
import numpy as np
import scipy
import scipy.special as sc
from scipy.integrate import quad
import mpmath as mp #mpmath can be made faster by installing *gmpy2*

#python/package versions on which this package was built:
#python: '3.8.1' (also tested on '3.7.0')
#numpy: '1.19.1' [latest version (at this writing)]
#scipy: '1.4.1' [also tested on the latest version '1.5.2']
#mpmath: '1.1.0' [latest version]

#Conversions from p to Z and Z to p
def Zfromp(p):
    '''
    Converts one-sided p-value to significance Z using:
    Z = Sqrt[2] InverseErfc[2 p]
    '''
    return np.sqrt(2)*sc.erfcinv(2*p)

def pfromZ(Z):
    '''
    Converts significance Z to one-sided p-value using:
    p = 0.5 Erfc[Z/Sqrt[2]]
    '''
    return 0.5*sc.erfc(Z/np.sqrt(2))

def DeltaP(n,m,tau,s,b):
    '''
    Parameters
    ----------

    n: number of Poisson events observed in the signal (on) region
    m: number of Poisson events observed in the supposed background-only (off) region
    tau: ratio of background means in off and on regions
    s: Poisson mean of the signal
    b: Poisson mean of the background in the on region

    Note: In the case of uncertain background,
    *m* and *tau* should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    dbhat: Poisson uncertainty of the background in the on region

    Returns
    -------

    Known background case (infinite *m* and *tau*): returns the Poisson probability of observing n-events given mean = s+b (signal+background events).

    Uncertain background case (finite *m* and *tau*): returns the probability of observing *n* events in the on region, given *m* events in the off-region
    and a number *tau* (ratio of Poisson means in off- and on-regions) in the on-off problem, and number of signal events *s*.
    (here, *b* is a dummy variable and should not affect the result)

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P. N. Bhattiprolu, S. P. Martin, J. D. Wells [arXiv: 2009.07249 [physics.data-an]]
    '''
    n = float(n)
    m = float(m)
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

def pvalue(n, m, tau, s, b, more_info=False):
    '''
    Parameters
    ----------

    n: number of Poisson events observed in the signal (on) region
    m: number of Poisson events observed in the supposed background-only (off) region
    tau: ratio of background means in off and on regions
    s: Poisson mean of the signal
    b(hat): Poisson mean of the background (in the on region)
    more_info: Set to False (or 0) by default. Set this to True to show more information
    about the computation of p-values.

    Note: In the case of uncertain background,
    m and tau should be thought of as proxy for
    bhat: Poisson mean of the background in the on region
    dbhat: Poisson uncertainty of the background in the on region

    Returns
    -------

    p-value for:

    Exclusion case - when s = 0
    Discovery case - when s > 0
               &
    Known background case - when *m* and *tau* are infinite
    Uncertain background case - when *m* and *tau* are finite
    (here *b* is a dummy variable and should not affect the result)

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P. N. Bhattiprolu, S. P. Martin, J. D. Wells [arXiv: 2009.07249 [physics.data-an]]
    '''
    n = float(n)
    m = float(m)
    if more_info:
        print('Inputs:')
        print('n = %s' %(n))
        print('m = %s' %(m))
        print('tau = %s' %(tau))
        print('s = %s' %(s))
        print('b = %s\n' %(b))
    if m == scipy.inf and tau == scipy.inf:
        if s == 0:
            #p-value for discovery when background is known exactly.
            pBi = sc.gammainc(n,b) #has round-off errors for large n and small b, for example, when (n,b)=(>=118,0.1)
            if more_info:
                print('p-value for discovery case when background is known exaclty:\n')
                print('\tp_disc = %s\n' %(pBi))
        else:
            #p-value for exclusion when background is known exactly.
            pBi = sc.gammaincc(n+1,s+b)
            if more_info:
                print('p-value for exclusion case when background is known exaclty:\n')
                print('\tp_excl = %s\n' %(pBi))
    else:
        if s == 0:
            #p-value for discovery when background is not known exactly.
            w = 1/(1+tau)
            pBi = sc.betainc(n,m+1,w) if n != 0 else 1.0
            if more_info:
                print('p-value for discovery case when background is not known exaclty:\n')
                print('\tp_disc = %s\n' %(pBi))
        else:
            if more_info: print('p-value for exclusion case when background is not known exaclty:\n')
            #p-value for exclusion when background is not known exactly.
            if n.is_integer():
                pBi = float((mp.power(tau,m+1)/mp.gamma(m+1))*sum([(mp.gamma(k+m+1)/(mp.factorial(k)*mp.power(1+tau,k+m+1)))*sc.gammaincc(n-float(k)+1,s) for k in mp.arange(n+1)]))
                if more_info:
                    print('Here, n is an integer. p_excl evaluated using the sum formula is:\n')
                    print('\tp_excl = %s\n' %(pBi))
            else:
                #only used for Asimov approximations in the exclusion case (if n is not an integer)
                pBilist = []

                #first computing pBisum = p_excl(l,m,tau,s) using sum formula with l = closest integer to n
                pBisum = float((mp.power(tau,m+1)/mp.gamma(m+1))*sum([(mp.gamma(k+m+1)/(mp.factorial(k)*mp.power(1+tau,k+m+1)))*sc.gammaincc(n-float(k)+1,s) for k in mp.arange(n+1)]))#always reliable
                pBilist.append(pBisum)


                #then computing pBiint = p_excl(n,m,tau,s) using the first integral formula
                #Computing with scipy's quad and gammaincc functions
                pBiint = float((mp.power(tau,m+1)/mp.gamma(m+1))*quad(lambda bb: mp.power(bb,m)*mp.exp(-bb*tau)*sc.gammaincc(n+1,s+bb),0,scipy.inf)[0])#in principle same as pBibyparts, but may have numerical instabilities

                #To compute the integral to arbitrary precision, use mpmath's quad (for numerical integration) and inc. gamma functions
                #To do so, uncomment the below two lines and comment out the line somewhere above, starting with *pBiint*, which uses scipy to do the integral*
                #mp.mp.dps = 15 #This sets the global *decimal precision* of *mpmath*. In case of errors, try increasing this number (which in turn increases the computation time).
                #pBiint = float((mp.power(tau,m+1)/mp.gamma(m+1))*mp.quad(lambda bb: mp.power(bb,m)*mp.exp(-bb*tau)*mp.gammainc(n+1,s+bb,regularized=True),[0,mp.inf]))

                if pBiint >= 0 and not (np.isnan(pBiint) or np.isinf(pBiint)):
                    percentdiffint = 100*2*abs(pBisum-pBiint)/(pBisum+pBiint)
                    if percentdiffint < 20:
                        pBilist.append(pBiint)

                #finally computing pBibyparts = p_excl(n,m,tau,s) using the second integral formula
                #Computing with scipy's quad and gammaincc functions
                pBibyparts = float(sc.gammaincc(n+1,s) - (quad(lambda bb: mp.power(s+bb,n)*mp.exp(-(s+bb))*sc.gammaincc(m+1,tau*bb),0,scipy.inf)[0]/mp.gamma(n+1)))#can be inaccurate when n >> b

                #To compute the integral to arbitrary precision, use mpmath's quad (for numerical integration) and inc. gamma functions
                #To do so, uncomment the below two lines and comment out the line somewhere above, starting with *pBibyparts*, which uses scipy to do the integral*
                #mp.mp.dps = 15 #This sets the global *decimal precision* of *mpmath*. In case of errors, try increasing this number (which in turn increases the computation time).
                #pBibyparts = float(mp.gammainc(n+1,s,regularized=True) - (mp.quad(lambda bb: mp.power(s+bb,n)*mp.exp(-(s+bb))*mp.gammainc(m+1,tau*bb,regularized=True),[0,mp.inf])/mp.gamma(n+1)))

                if pBibyparts >= 0 and not (np.isnan(pBibyparts) or np.isinf(pBibyparts)):
                    percentdiffbyparts = 100*2*abs(pBisum-pBibyparts)/(pBisum+pBibyparts)
                    if percentdiffbyparts < 20:
                        pBilist.append(pBibyparts)

                #Find the best solution among {pBisum (always reliable), pBibyparts, pBiint}
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
                    print('First, computing psum = p_excl(l,m,tau,s) using a sum formula with l = closest integer to n: psum = %s' %(pBisum))
                    print('(Note: psum is an approximate answer when n != integer, but can be computed more reliably)\n')
                    print('Then, computing pint1 = p_excl(n,m,tau,s) using the first integral formula: pint1 = %s\n' %(pBiint))
                    print('Finally, computing pint2 = p_excl(n,m,tau,s) using the second integral formula: pint2 = %s\n' %(pBibyparts))
                    print('Now, finding the closest of {pint1, pint2} to psum:')
                    print('\tp_excl = %s\n' %(pBi))
                    print('Note: In some cases only one of the integral formulas gives a sensible result. If so, this funtion returns that.')
                    print('And, in some other rare cases, both pint1 and pint2 are not sensible, where, this function returns psum instead, which is an approximation to the true result.\n')
    return pBi

def Zdisc(s, bhat, dbhat=0, asimov_only=True, Zcriteria=5.0, more_info=False):
    '''
    Parameters
    ----------

    s: Poisson mean of the signal

    bhat: Poisson mean of the background in the signal region

    dbhat: Poisson uncertainty of the background in the signal region.
    Set to *0* by default

    asimov_only: By default this is set to *True* (or *1*) and the function *Zdisc*
    only returns the exact Asimov discovery significance. Instead if this option
    is set to *False* (or *0*), the function *Zdisc* returns a list of
    {exact Asimov significance,
    Mean Z with Z != -Infinity,
    Mean Z with Z > 0,
    Median Z,
    Z obtained from Mean p-value,
    Probability of obtaining Zdisc > Zcriteria}

    Zcriteria: Set to 5.0 sigma by default. This is used only when *asimov_only* is
    set to *False* to compute the probability of obtaining Zdisc > Zcriteria

    more_info: Set to *False* (or *0*) by default. Set this to *True* to show more information
    about the computation of significance(s).

    Returns
    -------

    This function returns the expected significance(s) for discovery, given *s*,
    *b* and *dbhat*

    If the option *asimov_only* is set to *True* (default case): returns the exact
    Asimov discovery significance

    If the option *asimov_only* is set to *False*: returns a list of
    {exact Asimov significance,
    Mean Z with Z != -Infinity,
    Mean Z with Z > 0,
    Median Z,
    Z obtained from Mean p-value,
    Probability of obtaining Zdisc > Zcriteria}

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P. N. Bhattiprolu, S. P. Martin, J. D. Wells [arXiv: 2009.07249 [physics.data-an]]
    '''
    #Trading (bhat, Deltabhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    #Mean n
    nmean = s + bhat + dbhat**2/bhat
    if more_info:
        print('Inputs:')
        print('s = %s' %(s))
        print('b hat = %s' %(bhat))
        print('Deltabhat = %s\n' %(dbhat))
        print('Trading (bhat, Deltabhat) for (m, tau):')
        print('m = %s' %(m))
        print('tau = %s\n' %(tau))
    if asimov_only:
        #Computing exact Asimov significance
        Zasimov = Zfromp(pvalue(nmean,m,tau,0,bhat))
        if more_info:
            print('Mean n = %s\n' %(nmean))
            print('Recommended: (The exact Asimov discovery significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
        return Zasimov
    else:
        #Computing exact Asimov significance
        Zasimov = Zfromp(pvalue(nmean,m,tau,0,bhat))
        nobs = -1
        probtotal = 0
        foundMedian = False
        Zmeanlist = []
        pmeanlist = []
        pZlist = []
        if more_info:
            print('Listing all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability:\n')
            print('%-10s %-25s %-25s %-25s %-10s' %('n', 'Probability', 'Cumulative Probability', 'pvalue', 'Significance (Zdisc)'))
        #Considering all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability
        while probtotal < 0.99999999:
            nobs = nobs + 1
            probval = DeltaP(nobs,m,tau,s,bhat)
            pval = pvalue(nobs,m,tau,0,bhat)
            Zval = Zfromp(pval)
            Zmeanlist.append(probval*Zval)
            pmeanlist.append(probval*pval)
            if Zval > Zcriteria: pZlist.append(probval)
            probtotal = probtotal + probval
            if more_info: print('%-10s %-25s %-25s %-25s %-10s' %(nobs, probval, probtotal, pval, Zval))
            if probtotal > 0.5 and not foundMedian:
                #Finding Median Z
                nmedian = nobs
                Zmedian = Zval
                foundMedian = True
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
            print('Median n = %s\n' %(nmedian))
            print('Discovery significances:\n')
            print('Recommended: (The exact Asimov significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
            print('Reasonable:\n')
            print('\tMean Z (Z != -Infinity) = %s\n' %(Zmean_not_neginf))
            print('\tMean Z (Z > 0) = %s\n' %(Zmean_not_neg))
            print('Not Recommended:\n')
            print('\tMedian Z = %s\n' %(Zmedian))
            print('\tZ of Mean p-value = %s\n' %(Zpmean))
            print('Probability of obtaining a significance (Zdisc) greater than Zcriteria = %s is %s\n' %(Zcriteria, pZcriteria))
        return [Zasimov, Zmean_not_neginf, Zmean_not_neg, Zmedian, Zpmean, pZcriteria]

def Zexcl(s, bhat, dbhat=0, asimov_only=True, Zcriteria = 1.645, more_info=False):
    '''
    Parameters
    ----------

    s: Poisson mean of the signal

    bhat: Poisson mean of the background in the signal region

    dbhat: Poisson uncertainty of the background in the signal region.
    Set to *0* by default

    asimov_only: By default this is set to *True* (or *1*) and the function *Zexcl*
    only returns the exact Asimov exclusion significance. Instead if this option
    is set to *False* (or *0*), the function *Zexcl* returns a list of
    {exact Asimov significance,
    Mean Z,
    Mean Z with Z > 0,
    Median Z,
    Z obtained from Mean p-value,
    Probability of obtaining Zexcl > Zcriteria}

    Zcriteria: Set to 1.645 sigma by default. This is used only when *asimov_only* is
    set to *False* to compute the probability of obtaining Zexcl > Zcriteria

    more_info: Set to *False* (or *0*) by default. Set this to *True* to show more information
    about the computation of significance(s).

    Returns
    -------

    This function returns the expected significance(s) for exclusion, given *s*,
    *b* and *dbhat*

    If the option *asimov_only* is set to *True* (default case): returns the exact
    Asimov exclusion significance

    If the option *asimov_only* is set to *False*: returns a list of
    {exact Asimov significance,
    Mean Z,
    Mean Z with Z > 0,
    Median Z,
    Z obtained from Mean p-value,
    Probability of obtaining Zexcl > Zcriteria}

    Relations
    ---------

    bhat = m/tau
    dbhat = Sqrt[m]/tau

    References
    ----------

    P. N. Bhattiprolu, S. P. Martin, J. D. Wells [arXiv: 2009.07249 [physics.data-an]]
    '''
    #Trading (bhat, Deltabhat) for (m, tau)
    if dbhat != 0:
        m = bhat**2/dbhat**2
        tau = bhat/dbhat**2
    else:
        m = scipy.inf
        tau = scipy.inf
    #Mean n
    nmean = bhat + dbhat**2/bhat
    if more_info:
        print('Inputs:')
        print('s = %s' %(s))
        print('bhat = %s' %(bhat))
        print('Deltabhat = %s\n' %(dbhat))
        print('Trading (bhat, Deltabhat) for (m, tau):')
        print('m = %s' %(m))
        print('tau = %s\n' %(tau))
    if asimov_only:
        #Computing exact Asimov significance
        Zasimov = Zfromp(pvalue(nmean,m,tau,s,bhat))
        if more_info:
            print('Mean n = %s\n' %(nmean))
            print('Recommended: (The exact Asimov exclusion significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
        return Zasimov
    else:
        #Computing exact Asimov significance
        Zasimov = Zfromp(pvalue(nmean,m,tau,s,bhat))
        nobs = -1
        probtotal = 0
        foundMedian = False
        Zmeanlist = []
        pmeanlist = []
        pZlist = []
        if more_info:
            print('Listing all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability:\n')
            print('%-10s %-25s %-25s %-25s %-10s' %('n', 'Probability', 'Cumulative Probability', 'pvalue', 'Significance (Zexcl)'))
        #Considering all pseudo-experiment outcomes in the signal-region that can possibly have a significant probability
        while probtotal < 0.99999999:
            nobs = nobs + 1
            probval = DeltaP(nobs,m,tau,0,bhat)
            pval = pvalue(nobs,m,tau,s,bhat)
            Zval = Zfromp(pval)
            Zmeanlist.append(probval*Zval)
            pmeanlist.append(probval*pval)
            if Zval > Zcriteria: pZlist.append(probval)
            probtotal = probtotal + probval
            if more_info: print('%-10s %-25s %-25s %-25s %-10s' %(nobs, probval, probtotal, pval, Zval))
            if probtotal > 0.5 and not foundMedian:
                #Finding Median Z
                nmedian = nobs
                Zmedian = Zval
                foundMedian = True
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
            print('Median n = %s\n' %(nmedian))
            print('Exclusion significances:\n')
            print('Recommended: (The exact Asimov significance)\n')
            print('\tAsimov Z = %s\n' %(Zasimov))
            print('Reasonable:\n')
            print('\tMean Z = %s\n' %(Zmean))
            print('\tMean Z (Z > 0) = %s\n' %(Zmean_not_neg))
            print('Not Recommended:\n')
            print('\tMedian Z = %s\n' %(Zmedian))
            print('\tZ of Mean p-value = %s\n' %(Zpmean))
            print('Probability of obtaining a significance (Zexcl) greater than Zcriteria = %s is %s\n' %(Zcriteria, pZcriteria))
        return [Zasimov, Zmean, Zmean_not_neg, Zmedian, Zpmean, pZcriteria]

def ZdiscAsimovCCGV(s, b, db=0):
    '''
    Parameters
    ----------

    s: Poisson mean of signal

    b: Poisson mean of background

    db: Poisson uncertainty of background

    Returns
    -------

    Computes the expected significance for discovery based on a likelihood ratio method

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

def ZexclAsimovKM(s, b, db=0):
    '''
    Parameters
    ----------

    s: Poisson mean of signal

    b: Poisson mean of background

    db: Poisson uncertainty of background

    Returns
    -------

    Computes the expected significance for exclusion based on a likelihood ratio method

    References
    ----------

    N. Kumar, S. P. Martin [arXiv: 1510.03456 [hep-ph]]
    '''
    if db == 0:
        ZKM = np.sqrt(2 * (s - b * np.log(1 + s/b)))
    else:
        x = np.sqrt((s + b)**2 - (4 * s * b * db**2/(b + db**2)))
        ZKM = np.sqrt(2 * (s - b * np.log((b + s + x)/(2 * b)) - (b**2/db**2) * np.log((b - s + x)/(2 * b))) - (b + s - x) * (1 + (b/db**2)))
    return ZKM
