#ifndef BAYESIANTOOLS_H
#define BAYESIANTOOLS_H

//Standard libraries
#include <iostream>/*Standard i/o library*/
#include <fstream>/*to write output to a file*/
#include <cstdio>/*provides printf*/
#include <string>/*to use string data type*/
#include <cmath>/*provides erf, erfc, lgamma, log, exp, ...*/

//Boost C++ library
#include <boost/math/quadrature/gauss_kronrod.hpp>/*Boost library that provide Gauss-Kronrod integration method (compatible with C++ lambda functions)*/
#include <boost/math/tools/roots.hpp>/*Boost library to find roots of a function*/

//Number of protons per kiloton of water
#define NPROTONS 3.34e32

/*Integration params*/
//True signal selection efficiencies
#define EPS_LOW 0.0
#define EPS_HIGH 1.0

//True background
#define B_LOW 0.0
#define B_HIGH INFINITY

//Proton's partial decay width
#define GAMMA_LOW 0.0
#define GAMMA_HIGH INFINITY

//For maximum depth and tolerance for adaptive Gauss-Kronrod integration
#define MAX_DEPTH 10
#define REL_ERROR 1e-9

/*Root solver settings*/
//Maximum number of iterations while solving for root
#define MAXIT 25

////A scaling factor that is used to bracket the root
//#define ROOT_FACTOR 2

//To print the results of each DeltaP computed
#define PRINT_STEPS 0

//To print the results of integral for each choice of upper limit on proton width
#define PRINT_GAMMA_RESULT 0

using namespace std;
using namespace boost::math::quadrature;
using namespace boost::math::tools;

//Convert significance to a p-value
double pfromZ(double Z){
    return 0.5*erfc(Z/sqrt(2.0));
}

//Likelihood function that takes uncertainties in b and eps into account by using Poissonian background distribution (as in on-off problem)
double DeltaP(double n, double b, double db, double eps, double deps, double exposure, double gam){
    
    //Trading b and db for m (no. of events in off region) and tau (ratio of backgrounds in off and on regions)
    double m = b*b/(db*db);
    double tau = b/(db*db);
    
    //Integrand for triple integral
    auto fullintegrand = [=] (double b1, double eps1) {
        
        double logint;

        logint = - gam*exposure*eps1 - b1 + n*log(gam*exposure*eps1 + b1) - lgamma(n + 1.0);

        double logint_eps=0;

        if ( deps != 0 ) {
            logint_eps = (eps1-eps)/deps;
            logint_eps = -0.5*logint_eps*logint_eps;
            logint_eps = logint_eps + log(sqrt(2.0/M_PI)*(1.0/(deps*(erf((1-eps)/(sqrt(2.0)*deps)) + erf(eps/(sqrt(2.0)*deps))))));
        }
        
        double logint_b=0;

        if ( db != 0 ) {
            //logint_b = - tau*b1 + m*log(tau*b1) - lgamma(m + 1.0);
            logint_b = - tau*b1 + m*log(tau*b1) + log(tau) - lgamma(m + 1.0);
        }

        logint = logint + logint_eps + logint_b;

        return exp(logint);
    };
    
    //Integration: Gauss-Kronrod
    double result=0;
    double error=0; 
    
    auto integrand = [=] (double b1) {
        
        auto integrandeps = [=] (double eps1) {
            return fullintegrand(b1, eps1);
        };

        if (deps != 0) {
            return gauss_kronrod<double, 31>::integrate(integrandeps, EPS_LOW, EPS_HIGH, MAX_DEPTH, REL_ERROR);
        }
        else {
            return integrandeps(eps);
        }
    };
    
    if (db != 0) {
        if ( b <= 5.0 ) { result = gauss_kronrod<double, 41>::integrate(integrand, B_LOW, B_HIGH, MAX_DEPTH, REL_ERROR, &error); }
        else if ( b <= 50.0 ) { result = gauss_kronrod<double, 51>::integrate(integrand, B_LOW, B_HIGH, MAX_DEPTH, REL_ERROR, &error); }
        else if ( b <= 250.0 ) { result = gauss_kronrod<double, 61>::integrate(integrand, B_LOW, B_HIGH, MAX_DEPTH, REL_ERROR, &error); }
        else if ( b <= 3000.0 ) { result = gauss_kronrod<double, 701>::integrate(integrand, B_LOW, B_HIGH, MAX_DEPTH, REL_ERROR, &error); }
        else if ( b > 3000.0 ) { result = gauss_kronrod<double, 1401>::integrate(integrand, B_LOW, B_HIGH, MAX_DEPTH, REL_ERROR, &error); }
    }
    else {
        result = integrand(b);
    }

    if (PRINT_STEPS == 1) {
        printf("\tDeltaP(n = %-3g, b = %-8g, db = %-8g, eps = %-8g, deps = %-8g, exposure = %-6g, gam = %-10g) = %-12g +/- %-12g \n", n, b, db, eps, deps, exposure, gam, result, error);
    }
    
    return result;
}

//Expected number of events
double nExp(double b, double db, double eps, double deps, double exposure, double gam) {
    double nobs;
    
    //To find the signal contribution to mean *n* in the discovery case
    if ( deps != 0) {//Averaging eps with a gaussian distribution
        auto integrandeps = [eps, deps] (double eps1) {
            double logint_eps = 0;
            
            logint_eps = (eps1-eps)/deps;
            logint_eps = -0.5*logint_eps*logint_eps;
            
            return eps1*exp(logint_eps);
        };
        //The factors upfront in the below statement are simply to ensure that the truncated gaussian amounts to unity upon integrating eps1 from 0 to 1
        nobs = gam*exposure*sqrt(2.0/M_PI)*(1.0/(deps*(erf((1-eps)/(sqrt(2.0)*deps)) + erf(eps/(sqrt(2.0)*deps)))))*gauss_kronrod<double, 31>::integrate(integrandeps, EPS_LOW, EPS_HIGH, MAX_DEPTH, REL_ERROR);
    }
    else { nobs = gam*exposure*eps; }

    double temp = b*b + db*db;
    if ( b != 0) { nobs += temp/b; }
    
    return nobs;
}

//Product of likelihoods in each channel
double DeltaPN(const vector<double> &narray, const vector<double> &barray, const vector<double> &dbarray, const vector<double> &epsarray, const vector<double> &depsarray, const vector<double> &exposurearray, double gam) {
    
    int len = barray.size();
    
    double result=1;
    
    //Product of likelihoods obtained using all three methods at HyperK   
    for (int i=0; (i < len) && (result != 0); i++) {
        
        double n = narray[i]; 
        double b = barray[i];
        double db = dbarray[i];
        double eps = epsarray[i];
        double deps = depsarray[i];
        double exposure = exposurearray[i];
        
        result *= DeltaP(n, b, db, eps, deps, exposure, gam); 
            
    }
    
    //if (PRINT_STEPS == 1) {cout << "\n gam = " << gam << ", runtime = " << runtime << ", product = " << result << "\n" << endl;}
    
    return result;
}

double CLExcl(const vector<double> &narray, const vector<double> &barray, const vector<double> &dbarray, const vector<double> &epsarray, const vector<double> &depsarray, const vector<double> &exposurearray, double gam, bool Normalized = true) {

    double error;
    
    auto integrand = [=] (double gamprime) { return DeltaPN(narray, barray, dbarray, epsarray, depsarray, exposurearray, gamprime); };
    
    double result = gauss_kronrod<double, 41>::integrate(integrand, gam, GAMMA_HIGH, MAX_DEPTH, REL_ERROR, &error);

    if ( Normalized ) {
        
        double norm = gauss_kronrod<double, 41>::integrate(integrand, GAMMA_LOW, GAMMA_HIGH, MAX_DEPTH, REL_ERROR, &error);
        //result = result/norm;
        result = ( result != 0 ) ? result/norm : 0.0; 

    }
    
    //if (PRINT_GAMMA_RESULT == 1) {cout << "\nFor upper limit on width = " << gam << ", the integration result is " << result << ", with error " << error << "\n\n" << endl;}
    
    return result;
}

double CLDisc(const vector<double> &narray, const vector<double> &barray, const vector<double> &dbarray, const vector<double> &epsarray, const vector<double> &depsarray, const vector<double> &exposurearray, double gam) {

    int len = barray.size();
    vector<double> zeroarray(len, 0.0);
    
    double result = DeltaPN(narray, barray, dbarray, zeroarray, zeroarray, zeroarray, 0.0)/DeltaPN(narray, barray, dbarray, epsarray, depsarray, exposurearray, gam); 
    
    //if (PRINT_GAMMA_RESULT == 1) {cout << "\nFor upper limit on width = " << gam << ", the integration result is " << result << ", with error " << error << "\n\n" << endl;}
    
    return result;
}

double WidthExclObs(const vector<double> &narray, const vector<double> &barray, const vector<double> &dbarray, const vector<double> &epsarray, const vector<double> &depsarray, const vector<double> &exposurearray, double Np, double CL, double gam_guess, double factor, int niter) {

    double norm = CLExcl(narray, barray, dbarray, epsarray, depsarray, exposurearray, GAMMA_LOW, false);

    auto solve_gam = [=] (double gam) { 

        double CLExclNoNorm = CLExcl(narray, barray, dbarray, epsarray, depsarray, exposurearray, gam, false);
        
        if ( CLExclNoNorm != 0 ) {

            return CLExclNoNorm/norm - 1 + CL;

        }
        else {

            return - 1 + CL;

        }
    };

    ////Options for root solver
    //double factor = ROOT_FACTOR;/*better if root is in between {guess/factor, guess*factor}*/

    const boost::uintmax_t maxit = niter;/*max no. of iterations while finding root*/
    boost::uintmax_t it = maxit;/*actual no. of iterations, updated during root finding*/

    bool is_rising = false;

    int digits = std::numeric_limits<double>::digits;/*max. possible binary digits accuracy for double type*/
    int get_digits = (digits * 1) /2;/*solve for root that is accurate to chosen number of digits*/
    eps_tolerance<double> tol(get_digits);/*tolerance*/

    pair<double, double> r = bracket_and_solve_root(solve_gam, gam_guess, factor, is_rising, tol, it);

    if(it >= maxit){//if root solver is unable to find root with chosen number of max. iterations
        cout << "Unable to locate solution with " << maxit << " iterations, and current best guess is between " << r.first << " and " << r.second << endl;
    }

    double GAMMA = r.first + (r.second - r.first)/2;//root is midway between brackets

    return GAMMA/Np;
}

double WidthExclExp(const vector<double> &barray, const vector<double> &dbarray, const vector<double> &epsarray, const vector<double> &depsarray, const vector<double> &exposurearray, double Np, double CL, double gam_guess, double factor, int niter) {

    int len = barray.size();
    vector<double> narray;
    
    for (int i=0; i < len; i++) {
        
        double b = barray[i];
        double db = dbarray[i];

        narray.push_back(nExp(b, db, 0.0, 0.0, 0.0, 0.0));
            
    }

    double GAMMA = WidthExclObs(narray, barray, dbarray, epsarray, depsarray, exposurearray, Np, CL, gam_guess, factor, niter);

    return GAMMA;
}

double WidthDiscObs(const vector<double> &narray, const vector<double> &barray, const vector<double> &dbarray, const vector<double> &epsarray, const vector<double> &depsarray, const vector<double> &exposurearray, double Np, double Z, double gam_guess, double factor, int niter) {

    double pcriteria = pfromZ(Z);

    auto solve_gam = [=] (double gam) { 
        return CLDisc(narray, barray, dbarray, epsarray, depsarray, exposurearray, gam) - pcriteria;
    };

    ////Options for root solver
    //double factor = ROOT_FACTOR;/*better if root is in between {guess/factor, guess*factor}*/

    const boost::uintmax_t maxit = niter;/*max no. of iterations while finding root*/
    boost::uintmax_t it = maxit;/*actual no. of iterations, updated during root finding*/

    bool is_rising = false;

    int digits = std::numeric_limits<double>::digits;/*max. possible binary digits accuracy for double type*/
    int get_digits = (digits * 1) /2;/*solve for root that is accurate to chosen number of digits*/
    eps_tolerance<double> tol(get_digits);/*tolerance*/

    pair<double, double> r = bracket_and_solve_root(solve_gam, gam_guess, factor, is_rising, tol, it);

    if(it >= maxit){//if root solver is unable to find root with chosen number of max. iterations
        cout << "Unable to locate solution with " << maxit << " iterations, and current best guess is between " << r.first << " and " << r.second << endl;
    }

    double GAMMA = r.first + (r.second - r.first)/2;//root is midway between brackets

    return GAMMA/Np;
}

double WidthDiscExp(const vector<double> &barray, const vector<double> &dbarray, const vector<double> &epsarray, const vector<double> &depsarray, const vector<double> &exposurearray, double Np, double Z, double gam_guess, double factor, int niter) {

    int len = barray.size();
    double pcriteria = pfromZ(Z);

    auto solve_gam = [=] (double gam) { 
        vector<double> narray;
        
        for (int i=0; i < len; i++) {
            
            double b = barray[i];
            double db = dbarray[i];
            double eps = epsarray[i];
            double deps = depsarray[i];
            double exposure = exposurearray[i];

            narray.push_back(nExp(b, db, eps, deps, exposure, gam));
                
        }

        return CLDisc(narray, barray, dbarray, epsarray, depsarray, exposurearray, gam) - pcriteria;
    };

    ////Options for root solver
    //double factor = ROOT_FACTOR;/*better if root is in between {guess/factor, guess*factor}*/

    const boost::uintmax_t maxit = niter;/*max no. of iterations while finding root*/
    boost::uintmax_t it = maxit;/*actual no. of iterations, updated during root finding*/

    bool is_rising = false;

    int digits = std::numeric_limits<double>::digits;/*max. possible binary digits accuracy for double type*/
    int get_digits = (digits * 1) /2;/*solve for root that is accurate to chosen number of digits*/
    eps_tolerance<double> tol(get_digits);/*tolerance*/

    pair<double, double> r = bracket_and_solve_root(solve_gam, gam_guess, factor, is_rising, tol, it);

    if(it >= maxit){//if root solver is unable to find root with chosen number of max. iterations
        cout << "Unable to locate solution with " << maxit << " iterations, and current best guess is between " << r.first << " and " << r.second << endl;
        double GAMMA = r.first + (r.second - r.first)/2;//root is midway between brackets
    }

    double GAMMA = r.first + (r.second - r.first)/2;//root is midway between brackets

    return GAMMA/Np;
}

#endif
