#ifndef FCTOOLS_H
#define FCTOOLS_H

//C++ header file to compute the Feldman-Cousins (FC) upper limits, aacceptance regions, and experimental sensitivity
//for Poisson distributions with perfectly known background.
//(Based on the procedure listed in https://arxiv.org/pdf/physics/9711021.pdf)

//Some example programs (*.cc files) that show the usage of the functions in this header are in the *examples* folder

#include <iostream>/*standard i/o library*/
#include <fstream>/*to write output to a file*/
#include <string>/*to use string data type*/
#include <cstdio>/*provides printf*/
#include <cmath>/*provides erf, erfc, lgamma, log, exp, ...*/
#include <algorithm>/*provides max, sort, minmax_element*/
#include <vector>/*provides vector*/

#define S_STEP 0.0025/*stepsize by which signal is increased while computing upper limits*/
#define S_MAX 500.0/*max signal while computing (fake) upper limit*/
#define S_MAXSTEP 1.0/*max signal while computing true upper limit*/
#define PROB_MAX 0.99999999/*max cumulative probability considered*/

using namespace std;

//Poisson probability to observe *n* events given a mean *mu*
double Poisson(double n, double mu) {
    double logresult = -mu;
    if ( n != 0 ) { logresult += n*log(mu) - lgamma(n + 1.0); }
    return exp(logresult);
}

//Function to find the values of *n* that are in the acceptance region of a chosen *CL* confidence intervals, given signal *s* and background *b*
pair<int, int> AcceptanceInterval(double s, double b, double CL, bool PRINT_DETAILS=0) {

    if ( PRINT_DETAILS ) {
        printf("\n(s, b, CL) = (%g, %g, %g)\n\n", s, b, CL);
    }

    int nobs = -1;
    double probtotal = 0;
    double probval, sbest, probvalbest, R;
    
    struct FCparams
    {
        int NOBS;/*Observed number of events*/
        double PROBVAL;/*Poisson(n | s + b)*/
        double RVAL;/*R is the likelihood ratio defined in the FC paper*/
    };

    vector<FCparams> array;

    if ( PRINT_DETAILS ) {
        //To print data similar to that of Table 1 of the FC paper
        printf("Listing all pseudo-experiment outcomes that can possibly have a significant probability:\n");
        printf("%-10s %-25s %-10s %-25s %-25s %-25s\n", "n", "Prob", "Best s", "Prob with best s", "R", "Cumulative Probability");
    }

    //Considering all pseudo-experiment outcomes that can possibly have a significant probability
    while ( probtotal < PROB_MAX ) {
        nobs = nobs + 1;
        probval = Poisson(nobs, s + b);
        sbest = max(0.0, nobs-b);
        probvalbest = Poisson(nobs, sbest + b);
        probtotal = probtotal + probval;
        R = probval/probvalbest;//probabilities in probval with largest R are assigned the lowest rank
        if ( PRINT_DETAILS ) { 
            printf("%-10d %-25g %-10g %-25g %-25g %-25.8f\n", nobs, probval, sbest, probvalbest, R, probtotal);
        }
        array.push_back({nobs, probval, R});//collecting the probability *probval* and likelihood ratio *R* for each outcome *nobs*
    }

    //Sort the array in descending order of R 
    sort(array.begin(), array.end(), [] (const FCparams &i, const FCparams &j) { return i.RVAL > j.RVAL; } );
    
    //Now, making a list of all n with decreasing R until the sum of the corresponding probabilities adds up to or exceeds the chosen CL
    int i=-1;
    double CLcheck=0;
    vector<int> narray;
    
    while ( CLcheck < CL ) {
        i = i + 1;
        CLcheck =  CLcheck + array[i].PROBVAL;
        narray.push_back(array[i].NOBS);
    }
    
    if ( PRINT_DETAILS ) {
        printf("\n\n Acceptance Interval for n: %d to %d\n\n", *min_element(begin(narray), end(narray)), *max_element(begin(narray), end(narray)));
    }
    
    return {*min_element(begin(narray), end(narray)), *max_element(begin(narray), end(narray))};
    }

//Feldman-Cousins (FC) upper limits on the signal when background *b* and count *n* are known
double _FakeUpperLimit(int n, double b, double CL=0.9, double s_min=0.0, double s_max=S_MAX, bool PRINT_DETAILS=0) {
    //Finding conservative upper and lower limits on the signal, given background and observed integer n, by
    //matching the minimum and maximum of the acceptance interval [n_min, n_max] obtained using the above
    //procedure, by varying s in small steps, with the observed n and (n-1) respectively.
    
    //Also, this function does not necessarily give the true upper limit.
    //(See https://docs.gammapy.org/0.9/stats/feldman_cousins.html)
    //The true upper limit is computed below with the function UpperLimit()

    if (PRINT_DETAILS) {
        printf("\n(n, b, CL) = (%d, %g, %g)\n", n, b, CL);
    }
    
    vector<double> sUL_array;
    pair<int, int> nint;

    bool FoundSolution = 0;
    for ( double s = max(s_min, n-b); s < s_max; s += S_STEP ) {
        //Compute the acceptance region for n for signal mean *s* 
        nint = AcceptanceInterval(s, b, CL);
        
        if (PRINT_DETAILS) {
            printf("\nFor s = %g, the acceptance interval for n: %d to %d\n", s, nint.first, nint.second);
        }
        
        //Check if if the input *n* is in the acceptance region, repeat for various *s* till *n* is not in the acceptance region
        if ( nint.first == n ) {
            sUL_array.push_back(s);
            FoundSolution=1;
        }
        
        if ( (FoundSolution) && ( nint.first > n ) ) { break; }
    }
    
    if ( sUL_array.size() > 0 ) {

        if (PRINT_DETAILS) {
            printf("(Fake) Upper limit on s: %g\n\n", *max_element(begin(sUL_array), end(sUL_array)));
        }

        return *max_element(begin(sUL_array), end(sUL_array));
    }
    else {
        return -1.0;
    }
}

//Actual upper limit 
double UpperLimit(int n, double b, double CL=0.9, bool PRINT_DETAILS=0) {

    if ( PRINT_DETAILS ) {
        printf("\n(n, b, CL) = (%d, %g, %g)\n\n", n, b, CL);
    }

    double temp = _FakeUpperLimit(n, b, CL, 0.0, S_MAX, PRINT_DETAILS);
    double result = temp;
    
    if (PRINT_DETAILS) {
        printf("\n\n(Fake) upper limit on s: %g\n\n", result); 
    }

    while ( temp != -1.0 ) {
        result = temp;
        temp = _FakeUpperLimit(n, b, CL, result+S_STEP, result+S_MAXSTEP);

        if ( (PRINT_DETAILS) and (temp != 1.0) ) {
            printf("\nUpper limit on s: %g\n", result); 
        }
    }
    
    if (PRINT_DETAILS) {
        printf("\n\nActual upper limit on s: %g\n\n", result); 
    }
    return result;
}

//Function to calculate the average upper limit on the signal using the Feldman-Cousins method 
double Sensitivity(double b, double CL=0.9, bool PRINT_DETAILS=0) {

    if ( PRINT_DETAILS ) {
        printf("\n(b, CL) = (%g, %g)\n\n", b, CL);
    }

    int nobs = -1;
    double probval;
    double upperlimit;
    double probtotal = 0;
    double avg_FC = 0;

    if ( PRINT_DETAILS ) {
        printf("Listing all pseudo-experiment outcomes that can possibly have a significant probability:\n");
        printf("%-10s %-25s %-25s %-25s %-25s\n", "n", "Prob", "Upper limit", "Sensitivity (so far)", "Cumulative Probability");
    }

    while( probtotal < PROB_MAX ) {
        nobs = nobs + 1;
        probval = Poisson(nobs, b);
        upperlimit = UpperLimit(nobs, b, CL);
        avg_FC += probval*upperlimit;
        probtotal += probval;

        if ( PRINT_DETAILS ) { 
            printf("%-10d %-25g %-25g %-25g %-25.8f\n", nobs, probval, upperlimit, avg_FC, probtotal);
        }
    }
    
    if ( PRINT_DETAILS ) { 
        printf("\n\nExperimental sensitivity = %g\n\n", avg_FC);
    }

    return avg_FC;
    }

#endif
