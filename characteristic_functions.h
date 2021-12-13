#ifndef characteristic_functions_h
#define characteristic_functions_h
#include "global_variables.h"
#include <iostream>
#include <complex>
#include <string>
using namespace std;

complex<double>* chf_black_scholes(complex<double> *x_arr, int n, double S0, double tau, double r, double d, double sigma){
    cout << "it's the characteristic function of Black Scholes" << endl;
    complex<double> *result = new complex<double>[n];
    for(int idx=0; idx<n; idx++){
        // (log(S0) + r*tau - d*tau - 0.5*sigma*sigma*tau)*(1i)*x_arr
        complex<double> former_term = (log(S0) + r*tau - d*tau - 0.5*sigma*sigma*tau)*(1i)*x_arr[idx];
        // (-0.5*sigma*sigma*tau)*x_arr*x_arr
        complex<double> latter_term = (-0.5*sigma*sigma*tau)*x_arr[idx]*x_arr[idx];
        result[idx] = exp(former_term + latter_term);
    }
    return result;
}

void chf_heston(){
    cout << "it's the characteristic function of Heston" << endl;
}

void chf_bates(){
    cout << "it's the characteristic function of Bates" << endl;
}

#endif /* characteristic_functions_h */
