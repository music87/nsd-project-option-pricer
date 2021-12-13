#include <iostream>
#include <vector>
#include <numbers>
#include <x86intrin.h>
#include "global_variables.h"
//#include "basic_arithmetic_operations.h"
#include "characteristic_functions.h"
//#include "fourier_transform_modified_call_price.h"
#include "fast_fourier_transfrom_call_pricer.h"
using namespace std;

int main(int argc, const char * argv[]) {
    cout << "start ... " << endl;
    // hyper-parameter
    int N = pow(2, 20);
    double eta = 0.01;
    double alpha = 1.0;
    
    // option contract
    double r = 0.06;
    double tau = 0.4;
    double d = 0.01;
    double n_contract = 7;
    double sigma = 0.2;
    double* strike_prices = new double[n_contract];
    strike_prices[0]=5.0;
    strike_prices[1]=10.0;
    strike_prices[2]=30.0;
    strike_prices[3]=50.0;
    strike_prices[4]=55.0;
    strike_prices[5]=60.0;
    strike_prices[6]=100.0;
    
    // market information
    double S0 = 50.0;
    
    // test call option pricer
    double* call_prices = fft_call_pricer(N, eta, alpha, S0, tau, r, d, strike_prices, n_contract, &sigma);
    
    cout << call_prices[3] << endl;
}
