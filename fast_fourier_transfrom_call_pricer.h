#pragma once
#include "global_variables.h"
#include <iostream>
#include <complex>
#include <string>
#include <stdexcept>
#include <vector>
#include <array>
using namespace std;

int bit_reversal(int x, int bit_length) {
    int y=0;
    // for each bit of x
    for (int bit_loc=0; bit_loc < bit_length; bit_loc++) {
        y <<= 1; // shift 1 bit from right to left eg. 10->00, 01->10, 1011->0110
        y = y|(x&1); // get the last bit; note that 1 means 000...001
        x >>= 1; // shift 1 bit from left to right eg. 1011->0101, 110101->011010
    }
    return y;
}

template<class T>
complex<double>* fft(T* const x, int bit_length){
    int N = 1<<bit_length; // N = 2^{bit_length}
    complex<double> *y = new complex<double>[N];
    // bit reversal permutation
    for(int idx=0; idx<N; idx++){
        int tar_idx = bit_reversal(idx, bit_length);
        y[tar_idx] = static_cast<T>(x[idx]);
    }
    // O(nlogn)
    for(int bit_loc=1; bit_loc<=bit_length; bit_loc++){
        int n = 1<<bit_loc; // n = 2^{bit_loc}
        complex<double> w = exp(-1i*(2.0*PI/n));
        complex<double> wj = 1.0+0.0i;
        for (int k=0; k<n/2; k++) {
            for (int j=k; j<N; j+=n) {
                complex<double> even_term = y[j];
                complex<double> odd_term = y[j+n/2];
                y[j] = even_term + wj*odd_term;
                y[j+n/2] = even_term - wj*odd_term;
            }
            wj *= w;
        }
    }
    return y;
}

template<class T>
function<T*(T*, int)> interpolate(T const *known_x_arr, T const *known_y_arr, int N, string type){
    // linear interpolation
    auto linear_interp = [known_x_arr, known_y_arr, N](T* tar_x_arr, int n){
        T* tar_y_arr = new T[n];
        for(int idx=0; idx<n; idx++){
            // find lower bound index of known_x_arr which is closest to tar_x
            T tar_x = tar_x_arr[idx];
            int m;
            int l = 0;
            int h = N;
            while(l<h){
                m = (l+h)/2;
                if(tar_x > known_x_arr[m]){
                    // check known_x_arr[m+1:]
                    l=m+1;
                } else{
                    // check known_x_arr[:m]
                    h = m;
                }
            }
            // make sure that known_x_arr[l] <= tar_x
            if((l>0) && (known_x_arr[l]>tar_x)) {
               l--;
            }
            
            // interpolate
            tar_y_arr[idx] = known_y_arr[l] + (known_y_arr[l+1]-known_y_arr[l]) * (tar_x_arr[idx]-known_x_arr[l]) / (known_x_arr[l+1]-known_x_arr[l]);
            
        }
        return tar_y_arr;
    };
    // cubic spline interpolation's definition not yet
    auto cubic_interp = nullptr;
    
    if(type == "linear")
        return linear_interp;
    else if(type == "cubic")
        return cubic_interp;
    else
        throw logic_error("interpolation not implement error");
}

double* fft_call_pricer(int N, double eta, double alpha, double S0, double tau, double r, double d, double* strike_prices, int n_contract, double* chf_args){
    // price lots of options with different strike prices
    cout << "it's fast fourier transfrom call option pricer" << endl;
    
    // ========== initialization ==========
    double lambda = 2*PI / (N * eta);
    double beta = log(S0) - lambda * N / 2; // beta:= -b + log(S0)
    
    complex<double> *parameter = new complex<double>[N];
    
    for(int idx=0; idx<N; idx++){
        double v_j = eta*idx;
        parameter[idx] = v_j - (alpha+1)*1i;
    }
    complex<double> *result_cf = chf_black_scholes(parameter, N, S0, tau, r, d, chf_args[0]);
    delete[] parameter;
    
    // ========== calculation ==========
    complex<double> *x_arr = new complex<double>[N];
    // delta = 0, idx = 0. hence, v_j = eta*idx = 0.
    x_arr[0] = 1.0 * (exp(-r*tau) * result_cf[0] / (alpha*alpha + alpha)) * (eta/3.0) * (3.0 + pow(-1.0,1.0) - 1);
    // delta = 1
    for(int idx=1; idx<N; idx++){
        double v_j = eta*idx;
        complex<double> modified_call_price = exp(-r*tau) * result_cf[idx] / (alpha*alpha + alpha - v_j*v_j + 1i*(2*alpha+1)*v_j);
        x_arr[idx] = exp(-1i * beta * v_j) * modified_call_price * (eta/3.0) * (3.0 + pow(-1.0,idx+1.0) - 0);
    }
    
    complex<double>* fft_prices = fft(x_arr, log(N)/log(2)); // log_2(N) = ln(N)/ln(2)
    
    double* sim_c_arr = new double[N]; // sim_c_arr:= simulated call option price array
    double* sim_k_arr = new double[N]; // sim_k_arr:= simulated log strike price array
    for(int idx=0; idx<N; idx++){
        double log_strike = beta + lambda*idx;
        sim_k_arr[idx] = exp(log_strike);
        sim_c_arr[idx] = (exp(-alpha * log_strike) / PI) * fft_prices[idx].real();
    }
    
    // ========== interpolation ==========
    double* call_prices = nullptr;
    try{
        call_prices = interpolate(sim_k_arr, sim_c_arr, N, "linear")(strike_prices, n_contract);
    } catch(const exception& e){
        cerr << e.what() << endl;
    }
    return call_prices;
}
