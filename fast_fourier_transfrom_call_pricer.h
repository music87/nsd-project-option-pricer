//#pragma once
#include "global_variables.h"
#include "characteristic_functions.h"
#include "utils.h"
#include <iostream>
#include <complex>
#include <string>
#include <stdexcept>
#include <vector>
#include <array>
/*#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>*/
using namespace std;
class FFT_call_option_pricer{
public:
    FFT_call_option_pricer(int N, double eta, double alpha);
    void set_contract(double S0, double tau, double r, double d);
    void set_black_scholes_chf_param(double sigma);
    
    void print_setting();
    
    vector<double> calc_price(vector<double> strike_prices);
    
    
private:
    // hyper-parameters
    int dN;
    double deta;
    double dalpha;
    double dsigma;
    
    // option contract
    double dS0;
    double dtau;
    double dr;
    double dd;
};

FFT_call_option_pricer::FFT_call_option_pricer(int N, double eta, double alpha){
    dN=N;
    deta=eta;
    dalpha=alpha;
    dsigma = -1.0;
    dS0 = -1.0;
    dtau = -1.0;
    dr = -1.0;
    dd = -1.0;
}

void FFT_call_option_pricer::set_contract(double S0, double tau, double r, double d){
    dS0 = S0;
    dtau = tau;
    dr = r;
    dd = d;
}

void FFT_call_option_pricer::set_black_scholes_chf_param(double sigma){
    dsigma = sigma;
}

void FFT_call_option_pricer::print_setting(){
    cout << "========= hyper-parameters for fft =========" << endl;
    cout << "N: " << dN << endl;
    cout << "eta: " << deta << endl;
    cout << "alpha: " << dalpha << endl;
    cout << "========= hyper-parameters for chf =========" << endl;
    cout << "sigma: " << dsigma << endl;
    cout << "========= contract =========" << endl;
    cout << "underlying asset price: " << dS0 << endl;
    cout << "time to maturity: " << dtau << endl;
    cout << "risk free rate: " << dr << endl;
    cout << "dividend rate: " << dd << endl;
}

vector<double> FFT_call_option_pricer::calc_price(vector<double> strike_prices){
    // price lots of options with different strike prices
    
    // ========== initialization ==========
    double lambda = 2*PI / (dN * deta);
    double beta = log(dS0) - lambda * dN / 2; // beta:= -b + log(S0)
    
    complex<double> *parameter = new complex<double>[dN];
    
    for(int idx=0; idx<dN; idx++){
        double v_j = deta*idx;
        parameter[idx] = v_j - (dalpha+1)*1i;
    }
    complex<double> *result_cf = chf_black_scholes(parameter, dN, dS0, dtau, dr, dd, dsigma);
    delete[] parameter;
    
    // ========== calculation ==========
    complex<double> *x_arr = new complex<double>[dN];
    // delta = 0, idx = 0. hence, v_j = eta*idx = 0.
    x_arr[0] = 1.0 * (exp(-dr*dtau) * result_cf[0] / (dalpha*dalpha + dalpha)) * (deta/3.0) * (3.0 + pow(-1.0,1.0) - 1);
    // delta = 1
    for(int idx=1; idx<dN; idx++){
        double v_j = deta*idx;
        complex<double> modified_call_price = exp(-dr*dtau) * result_cf[idx] / (dalpha*dalpha + dalpha - v_j*v_j + 1i*(2*dalpha+1)*v_j);
        x_arr[idx] = exp(-1i * beta * v_j) * modified_call_price * (deta/3.0) * (3.0 + pow(-1.0,idx+1.0) - 0);
    }
    
    complex<double>* fft_prices = fft(x_arr, log(dN)/log(2)); // log_2(N) = ln(N)/ln(2)
    
    double* sim_c_arr = new double[dN]; // sim_c_arr:= simulated call option price array
    double* sim_k_arr = new double[dN]; // sim_k_arr:= simulated strike price array
    for(int idx=0; idx<dN; idx++){
        double log_strike = beta + lambda*idx;
        sim_k_arr[idx] = exp(log_strike);
        sim_c_arr[idx] = (exp(-dalpha * log_strike) / PI) * fft_prices[idx].real();
    }
    
    // ========== interpolation ==========
    vector<double> call_prices;
    try{
        call_prices = interpolate(sim_k_arr, sim_c_arr, dN, "linear")(strike_prices);
    } catch(const exception& e){
        cerr << e.what() << endl;
    }
    return call_prices;
}

/*PYBIND11_MODULE(fast_fourier_transfrom_call_pricer, m){
    m.doc() = "compute the call option price over a series of strike prices";
    pybind11::class_<FFT_call_option_pricer>(m,"fftp")
        .def(pybind11::init<int, double, double>())
        .def("set_contract", &FFT_call_option_pricer::set_contract)
        .def("set_chf_param", &FFT_call_option_pricer::set_black_scholes_chf_param)
        .def("calc_price", &FFT_call_option_pricer::calc_price)
        .def("print_setting", &FFT_call_option_pricer::print_setting);
}*/
