#pragma once
#include <functional>
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
std::function<vector<T>(vector<T>)> interpolate(T const *known_x_arr, T const *known_y_arr, int N, string type){
    // linear interpolation
    auto linear_interp = [known_x_arr, known_y_arr, N](vector<T> tar_x_arr){
        vector<T> tar_y_arr;
        for(size_t idx=0; idx<tar_x_arr.size(); idx++){
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
            tar_y_arr.push_back(known_y_arr[l] + (known_y_arr[l+1]-known_y_arr[l]) * (tar_x_arr[idx]-known_x_arr[l]) / (known_x_arr[l+1]-known_x_arr[l]));
            
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

