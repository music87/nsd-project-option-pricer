import pytest
import numpy as np
import timeit
import time
from fast_fourier_transfrom_call_pricer import fftp
class TestClass:
    def test_price(self):
        # hyper-parameters
        N = 2**20
        eta = 0.01
        alpha = 1.0
        
        # option contract
        r = 0.06
        tau = 0.4
        d = 0.01
        sigma = 0.2
        
        # market information
        S0 = 50.0
        
        # set up
        engine = fftp(N, eta, alpha)
        engine.set_chf_param(sigma)
        engine.set_contract(S0, tau, r, d)
        # engine.print_setting()
        
        # pricing
        strike_prices = [5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 65.0]
        gold_option_prices = [45.0, 40.0, 30.0, 20.0, 10.01, 0, 0] 
        calc_option_prices = engine.calc_price(strike_prices)
        
        assert np.all(np.abs((np.array(gold_option_prices) - np.array(calc_option_prices))) < 1)
    def test_time(self):
        setup='''
from fast_fourier_transfrom_call_pricer import fftp
from fftoptionlib.process_class import BlackScholes
from fftoptionlib.helper import spline_fitting
from fftoptionlib.fourier_pricer import carr_madan_fft_call_pricer
import numpy as np
N = 1048576
eta = 0.01
alpha = 1.0

# option contract
r = 0.06
tau = 0.4
d = 0.01
sigma = 0.2

# market information
S0 = 50.0

# golden
sim_k_arr, sim_c_arr = carr_madan_fft_call_pricer(N, eta, alpha, r, tau, S0, d, BlackScholes(sigma).set_type('chf'))
ffn_pricer = spline_fitting(sim_k_arr, sim_c_arr, 3)

# set up
engine = fftp(N, eta, alpha)
engine.set_chf_param(sigma)
engine.set_contract(S0, tau, r, d)
'''
        repeat = 5
        tcalc = timeit.Timer('calc_option_prices = engine.calc_price(list(np.arange(0,10000000)))', setup=setup)
        calc_sec = min(tcalc.repeat(repeat=repeat, number=1))
	
        tgold = timeit.Timer('gold_prices = ffn_pricer(np.arange(0,10000000))', setup=setup)
        gold_sec = min(tcalc.repeat(repeat=repeat, number=1))
        
        print(f"self-defined function takes {calc_sec} secounds")
        print(f"golden function takes {gold_sec} secounds")
        assert ((calc_sec/gold_sec) < 1)

