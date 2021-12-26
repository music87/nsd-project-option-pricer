import pytest
import numpy as np
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
        strike_prices = [5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 65.0]
        
        # market information
        S0 = 50.0
        
        # ground truth
        gold_option_prices = [45.0, 40.0, 30.0, 20.0, 10.01, 0, 0] 

        engine = fftp(N, eta, alpha)
        engine.set_chf_param(sigma)
        engine.set_contract(S0, tau, r, d)
        # engine.print_setting()
        
        # pricing
        calc_option_prices = engine.calc_price(strike_prices)
        
        assert np.all(np.abs((np.array(gold_option_prices) - np.array(calc_option_prices))) < 1)

