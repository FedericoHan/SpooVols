import pandas as pd
import numpy as np
import datetime
from datetime import timedelta  #to shift datetime days up or down
import matplotlib.pyplot as plt
#%matplotlib inline
#import blpapi

import blpapi
from FXOption import FXOption
import DownloadData_v4



def bisection(f, a, b, tol = 0.1, maxiter = 10):
    """
    :param f: The function to solve
    :param a: The x-axis value where f(a)<0
    :param b: The x-axis value where f(b)>0
    :param tol: The precision of the solution
    :param maxiter: Maximum number of iterations
    :return: The x-axis value of the root,
                number of iterations used
    """
    c = (a+b)*0.5  # Declare c as the midpoint ab
    n = 1  # Start with 1 iteration
    while n <= maxiter:
        c = (a+b)*0.5
        if f(c) == 0 or abs(a-b)*0.5 < tol:
            # Root is found or is very close
            return c, n

        n += 1
        if f(c) < 0:
            a = c
        else: 
            b = c
            
    return c, n
    

#a = FXOption('Put', 4400, 4000, 0.2822, 0.015, 0.01, date(2022,6,17), 4400*50, 'SIF')
#print('root is ', root)
#print('iteratoins is ', iterations)

class ImpliedVolatilityModel(object):
    
    def __init__(self, S0, r, T, div, is_call = False):
        
        self.S0 = S0
        self.r  = r
        self.T = T
        self.div = div
        #self.N = N
        self.is_call = is_call
         
    def _option_valuation_(self, K , sigma):
        #use FX option  class
        option = FXOption('Put', self.S0, K, sigma, self.div, self.r, exp, 1, 
                           'DomesticPips')
        print(self.T)
        return option.price()
    
    def get_implied_volatilities(self, Ks, opt_prices):
        impvols = []
        for i in range(len(Ks)):
            f = lambda sigma: self._option_valuation_(Ks[i], sigma)- \
                opt_prices[i]
            print(f)
            impv = bisection(f, 0.01, 0.99, 0.0001, 100)[0]
            impvols.append(impv)
            print('strike'+str(Ks[i]), 'opt price is '+str(opt_prices[i]), impv)
        return impvols
    

#if __name__ == '__main__':
strikes= np.arange(2800, 5200, 100)
start = datetime.datetime(2022,2,11)
exp = datetime.date(2022,6, 18)
T = 120/365#(exp- start).days / 365.0


       

vals = []
for k in strikes:
    #this crap returns a list of dfs , shape (1,1) --> (date as index, price)
    vals.append(DownloadData_v4.DownloadData('ESM2P '+str(k)+ ' Index', 'PX_LAST',\
                                    datetime.datetime(2022,2,11),datetime.datetime(2022,2,11),'DAILY', 'blp').get_data_blp_historical(0))

vals_neater = []
for val in vals:
    vals_neater.append(val.iloc[0,0])
    
test2 = np.array(vals_neater)
model = ImpliedVolatilityModel(4400, 0.00125, T, 0.015,  is_call  = False)

impvols_put = model.get_implied_volatilities(strikes, test2)


plt.plot(strikes, impvols_put)
df = pd.DataFrame({'strikes': strikes,
                   'prems': test2,
                   'vols': impvols_put})