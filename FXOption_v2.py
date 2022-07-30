import math 
import numpy as np
from scipy.stats import norm
from datetime import date
import matplotlib.pyplot as plt
import numpy as np
   
#store common attributes of fx option
#these are split from market inputs like spot and vol, which appearing in the pricing and greek methods


class FXOption(object):   #for NonDeliverable need to overwrite outright forward
    '''
        For FXOptions make sure notional is in the LHS currency, as it's most commonly asked, e.g. AUD 100mio vs USD, 
        EUR 35mio vs USD, USD 50mio vs JPY
        
        premConvention:  
        
        %Foreign: default. This is the case for USDJPY, USDNJA etc.  You take the domestic amount (JPY), 
        divide by spot and LHS notional

        DomesticPips: this is the case in the brokers for AUDUSD, EURUSD.  LHS_notional x pips converts
        into domestic (USD) units,   as well as for options on Futures and Equities! So use this for all that stuff
        
        %Domestic: rare
        Foreign_Pips: very rare
        
        dateEnd format should be date(YYYY, MM, DD), NOT datetime.date() else cannot subtract  from date.today()
        
    ''' 

    
    def __init__(self, CallPutFlag, K , dateEnd, notional, premConvention):
        #date format: date(2017, 6, 30)  NOT datetime.date()
        
        #constructors 
        self.CallPutFlag = CallPutFlag
        self.K = K
        ###-------------------------------------
        # Convert date to unit of years
        ###--------------------------------
        
        #if option expiry is set in interger format (number of days)
        if isinstance(dateEnd, np.int64) or isinstance(dateEnd, int):
            self.Tau = dateEnd/365
            
        #if option expiry is set in date format 
        else:
            today = date.today()    
            self.dateEnd = dateEnd
            self.Tau = (self.dateEnd - today).days / 365.0
       
        if self.Tau <= 0:
            raise ValueError('Option Expired')
            
        #sanity check #print('Time is :'+str(self.Tau)) 
        
        self.notional = notional
        self.premConvention = premConvention
   
        '''
        currency_basis = 0.0
        self.f =  self.S * math.exp( ((self.rDom - self.rFor) +currency_basis)*self.Tau) 
        #print('f: '+str(self.f))
        
        pi = 3.141592654
        #density function of d1
        self.nd_1 = (1 / pow(2 * pi,0.5)) * math.exp(-self.d1 * self.d1 / 2)'''
    
    def price(self, S, sigma, rFor, rDom):

        #rename constructors for readibility 
        CallPutFlag, K, Tau, notional, premConvention = \
        self.CallPutFlag,  self.K, self.Tau, self.notional, self.premConvention
        
        ###-------------------------------------
        # d1 and d2 from standard BlackScholes
        ###--------------------------------
        d1 = (math.log(S/ self.K) + (rDom - rFor + 0.5*pow(sigma,2))*Tau)/(sigma*math.sqrt(self.Tau))
        d2 = d1 - sigma * math.sqrt(self.Tau)
    
        if CallPutFlag == "Call":

            dummy = S * math.exp(-rFor * self.Tau) * norm.cdf(d1) - K * math.exp(-rDom * 
                                           self.Tau)* norm.cdf(d2)    
        else:
            dummy = K * math.exp(-rDom * self.Tau) * norm.cdf(-d2) - S * math.exp(-rFor * \
                                          self.Tau)* norm.cdf(-d1)
        
        if premConvention == 'DomesticPips':
            premium = notional *dummy
        else:
            ## default setting:  USD10mio x 100yen/usd converted back to LHS ie USD
            premium = notional * (dummy / S)
            
        #return '{:,}'.format(premium)  #format method to have commas on long numbers
        return premium
    
    def delta(self, S, sigma, rFor, rDom, convention):
        
        '''self.S = S
        self.sigma = sigma
        self.rFor = rFor
        self.rDom = rDom'''
        '''
        SEF: spot delta Excluding premium in Foreign, this is the standard BS delta, tells how much foreign currency (AUD, EUR etc) 
        needed today to hedge, excluding premium
        SIF: includes premium, used when premium is in foreign units like USDJPY ..or AUDUSD with premium in AUD% instead of standard USD pips
        'like in the brokers
        SED: rare
        SID: rare
        FEF: normal blackscholes forward delta
        FIF: this is what we use for long-dated CrossYEN...BE CAREFUL with the forward as we got basis points now on top
        of usual InterestRates parity
        FED: rare
        FID: rare
        Notional should be in Foreign, LHS, currency, ie EUR...if not, convert by spot rate
        The extreme example of a deep in the money,  0.6000 AUD call / USD put: delta is 40ish% if premium is in AUD
        (as the buyer pays AUD 6mio or so, so the remaining risk is only on AUD 4mio) vs 100% delta if it were in USD. '''
        
        d1 = (math.log(S/ self.K) + (rDom - rFor + 0.5*pow(sigma,2))*self.Tau)/(sigma*math.sqrt(self.Tau))
        d2 = d1 - sigma * math.sqrt(self.Tau)
        #print(d1, d2) 
        dummy = 0 
        if convention == 'SEF':
            if self.CallPutFlag == 'Call':
                dummy = math.exp(-rFor * self.Tau) * self.notional * self.CND(d1)
                print(dummy)
            else:
                dummy = -math.exp(-rFor * self.Tau) * self.notional * self.CND(-d1)
                
        elif convention == 'SIF':
            if self.CallPutFlag == 'Call':

                dummy = math.exp(-rDom * self.Tau) * self.notional * self.K * self.CND(d2)/ S
            else:
                 dummy = -math.exp(-rDom * self.Tau) * self.notional * self.K * self.CND(-d2)/ S
                    
        elif convention == 'FEF':
            if self.CallPutFlag == 'Call':
                dummy = self.notional * self.CND(d1)
            else:
                dummy = -self.notional* self.CND(-d1)
                
        elif convention == 'FIF':
            if self.CallPutFlag == 'Call':
            #be careful your forward is correct!
               dummy = self.notional * self.K * self.CND(d2) / self.f
            else: 
                dummy = -self.notional * self.K * self.CND(-d2)/ self.f

        return dummy #'{:,}'.format(dummy)
    
    def gamma(self, S, sigma, rFor, rDom , convention):

        spot_bump_up = S * (1+ 0.05/100) #bump 1bp higher
        delta_bump_up = self.delta(spot_bump_up, sigma, rFor, rDom, convention)
       
        spot_bump_down = S * (1- 0.05/100) #bump 1bp higher
        delta_bump_down = self.delta(spot_bump_down, sigma, rFor, rDom, convention)

        numerator = (delta_bump_up - delta_bump_down)
        denominator = (spot_bump_up - spot_bump_down)
        
        gamma = (numerator / denominator) * S / 100 
        
        '''  alternatively second derivative of price:
        price_bump_up = self.price(spot_bump_up, sigma, rFor, rDom)
        price_no_bump = self.price(S, sigma, rFor, rDom)
        spot_bump_down = S * (1- 0.01/1000) #bump 1bp higher
        print(spot_bump_down)
        price_bump_down = self.price(spot_bump_down, sigma, rFor, rDom)
        denominator = (spot_bump_up - S)**2
        print(denominator)
        gamma2 = (price_bump_up + price_bump_down - 2*price_no_bump)/denominator * S/ 100'''
        
        return gamma # '{:,}'.format(gamma)

    def vega(self):
        b = self.rDom - self.rFor
        #vega = self.notional * (1/100) * (self.S * math.exp(-self.rFor * self.Tau) * pow(self.Tau, 0.5) * self.nd_1
        
        vega = self.notional * math.exp((b - self.rDom) * self.Tau) * self.nd_1 * pow(self.Tau, 0.5) / 100
        return '{:,}'.format(vega)
    
    '''#def vanna(self):
        
     #   vanna= -Exp(-riskfree_foreign * Tau) * nd_1 * (d2 / Sigma)
      #  return vanna
    
    #def volga(self):
        
     #   volga = S * Exp(-riskfree_foreign * Tau) * Sqr(Tau) * nd_1 * (d1 * d2 / Sigma)
      #  return volga'''
    
    def strike_from_delta(self, delta_percent, convention):
        
        '''For DeltaNeutralStraddle, we want delta_call + delta_put = 0, or delta_call = -delta_put
        ie, exp(-rf*Tau)*CND(d1) = --exp(-rf*Tau)*CND(-d1)---> CND(d1) = CND(-d1)---> d1 = -d1
        d1 = (log(F/K) + (sigma *sigma * 0.5* Tau ))/(sigma*sqrt(Tau) ) = 0
        log (F/K) = - 0.5 * sigma * sigma * Tau  --> K = F * exp(0.5 * sigma * sigma * Tau)...premium unadjusted like EURUSD
        K =  F* exp(-0.5 *sigma *sigma*Tau) ...premium adjusted like USDXYZ'''
        
        delta_target = delta_percent
        increment = 0.0001
        strike_guess = self.f #initial guess
        debugger = 0
        k_SEF = 0 #initialize 

        if delta_percent == "DNS":
            if convention == "SEF":
                dummy = self.f * math.exp(0.5 * self.sigma *self.sigma * self.Tau)  #premium unadjusted
            else:
                dummy = self * math.exp(-0.5 * self.sigma *self.sigma * self.Tau) #premium adjusted
            return dummy
    
       #Below finds strikes for non-straddles. If EXCLUDING premium, we can easily invert the Delta Function
   
        if convention == 'SEF':
            if self.CallPutFlag == 'Call':
                k_SEF = -self.sigma * pow(self.Tau, 0.5) * norm.ppf(delta_percent / math.exp(-self.rFor  * self.Tau)) + \
                0.5 * self.sigma * self.sigma * self.Tau
            else:
                delta_percent = -1.0 * delta_percent #ensure it's a put, delta is negative
                k_SEF = self.sigma * pow(self.Tau, 0.5) * norm.ppf(-delta_percent / math.exp(-self.rFor  * self.Tau)) + \
                0.5 * self.sigma * self.sigma * self.Tau
            dummy = self.f * math.exp(k_SEF)
            return dummy
        
        #If INCLUDING premium, we need a numerical procedure to find the strike. use below only for <1y, else need forward 
        if convention == 'SIF':
            if self.CallPutFlag == 'Call':
                
                delta_guess = 0.5
                
                while delta_guess > delta_target:  #vs do until delta_guess < delta_target
                    
                    d1_guess = (math.log(self.S/ strike_guess) + (self.rDom - self.rFor + 0.5*pow(self.sigma,2))*self.Tau) \
                    /(self.sigma*math.sqrt(self.Tau))
                    d2_guess = d1_guess - self.sigma * math.sqrt(self.Tau)
                    
                    delta_guess = math.exp(-self.rDom * self.Tau) * strike_guess * self.CND(d2_guess) /self.S
                
                    strike_guess += increment 
                    #print(strike_guess, delta_guess)
                    debugger += 1
                    
                         
            else: #for puts
                delta_guess = -0.5
                delta_target = -1 *delta_target
                while delta_guess < delta_target:
                    
                    d1_guess = (math.log(self.S/ strike_guess) + (self.rDom - self.rFor + 0.5*pow(self.sigma,2))*self.Tau) \
                    /(self.sigma*math.sqrt(self.Tau))
                    d2_guess = d1_guess - self.sigma * math.sqrt(self.Tau)
                    
                    delta_guess = math.exp(-self.rDom * self.Tau) * -1 *strike_guess * self.CND(-d2_guess)/ self.S
                    strike_guess -= increment
                    debugger += 1
            
            dummy = strike_guess
                    
            return dummy
                                 
    def CND(self, x):
        
        '''generally x will be d1 or d2
        CND is the cumulative normal distribution of a random variable, ie, the integral of the density function_
        from negative infinity up to the random variable. Describes probability will be found at a value <= to x.
        it's the NORMSDIST function in excel
        for negative values, like -1.96 it will be asymptotic to 0 (lower probabilities), at 0 it's 0.50 (fifty/fifty, _
        for large quantities (like > 1.96) asymptotic to 1 (certain probability).
        Graphed ti's the shape ofcall-spread. '''
    
        pi = 3.141592654
        l = abs(x)
        k = 1 / (1 + 0.2316419 * l)
        a1 = 0.31938153
        a2 = -0.356563782
        a3 = 1.781477937
        a4 = -1.821255978
        a5 = 1.330274429
        
        n = 1 / math.sqrt(2 * pi) * math.exp(-l * l / 2)
        
        CND = 1 - n * (a1 * k + \
                       a2 * pow(k, 2) + \
                       a3 * pow(k, 3) + \
                       a4 * pow(k, 4) + \
                       a5 * pow(k, 5))
        if x < 0:
            CND = 1 - CND
        return CND

    def graph(self, dS):
 
        'generate value given S'
        CallPutFlag, S, K, sigma, rFor, rDom, Tau, notional, premConvention = \
        self.CallPutFlag, self.S, self.K, self.sigma, self.rFor, self.rDom,\
        self.Tau, self.notional, self.premConvention
        
        #dS = 0.001
        path = []
        values = [] 
        #upper_bound = int(round(S/dS)+1)
        upper_bound = int(S * 1.15)
        for i in range(upper_bound):
            path.append(S)
            opt_price = FXOption(CallPutFlag, S, K, sigma,  rFor, rDom, self.dateEnd, notional,\
                                 premConvention).price()
            values.append(opt_price)
            S = S + dS
#http://localhost:8888/notebooks/FXOption_class.ipynb#
        #print(values)
        for s, payoff in zip(path, values):
            # print ('%.2f  %.4f' % (s, payoff_values)) #TypeError: a float is required
            print('%.2f  %.4f' %(s, payoff))
    