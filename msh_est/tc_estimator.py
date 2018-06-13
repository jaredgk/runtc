"""
    tc_estimator.py  - copyright Jody Hey 2018
    trimmed down version of AlleleAgeEstimator.py

    requires math,scipy

    Provides estimates for tc, the  time to the bottom of an edge carrying a mutation shared by one or more sampled chromosomes.
        can use likelihood or bayesian estimators

        for k=1,   options are dobayes and nosnp
            dobayes = False is the default  ** bayesian estimate not recommended

            set nosnp=True if the data come from a base position on a chromosome
                that does not carry a derived allele at that position
                nosnp = false is the default

        for k>1 gene copes, options are dobayes and estmethod  (nosnp does not apply for k>1)
            dobayes = False is the default
            there are 3 estimators:
            'a' : average, mean of k maximum likelihood or bayesian estimates
            'b' : product of exponentials taken before integration over phi (pex method)
            'c' : composite likelihood (default)

    three classes:

        fitmodel()
            used to specify a population size model and the estimation procedure
            create at least one of these before creating data or datalist

        data()
            create a data object
            each data object has a model object
            estimation functions belong to data objects
                call: estimate_tc(self,dobayes = None):

        datalist()
            create a datalist object
            a special kind of list of data objects
            all data objects must have the same model
            estimation functions are run over all data objects in the data list:
                    estimate_tc(self,dobayes = None):


    For k==1 use a single instance of data, e.g. d the main function is d.estimate_tc()

    For k > 1 use a datalist, e.g. dl, the main function is dl.estimate_tc()

    Examples at bottom of this file (comment out to run as a module)

"""
import math
import scipy.optimize as opt
import scipy.integrate as integrate
import numpy as np
import sys

from scipy.special import gammainc
from scipy.special import gamma

largepos = 1e100 ##use for invalid likelihood values during optimization
fact=[1]  ## used by gamma function

def gammafunc(*arg):
    """
        Calculates the gamma function, including factorials for positive integers.
        Also this function can calculate what Mathematica calls  Gamma[a,x]  which is the Gamma integration from x to infinity

        first element of arg is value a
        if there is a second element it is x

        For Mathematica Gamma[a,x] the integration is done from x to infinity.
        Some other sources call this 'incomplete gamma', there is a fair bit of confusion
        Scipy has a gammainc function but this returns what is sometimes called a regularized incomplete gamma
        The conversion, to use scipy in a function that returns the same as what Mathematica does with Gamma[a,x] is explained here:
        see http://stackoverflow.com/questions/38713199/incomplete-gamma-function-in-scipy

        to use scipy gammainc, a must be positive and x must be >= 0.
        if x=0 and a is an integer this function just returns (a-1)!
    """
    def factorial(a):
        """
            return the factorial of an integer
            build up array fact on as needed basis to save time
        """
        global fact  ## must be initialized somewhere else with fact = [1]
##        if not isinstance(a,int):
##            raise TypeError ('a is not an integer')
        while len(fact)<=a:
            fact.append(len(fact)*fact[-1])
        return fact[a]  ## fact[a] has the factorial for number a

    a = arg[0]
    if len(arg)>1:
        x = arg[1]
    else:
        x = 0
    if a<= 0 or x < 0:
        print("failure of gammafunc ",a,x)
        exit()
    if isinstance(a,int) and a>= 1:
        if x > 0:
            return factorial(a-1) * (1-gammainc(a,x))
        else:
            return factorial(a-1)
    else:
        return gamma(a)* (1-gammainc(a,x))

def roundChi(f,n):
    """
    Rounds chi to 14 digits after decimal, which fixes some caching errors
    """
    if n < 1:
        raise Exception("N value %d is not valid" % (n))
    return round(f,-int(math.floor(math.log10(f)))+(n-1))

class fitmodel:
    """
        makes a structure to hold information about a population size model, and optimization parameters

        can run one of three population size model types
            model c - constant population size
            model e - exponentially growing population
            model h - two phase history, with constant population size followed by exponential growth until time of sampling

        demography parameters:
            n : sample size
            N0 : final population size
            g : population growth rate
            te : time when exponential growth began

        All models require n and N0

        The model type can be specified when created, if it is not then
            if neither growth rate nor te are given, model c (constant) is invoked
            If growth rate g is given, then model e (exponential) or h (two phase) is invoked
            if growth rate and te are given , then model h (two phase is invoked)


        Optimization uses Brent, but user can specify optmodel = "Golden" to try  (usually it gives the same answer)

        User can specify the bracket intervals for the initial part of optimization
        not sure how much this matters.
            lbb - lowerbound of starting bracket
            ubb - upperbound of starting bracket
            defaults are 10,100 for models c and h and 1,maxtc for model e

        A second bracket is also defined for use with datalists.  In this case optimization is tried with both brackets and the best results reported.

    """
    def __init__(self, n=None, N0=None, g=None, te=None, popmodel=None, lbb = None,ubb = None, optmethod = "Brent", nosnp = None):
        MINTIME = 0.01
        self.n = n
        self.N0 = N0
        self.g = g
        self.te = te
        self.nosnp = nosnp == True
        assert n > 1 and N0 > 1
        if popmodel==None:
            if g == None:
                if te != None:
                    print("te is used without g")
                    exit()
                self.popmodel = "c"
                self.maxtc = None
                self.mintc = MINTIME
                self.zeropoint = None
                self.bracket = [self.mintc,100]
                self.bracket2 = [100,1000]
            else:
                assert g > 0
                if te == None:
                    self.popmodel = "e"
                    self.mintc = MINTIME  ##  4 seems to be necessary when using node2
                    self.maxtc = math.log(N0)/g
                    self.zeropoint = None
                    self.bracket = [self.mintc,self.maxtc/2]
                    self.bracket2 = [self.maxtc/2,self.maxtc]
                else:
                    assert te > 0
                    self.popmodel = "h"
                    self.mintc = MINTIME
                    self.maxtc = None
                    self.zeropoint = (math.exp(-(g*te))*(n - 4*g*N0 - n*math.exp(g*te) + g*n*te*math.exp(g*te)))/(g*n)
##                    self.bracket = [self.mintc,100]
##                    self.bracket2 = [100,1000]
                    self.bracket = [self.mintc,min(self.zeropoint/2,self.mintc*10)]
                    self.bracket2 = [5*self.zeropoint,10000]
        else:
            self.popmodel = popmodel
            if popmodel=="c":
                self.maxtc = None
                self.mintc = MINTIME
                self.zeropoint = None
                self.bracket = [self.mintc,100]
                self.bracket2 = [100,1000]
            if popmodel == "e":
                self.maxtc = math.log(N0)/g
                self.mintc = MINTIME ## 4seems to be necessary for node2
                self.zeropoint = None
                self.bracket = [self.mintc,self.maxtc/2]
                self.bracket2 = [self.maxtc/2,self.maxtc]
            if popmodel == "h":
                self.maxtc = None
                self.mintc = MINTIME
                try:
                    self.zeropoint = te + (-1 + math.exp(-(g*te)))/g - (4*N0* math.exp(-(g*te)))/n
##                    self.zeropoint = (math.exp(-(g*te))*(n - 4*g*N0 - n*math.exp(g*te) + g*n*te*math.exp(g*te)))/(g*n)
                except Exception as ex:
                    print("n,No,g,Te:",n,N0,g,te)
                    print(ex)
                    exit()
                self.bracket = [self.mintc,min(self.zeropoint/2,self.mintc*10)]
                self.bracket2 = [5*self.zeropoint,10000]
        if n==None or N0 == None:
            print("n or N0 is none")
            exit()
        if lbb != None:
            self.bracket[0] = lbb
            self.bracket2[0] = lbb
        if ubb != None:
            self.bracket[1] = ubb
            self.bracket[2] = ubb
        self.optmethod = optmethod


class data:
    """
        user creates an instance of data using physical and genetic map distances for one or both sides of the focal snp

        side 1 and side 2 are arbitrary.  If data comes from only one side, it is side 1

        dis1 : physical distance to side 1, REQUIRED
        dis2 : physical distance to side 2
        mu : mutation rate per base pair per generation (e.g. 1.25e-8), REQUIRED
        rho1 :  recombination rate per base pair per generation at last base on side 1, REQUIRED for tc estimates (not required for ta estimates)
        rho2 :  recombination rate per base pair per generation at last base on side 2
        morgans1 : genetic map distance (in Morgans, not centimorgans) on side 1
        morgans2 : genetic map distance on side 2

        if dis2 == None, then estimate is done for side 1
        if morgans1== None  (morgans2==None)  then map distance to side 1 is rho1*dis1

        The mutation rate is constant on both sides.

        e.g.
        d = data(dis1=103100,dis2=421113,mu=1e-8,rho1 = 1e-8)   two distances, constant recombination rate, use the same mu and rho
        d = data(dis1=103100,mu=1e-8,rho1 = 1e-8)   one distance
        d = data(dis1 = 1000000, dis2=1000000, rho1 = 1e-8, rho2 = 4.5e-8, mu = 1e-8, morgans1 = 0.02,morgans2 = 0.1)
            two distances, each with variable recombination rates

        To generate an estimate after an instance has been created :
            data.estimate_tc (model, dobayes = None)
                where model is an instance of fitmodel()
                dobayes is a boolean for whether the prior should be used or not (default is False)
    """

    def __init__(self,dis1 = None, dis2=None, rho1 = None, rho2 = None, mu = None, morgans1 = None, morgans2 = None, model = None, side = None, chi = None):
        """
            must have either   rho1,morgans1,dis1 all defined, or rho2,morgans2,dis2 all defined, or both
        """
        if model != None:
            if not isinstance(model, fitmodel):
                raise TypeError('item is not of type fitmodel')
            self.model = model
        else:
            self.model = None
        self.rho1 = rho1
        self.rho2 = rho2
        self.morgans1 = morgans1
        self.morgans2 = morgans2
        self.dis1 = dis1
        self.dis2 = dis2
        if mu == None:
            mu = 0.0  # do not use mutation distance in the method
        self.mu = mu
        if chi is not None:
            self.chi = chi
            self.side = side
            if side == 3:
                self.singlex = False
            else:
                self.singlex = True
            return
        if self.morgans1 == None and self.rho1 != None and  self.dis1 != None:
            self.morgans1 = self.rho1 * self.dis1
        if self.morgans2 == None and self.rho2 != None and  self.dis2 != None:
            self.morgans2 = self.rho2 * self.dis2
        if (self.morgans1 == None and self.dis1 != None) or (self.morgans1 != None and self.dis1 == None):
            print("morgans1 and dis1 problem, can't have just one be None")
            exit()
        if (self.morgans2 == None and self.dis2 != None) or (self.morgans2 != None and self.dis2 == None):
            print("morgans2 and dis2 problem, can't have just one be None")
            exit()

        if self.morgans1 == None or self.dis1 == None:  ## data only from side 2
            #assert self.morgans2 != None and self.dis2 != None and self.morgans2 > 0 and self.dis2 > 0
            self.side = 2
            self.singlex = True
            self.chi = (mu * (self.dis2) + self.morgans2 )
        else:
            if self.morgans2 == None or self.dis2 == None: ## data only from side 1
                #assert self.morgans1 != None and self.dis1 != None and self.morgans1 > 0 and self.dis1 > 0
                self.side = 1
                self.singlex = True
                self.chi = (mu * (self.dis1) + self.morgans1 )
            else:
                self.side = 3
                self.singlex = False
                self.chi = (mu * (self.dis1 + self.dis2) + self.morgans1 + self.morgans2)

    def setmodel(self,model):
        if not isinstance(model, fitmodel):
            raise TypeError('item is not of type fitmodel')
        self.model = model

    def calcprior(self,tc):
        n = self.model.n
        N0 = self.model.N0
        if self.model.popmodel == "h":
            te = self.model.te
        if self.model.popmodel == "c":
            try:
                self.logprior = math.log((4*n)/(N0*(2 + (n*tc)/(2.*N0))**3))
            except Exception as ex:
                self.logprior = 0.0
        elif self.model.popmodel =="e" or ( self.model.popmodel =="h"  and tc < te ):
            g = self.model.g
            self.logprior = math.log( (4*n* math.exp(g*tc))/(N0*(2 + (n*(-1 + math.exp(g*tc)))/(2*g*N0))**3) )
        else:
            g = self.model.g
            self.logprior = math.log( (4*n* math.exp(g*te))/(N0*(2 + n*((-1 + math.exp(g*te))/(2 *g*N0) + ((tc - te)* math.exp(g*te))/(2 * N0)))**3) )
        return self.logprior

    def tc_likelihood_neg(self,tc, *args):
        if self.model == None:
            print("no model specified for data")
            exit()
        try:
            dobayes = args[0]
        except Exception as ex:
            print("problem with dobayes = args[0]", args)
            print(ex)
            exit()
        try:
            if len(args) > 1:
                nosnp = args[1]
            else:
                nosnp = False
        except Exception as ex:
            print("problem with nosnp = args[1]", args)
            print(ex)
            exit()
        if tc < self.model.mintc or (self.model.maxtc != None and tc > self.model.maxtc):
            return largepos
        if self.model.zeropoint != None and abs(tc - self.model.zeropoint)< 0.001:
            tc = self.model.zeropoint + 0.001
        mu = self.mu
        chi = self.chi
        n = self.model.n
        N0 = self.model.N0
        if self.model.popmodel == "h":
            te = self.model.te
        if self.model.popmodel=="c":## or model == "f":
            if self.singlex:
                try:
                    temp = (8 * N0 * tc *  math.exp(-2 * tc * chi))/(4 * N0 + n * tc) + (-(n *  math.exp(-2 * tc * chi)) - 2 * n * tc * chi * math.exp(-2 * tc * chi) + n *  (1 + tc * chi) * math.exp(-(tc * chi)))/((4 * N0 + n * tc) * chi**2)
                except Exception as ex:
                    return largepos
            else:
                try:
                    temp = (16 * N0 * tc**2 *  math.exp(-2 * tc * chi))/(4 * N0 + n * tc) + (-2 * n *  math.exp(-2 * tc * chi) - 4 * n * tc * chi * math.exp(-2 * tc * chi) - 4 * n * tc**2 *  chi**2 * math.exp(-2 * tc * chi) + 2 * n *  math.exp(-(tc * chi)) + 2 * n * tc * chi * math.exp(-(tc * chi)) + n * tc**2 * chi**2 * math.exp(-(tc * chi)))/((4 * N0 + n * tc) * chi**3)
                except Exception as ex:
                    return largepos
        elif self.model.popmodel =="e" or ( self.model.popmodel =="h"  and tc < te ):
            maxtc = self.model.maxtc
            g = self.model.g
            if (tc*(g +  chi)) > 500:
                return largepos
            if self.singlex:
                try:
                    temp= (8*g*N0*tc*  math.exp(-2*tc* chi ))/(4*g*N0 + n*(-1 +  math.exp(g*tc))) + (g*n*  math.exp(-2*tc* chi )*(-1 - 2*tc*(g +  chi ) + (1 + tc*(g +  chi ))* math.exp(tc*(g +  chi ))))/((g +  chi )**2*(4*g*N0 + n*(-1 +  math.exp(g*tc))))
                except Exception as ex:
                    return largepos
            else:
                try:
                    temp1 = (16*g*N0*tc**2*  math.exp(-2*tc* chi))/(4*g*N0 + n*(-1 +  math.exp(g*tc)))
                    temp2 = (g*n* math.exp(-2*tc* chi)*(-2 - 4*tc*(g +  chi)*(1 + g*tc + tc* chi) + (2 + tc*(g +  chi)*(2 + tc*(g +  chi)))* math.exp(tc*(g +  chi))))/((g +  chi)**3*(4*g*N0 + n*(-1 +  math.exp(g*tc))))
                    temp = temp1 + temp2
                except Exception as ex:
                    return largepos
        else:  ## model "h" and tc > te
            assert self.model.popmodel=="h" and tc >= te
            g = self.model.g
##            if  (te*(g +  chi)) > 500:
##                return largepos
            if self.singlex:
                try:
                    temp = (g* math.exp(-3*tc* chi )*(8*N0*tc* math.exp(tc* chi ) + (n*((-1 - 2*tc* chi  + te* chi )* math.exp(tc* chi ) + (1 + tc* chi )* math.exp((2*tc - te)* chi ))* math.exp(te*(g +  chi )))/ chi **2 + (n* math.exp(tc* chi )*(-1 - 2*tc*(g +  chi ) + (1 + (2*tc - te)*(g +  chi ))* math.exp(te*(g +  chi ))))/(g +  chi )**2))/(4*g*N0 + n*(-1 + (1 + g*tc - g*te)* math.exp(g*te)))
                except Exception as ex:
                    return largepos
            else:
                try:
                    temp1 = (g*n*  math.exp(g*te)*((2 + tc* chi *(2 + tc* chi ))* math.exp(-(tc* chi )) - (2 + (2*tc - te)* chi *(2 + 2*tc* chi  - te* chi ))* math.exp((-2*tc + te)* chi )))/( chi **3*(4*g*N0 + n*(-1 + (1 + g*tc - g*te)* math.exp(g*te))))
                    temp2 = (g*n* math.exp(-2*tc* chi )*(-2 - 4*tc*(g +  chi )*(1 + g*tc + tc* chi ) + (2 + (2*tc - te)*(g +  chi )*(2 + (2*tc - te)*(g +  chi )))* math.exp(te*(g +  chi ))))/((g +  chi )**3*(4*g*N0 + n*(-1 + (1 + g*tc - g*te)* math.exp(g*te))))
                    temp3 = (16*g*N0*tc**2*  math.exp(-2*tc* chi ))/(4*g*N0 + n*(-1 + (1 + g*tc - g*te)* math.exp(g*te)))
                    temp = temp1 +  temp2 + temp3
                except Exception as ex:
                    return largepos
        if temp > 0.0 and temp < float("inf"):
            if self.model.nosnp:
                mterm = 0.0
            else:
                mterm = math.log ( 1 - math.exp(-(tc*mu)))
            l = math.log(temp) + mterm
        else:
            return largepos
        if dobayes == True:
            self.calcprior(tc)
            l += self.logprior
        return -l

    def estimate_tc(self,dobayes = None, nosnp = None):
        """
            estimate_tc for k=1
        """
        if nosnp:
            args = (dobayes,nosnp)
        else:
            args = (dobayes,)
        rval = opt.minimize_scalar(fun=self.tc_likelihood_neg,args=args,method=self.model.optmethod,bracket=self.model.bracket)
        self.tcest = rval.x
        self.loglike = -rval.fun


class datalist(list):
    """
        a class that is a list, plus some other things
        to create an instance, a = datalist()
        can be filled only with instances of type data
        e.g.
            a = datalist()   to create an empty datalist
            d = data(...)     to create an instance of data
            a.append(d)     to put things in the datalist

        To generate an estimate of tc after an instance has been created :

            estimates can be generated using a composite likelihood estimator (default)
            or a mean likelihood estimator

            estimate_tc (self,model, dobayes = None, estmethod = None)
                where model is an instance of fitmodel()
                dobayes is a boolean for whether the prior should be used or not (default is to not use it)
                estmethod sets the type of estimator:
                    estmethod = 'a' : use mean of likelihood estimates
                    estmethod = 'c' : use composite likelihood (default)

        To generate an estimate of ta :
            estimate_tmrca (self,mod,k,dobayes = None)
                where model is an instance of fitmodel()
                k is the number of chromosomes that share the mutation
                    Note, there should be k(k-1)/2 data points in the list

    """
    def __init__(self):
        self.cache = {}
        self.cache_single = {}
        self.bins = []
        self.bins_single = []
        self.estimates_bin = []
        self.estimates_bin_single = []
        self.cache_hits = 0
        self.cache_total = 0

    def append(self, item):
        if not isinstance(item, data):
            raise TypeError('item is not of type data')
        super(datalist, self).append(item) ## not sure why this works, but it overrides append
        if len(self) > 1:
            assert ( self[-1].model == None and self[-2].model == None) or self[-1].model.popmodel == self[-2].model.popmodel

    def setmodel(model):
        """
            change the model for all the data elements in self
        """
        if not isinstance(model, fitmodel):
            raise TypeError('item is not of type fitmodel')
        for dval in self:
            dval.setmodel(model)

    def tc_likelihood_composite_neg(self,tc, args):
        """
            calculate the composite likelihood
            if dobayes, calculate the prior just once
        """
        sum = 0.0
        dobayes = args
        args = False
        for dval in self:
            sum += dval.tc_likelihood_neg(tc, args)
        if dobayes == True:
            self[0].calcprior(tc)
            sum -= self[0].logprior  ## subtract because tc_likelihood_neg passes the negative of the likelihood
##        if sum==float("inf") or sum==float("-inf"):
##            print "inf"
        return sum


    def estimate_tc(self,dobayes = None,estmethod = None):
        """
            estimate_tc for k>1,  arguments are dobayes and estmethod

            unlike estimate_tc for k==1,  there is no 'nosnp' argument

            estmethod can be
            'a' : average, mean of k maximum likelihood or bayesian estimates
            'b' : product of exponentials taken before integration over phi (pex method)
            'c' : composite likelihood (default)
        """
        args = (dobayes,)
        optmethod = self[0].model.optmethod
        bracket = self[0].model.bracket
        bracket2 = self[0].model.bracket2
        if estmethod == None or estmethod == 'c':
            bracket = self[0].model.bracket
            rval1 = opt.minimize_scalar(fun=self.tc_likelihood_composite_neg,args=args,method=optmethod,bracket=bracket)
            rval2 = opt.minimize_scalar(fun=self.tc_likelihood_composite_neg,args=args,method=optmethod,bracket=bracket2)
            if  rval1.fun < rval2.fun or rval2.fun == float("inf") or rval2.fun == float("-inf"):
                self.tcest = rval1.x
                self.loglike = -rval1.fun
            else:
                self.tcest = rval2.x
                self.loglike = -rval2.fun
        elif estmethod=='a':
            sumest = 0.0
            sumlike = 0.0
            for dval in self:
                dval.estimate_tc(dobayes)
                sumest += dval.tcest
                sumlike += dval.loglike
            self.tcest = sumest/len(self)
            self.loglike = sumlike/len(self)
        else:
            print("estimation method ",estmethod," not recognized")
            exit()

    #Cache options: basic, bin, round, none
    def estimate_tc_cache(self,cache=True,round=-1,bin=False):
        chi = self[0].chi
        if round != -1:
            chi = roundChi(chi,round)
        fill_cache = False
        if cache:
            if self[0].singlex and str(chi) in self.cache_single:
                self.tcest = self.cache_single[str(chi)]
                self.cache_hits += 1
            elif not self[0].singlex and str(chi) in self.cache:
                self.tcest = self.cache[str(chi)]
                self.cache_hits += 1
            else:
                fill_cache = True
            self.cache_total += 1
        if bin:
            self.estimate_tc_bin()
        else:
            self.estimate_tc(dobayes=False)
        if fill_cache:
            if self[0].singlex:
                self.cache_single[str(chi)] = self.tcest
            else:
                self.cache[str(chi)] = self.tcest

    def calc_tc_bins(self,start,stop,mu,model,num=10000):
        self.bins = list(np.geomspace(start,stop,num))
        self.bins_single = list(np.geomspace(start,stop,num))
        for val in self.bins:
            d = data(mu=mu,model=model,chi=val,side=3)
            d.estimate_tc(dobayes=False)
            self.estimates_bin.append(d.tcest)
        for val in self.bins_single:
            d = data(mu=mu,model=model,chi=val,side=2)
            d.estimate_tc(dobayes=False)
            self.estimates_bin_single.append(d.tcest)

    def estimate_tc_bin(self):
        chi = self[0].chi
        idx = self.getBinIdx(chi)
        self.tcest = self.getIntEst(chi,idx,self[0].singlex)


    def getBinIdx(self,chi):
        if len(self.bins) == 0:
            raise Exception("Bins have not been initiated")
        if chi < self.bins[0]:
            raise Exception("Chi value %f is lower than minimum in bins %f" % (chi,self.bins[0]))
        if chi > self.bins[-1]:
            raise Exception("Chi value %f is higher than maximum in bins %f" % (chi,self.bins[-1]))
        min_idx = 0
        max_idx = len(self.bins)-2
        while True:
            if max_idx < min_idx:
                raise Exception("Error finding region for chi")
            m = (min_idx+max_idx)//2
            if self.bins[m] <= chi < self.bins[m+1]:
                return m
            elif self.bins[m] > chi:
                max_idx = m-1
            else:
                min_idx = m + 1


    def getIntEst(self,chi,idx,singlex):
        s = chi - self.bins[idx]
        s /= (self.bins[idx+1]-self.bins[idx])
        if singlex:
            return s*self.estimates_bin_single[idx]+(1-s)*self.estimates_bin_single[idx+1]
        return s*self.estimates_bin[idx]+(1-s)*self.estimates_bin[idx+1]


##some examples  comment out when using as a module

##print ("\nk=1. Constant population model")
##m = fitmodel(n=7100,N0 = 10000)
##
##d = data(dis1 = 10000, dis2=10000, rho1 = 1.2e-8, rho2 = 1e-8, morgans1 = 0.0002, morgans2 = 0.0004,mu = 1e-8,model = m)
##d.estimate_tc()
##print ("tc likelihood estimate ", d.tcest)
##d.estimate_tc(nosnp = True)
##print ("tc likelihood estimate, no snp ", d.tcest)
##
##print ("\nk=1 (same data point).  Exponential growth model")
##m = fitmodel(n=7100,N0 = 552808,g = 0.03119)
##d.setmodel(m)
##d.estimate_tc()
##print ("tc likelihood estimate ", d.tcest)
##d.estimate_tc(nosnp = True)
##print ("tc likelihood estimate, no snp ", d.tcest)
##
##
##print ("\nk=1 (same data point)  Two phase population model")
##m = fitmodel(n=7100,N0 = 552808,g = 0.03119,te = 129.5,popmodel="h")
##d.setmodel(m)
##d.estimate_tc()
##print ("tc likelihood estimate ", d.tcest)
##d.estimate_tc(nosnp = True)
##print ("tc likelihood estimate, no snp ", d.tcest)
##
##print ("\nk=3, exponential growth model")
##m = fitmodel(n=7100,N0 = 552808,g = 0.03119)
##dl = datalist()
##d = data(dis1 = 10000000, dis2=100000, rho1 = 1.2e-8, rho2 = 1e-8, morgans1 = 0.02, morgans2 = 0.04,mu = 1e-8,model = m)
##dl.append(d)
##d = data(dis1 = 170000, dis2=90000, rho1 = 1e-8, rho2 = 2e-8, morgans1 = 0.003, morgans2 = 0.005,mu = 1e-8,model = m)
##dl.append(d)
##d = data(dis1 = 50000, dis2=230000, rho1 = 3e-8, rho2 = 2.2e-8, morgans1 = 0.0035, morgans2 = 0.0052,mu = 1e-8,model = m)
##dl.append(d)
##
##dl.estimate_tc(estmethod = 'a')
##print ("tc average likelihood estimate ", dl.tcest)
##dl.estimate_tc(estmethod = 'c')
##print ("tc composite likelihood estimate ", dl.tcest)
