# See https://sites.google.com/site/applegatearchive/software/python-model-fitting

###################################
# Utilities for fitting models
# 
# Based on solution by abeardmore found on http://code.google.com/p/pyminuit/issues/detail?id=6
#
# Modified and extended by Douglas Applegate
###################################

import numpy
import minuit
import math, inspect
import scipy.stats as stats
import numpy as np

###############################

__cvs_id__ = "$Id: fitmodel.py,v 1.2 2010-07-02 23:08:47 dapple Exp $"

###############################

###############################
# Statistical Distribution Look-up functions
###############################

def chisq_exceeds_prob(chisq, dof):
    '''
    Probability that chisq exceeds value, given degrees of freedom dof
    '''
    return stats.chi2.sf(chisq, dof)

###

def f_test_exceeds_prob(chisq_old, dof_old, chisq_new, dof_new):
    '''
    Probability that the improvement in a fit by adding extra parameters is random
    '''

    deltaDOF = dof_old - dof_new

    F = (chisq_old - chisq_new)/(deltaDOF*chisq_new/dof_new)

    return stats.f.sf(F, deltaDOF, dof_new)

###############################
# Common Models
###############################

def ConstantModel(x, a0):
    
    return a0

#######

def LinearModel(x, a0, a1):

    return a0 + a1*x

########

def QuadraticModel(x, a0, a1, a2):

    return a0 + a1*x + a2*x**2

########

class PolynomialModel(object):
    """
    Creates a polynomial model of the form
    a0 + a1*x + a2*x**2 + ... 
    where the order parameter controls which orders are included
    """

    def __init__(self, order):
        '''
        order is a list of positive integers specifying polynomial order to include
        0: constant, 1: linear, 2: quadratic, etc.
        Does not include lower order terms implicitly (ie specify [0,1,2], etc
        '''

        self.order = order

        self.basis = {}
        for o in order:
            param = 'a%d' % o
            def base(x, a, order=o):
                return a*(x**order)
            self.basis[param] = base

        self.params=self.basis.keys()

    def __call__(self, x, *params, **keyword_params):

        for key, val in zip(self.params, params):
            keyword_params[key] = val
	
        sum = 0.
        for key, val in keyword_params.iteritems():
            sum += self.basis[key](x, val)

        return sum
		
def FringeModel(x, P0, P1, P2, P3):
	
	K = 46.3059
	a0 = 7.4636846162628681
	an = np.array([  8.30746496e-01,   6.07104132e-01,   4.09709205e-01,
         2.59069834e-01,   1.50693753e-01,   7.43633937e-02,
         2.83097700e-02,  -2.21975469e-03,  -1.87698388e-02,
        -2.47989956e-02,  -2.35531438e-02,  -2.06760491e-02,
        -1.66221750e-02,  -1.22393197e-02,  -1.02114038e-02,
        -7.51736988e-03,  -5.69295385e-03,  -3.19839702e-03,
        -1.58987130e-03,  -1.47961862e-03,   8.71519206e-05,
        -5.99893671e-04,  -2.61315191e-04,   5.51034636e-05,
        -6.89555948e-04,  -8.32234569e-04,  -1.87782953e-04,
         1.81491780e-03,   5.20437031e-03,   4.77993879e-03,
         7.89033683e-03,   9.94192145e-03,   1.18143615e-02,
         1.98280552e-02,   2.90215432e-02,   4.04890778e-02,
         5.59149190e-02,   6.94247417e-02,   8.02758507e-02,
         8.27707278e-02,   7.88008065e-02,   5.30145204e-02,
        -3.52681929e-03,  -1.02587410e-01,  -2.62361423e-01])
	bn = np.array([  1.71304710e-01,   2.52290171e-01,   2.77629497e-01,
         2.63994694e-01,   2.32952157e-01,   1.90851079e-01,
         1.48741865e-01,   1.13435660e-01,   8.06010666e-02,
         5.48605872e-02,   3.39137017e-02,   2.03681459e-02,
         1.13235426e-02,   4.32165070e-03,   3.72824870e-03,
         3.56841239e-03,   9.78768315e-04,   4.18165104e-03,
         9.65436338e-04,  -1.80753154e-03,  -1.10939094e-03,
        -1.61726916e-03,  -1.91963947e-04,  -3.19529958e-04,
        -1.58112130e-03,  -7.38734958e-04,  -2.32834256e-03,
        -4.07224082e-04,   7.94026775e-04,  -3.24343894e-03,
        -2.64329371e-03,  -4.39668096e-03,  -5.37573007e-03,
        -3.37157863e-03,   3.32274059e-04,   8.09470952e-03,
         2.23203056e-02,   4.50449112e-02,   8.01770794e-02,
         1.26785615e-01,   1.89062110e-01,   2.72332104e-01,
         3.69862054e-01,   4.84164841e-01,   6.02820240e-01])
	
	fringe = np.zeros(np.size(x))
	
	for i in range(0,np.size(x)):
		sumN = 0
		U2p = P1*a0
		for n in range(0,np.size(an)):
			# Note n is 1 indexed below (e.g., n+1) because the 1st FC is in an[0], bn[0]
			sumN = sumN + (an[n]*np.cos(2*math.pi*(n+1)/K*(x[i]-K*P2)) +
						   bn[n]*np.sin(2*math.pi*(n+1)/K*(x[i]-K*P2))) * np.exp(-((n+1)*P3)**2)
			
		fringe[i] = P0 + (U2p*(1+sumN))#*(x[i]*P[4]+P[5]))
	
	# Return the modeled fringe
	return fringe

def PowerLawModel(x, alpha, beta):

    return alpha*x**beta

###########

def GaussianModel(x, A, mu, sigma):

    z = (x - mu) / sigma
    return A*numpy.exp(-0.5*z**2)

###############################
# Statistical Fuctions for Minimization
###############################

def ChiSqStat(ydata, yerr, ymodel):
    """
    Returns the chi-square given arrays of ydata, yerr, and ymodel values.
    """
    chisquared = ((ydata - ymodel)/yerr)**2 
    stat = chisquared.sum()
    return stat

####################

def CStat(ydata, yerr, ymodel):
    """
    Returns the cstat a la xspec given arrays of data and model values.
    This is a -2.0 log likelihood statistic.
    """

    lmodel = numpy.zeros(ymodel.size)

    lmodel[ymodel <= 0.0] = -32.

    lmodel[ymodel > 0.0] = numpy.log(ymodel[ymodel > 0.0])

    ldata = numpy.zeros(ydata.size)

    ldata[ydata <= 0.0] = -32.0

    ldata[ydata > 0.0] = numpy.log(ydata[ydata > 0.0])

    # fitstat = ymodel - ydata  + ydata * (ldata - lmodel)

    fitstat = ymodel + ydata  * ((ldata - lmodel) - 1.0)

    stat = 2.0* fitstat.sum()

    return stat



###############################
# Fitting Class -- Use to perform minimizations
###############################

class FitModel:
    """
    Fits a generic model (provided by the class Model to data (numpy arrays
    xdata and ydata), with a fit statistic provided by StatFunc.
    """
    def __init__(self, xdata, ydata, yerr, model, 
                 statfunc = ChiSqStat, guess = []):

        self.xdata = numpy.array(xdata, dtype=numpy.float64)
        self.ydata = numpy.array(ydata, dtype=numpy.float64)
        self.yerr = numpy.array(yerr, dtype=numpy.float64)

        self.model = model
        
        self.statfunc = statfunc

        self.guess = guess

        self.fcn = FCN(self.xdata, self.ydata, self.yerr, model, statfunc)

        self.m = minuit.Minuit( self.fcn )

        self.params = self.m.parameters

        if self.guess == []:
            self.guess = numpy.ones(len(self.params))

        for param, value in zip(self.params, self.guess):
            self.m.values[param] = value
            self.m.errors[param] = math.fabs(value) * 0.05

        self.m.strategy = 1
        self.m.tol = 1.0

        self.have_fit = False

    def fixed(self, fparams):
        """
        Fix or unfix the parameters specified in the dictionary fparams, which
        contain True or False values.
        """
        for key in fparams.keys():
            self.m.fixed[key] = fparams[key]

    def limits(self, lparams):
        """
        Set limits given by the parameters in the dictionary lparams.
        """
        for key in lparams.keys():
            self.m.limits[key] = lparams[key]

    def fit(self, printmode = 0):
        """
        Call migrad to fit the model to the data.
        Set printmode = 1 to monitor the progress of the fitting.
        """

        self.m.printMode = printmode

        self.par_vals = {}

        self.ymodel = None

        try :

            self.m.migrad()

            print("fval = %g, nfcn %d" % (self.m.fval, self.m.ncalls))

            self.m.migrad()

            print("fval = %g, nfcn %d" % (self.m.fval, self.m.ncalls))

            print("Fit parameters : ")
            print(self.m.values)

            self.par_vals = self.m.values

            # calculate the best fit model
            self.ymodel = self.model( self.xdata, **self.m.values )

            self.statval = self.m.fval

            self.have_fit = True

        except minuit.MinuitError :

            # reset have_fit if migrad fails
            self.have_fit = False



    def uncert(self, nsigma = 1.0):
        """
        Calculate the parameter uncertainties at the nsigma**2
        confidence level. E.g. for one parameter of interest
        nsigma = 1.0   for 68%
                 1.645 for 90%
                 2.0   for 95.45%
                 3.0   for 99.73%
        """

        if not(self.have_fit) :
            print("Warning: uncert requires a valid fit.")
            return

        # in case minos fails
        self.m.hesse()

        print("Hesse errors : ")
        print(self.m.errors)

        self.par_err = {}

        for key in self.m.values.keys():

            if (self.m.fixed[key] == True):
                continue

            try:

                self.m.minos(key, -nsigma)
                self.m.minos(key,  nsigma)


                error = (self.m.merrors[key, -nsigma], 
                         self.m.merrors[key,  nsigma])

            except minuit.MinuitError :

                print("Caught MinuitError: Minos failed. using Hesse error.")
                print("Only really valid for a well behaved fitting FCN !")
                error = self.m.errors[key] * nsigma
    

            self.par_err[key] = error

            
        print("Parameter errors :")
        print(self.par_err)



    def corr_matrix(self):
        """
        Display the fit parameter correlation matrix."
        """
        if not(self.have_fit) :
            print("Warning: uncert requires a valid fit.")
            return

        print("Correlation matrix :")
        print(numpy.array(self.m.matrix(correlation=True)))

    



#####################################
# Utilities
####################################

    
def FCN(x,y,yerr, model, statfunc):
    """
    Calculates the fitting FCN for pyMinuit(2) given the data (xdata & ydata)
    and model (class Model, with a tuple of initial parameters, params),
    using the class StatFunc to calculate the statistic.
    """

    #assumes model is a function with first arg being X values
    if inspect.isfunction(model):
        params = inspect.getargspec(model)[0][1:]
    elif hasattr(model, '__call__'):
        args = inspect.getargspec(model.__call__)[0]
        if len(args) < 3:
            paramAttr = inspect.getargspec(model.__call__)[1]
            params = getattr(model, paramAttr)
        else:
            params = args[2:]

    paramstring = ','.join(params)

    class_template = '''class fitclass(object):
    def __init__(self, x, y, yerr, model, statfunc):
        self.x = x
        self.y = y
        self.yerr = yerr
        self.model = model
        self.statfunc = statfunc

    def __call__(self, %s):
        return self.statfunc(self.y, self.yerr, self.model(self.x, %s))
''' % (paramstring, paramstring)

    
    exec(class_template)


    return fitclass(x,y,yerr,model,statfunc)


