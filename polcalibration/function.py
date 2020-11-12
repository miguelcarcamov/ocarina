import numpy as np
from scipy.optimize import curve_fit
import abc

class Function(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self, xdata=[], x_0=0.0):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.coeff = []

    def getCoeffs(self):
        return self.coeff

    def getxdata(self):
        return self.xdata

    def getX_0(self):
        return self.x_0

    def setCoeffs(self, coeff):
        self.coeff = coeff

    def setxdata(self, xdata):
        self.xdata = xdata

    def setX_0(self, x_0):
        self.x_0 = 0

    @abc.abstractmethod
    def f(self, xdata, *args):
        return

    @abc.abstractmethod
    def f_eval(self, xdata, coeffs):
        return

    def fit(self, xdata, data, initial_coeffs=None, lowerbound=-np.inf, upperbound=np.inf):
        lowerbounds = np.ones(len(initial_coeffs))*lowerbound
        upperbounds = np.ones(len(initial_coeffs))*upperbound
        popt, pcov = curve_fit(self.f, xdata, data, p0=initial_coeffs, check_finite=True, bounds=(lowerbounds,upperbounds))
        self.coeff = popt
        return popt, pcov

class FluxFunction(Function):
    def __init__(self, flux_0=0.0, **kwargs):
        super(FluxFunction, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

    def f(self, xdata, *args):
        nu_div = xdata/self.x_0
        exp = args[0] + args[1] * np.log(nu_div)
        s_flux = self.flux_0 * (nu_div)**exp
        return s_flux

    def f_eval(self, xdata, coeffs):
        nu_div = xdata/self.x_0
        exp = coeffs[0] + coeffs[1] * np.log(nu_div)
        s_flux = self.flux_0 * (nu_div)**exp
        return s_flux

class PolFunction(Function):
    def __init__(self, nterms=0, **kwargs):
        super(PolFunction, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        self.coeff = []

    def getNTerms(self):
        return self.nterms

    def f(self, xdata, *args):
        y = np.zeros(len(xdata))
        for i in range(0, self.nterms):
            y += args[i] * np.power((xdata - self.x_0) / self.x_0, i)
        return y

    def f_eval(self, xdata, coeffs):
        y = np.zeros(len(xdata))
        for i in range(0, self.nterms):
            y += coeffs[i] * np.power((xdata - self.x_0) / self.x_0, i)
        return y
