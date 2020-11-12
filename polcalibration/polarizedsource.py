import numpy as np
import os
from __casac__.table import table as tb
from utils import degreesToRadians, radiansToDegrees, queryTable
from function import Function, FluxFunction, PolFunction

# Object that takes information of different known polarized sources from https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol
class PolarizedSource(object):
    def __init__(self, nu=np.array([]), polangle=np.array([]), polfrac=np.array([]), name="", knownSource="", **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.name = ""
        self.spix_dict = {'Perley-Butler 2013' : 'PerleyButler2013Coeffs', 'Perley-Butler 2017' : 'PerleyButler2017Coeffs', 'Scaife-Heald 2012' : 'ScaifeHeald2012Coeffs'}
        switcher = {"3c48" : "p3c48", "3c48_2019": "p3c48_2019", "3c138" : "p3c138", "3c138_2019" : "p3c138_2019", "3c286": "p3c286", "3c286_2019": "p3c286_2019", "3c147" : "p3c147", "3c147_2019": "p3c147_2019"}
        if self.knownSource != "":
            method = getattr(self, switcher.get(self.knownSource, "init_empty"))
            method()

        # convert frequencies to Hz
        self.nu *= 1e9

        # Pol angle from degrees to radians
        self.polangle = degreesToRadians(self.polangle)

        # Pol fraction in percentage to fraction
        self.polfrac /= 100.0
        self.spidx_coeffs = []

    def p3c48(self):
        self.nu = np.array([1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5])
        self.polangle = np.array([25.0, 140.0, -5.0, -150.0, -120.0, -100.0, -92.0, -84.0, -75.0, -72.0, -68.0, -67.0, -64.0, -62.0, -62.0, -62.0, -63.0, -64.0, -66.0, -67.0, -70.0, -70.0, -77.0, -85.0])
        self.polfrac = np.array([0.3, 0.5, 0.7, 0.9, 1.4, 2.0, 2.5, 3.2, 3.8, 4.2, 5.2, 5.2, 5.3, 5.4, 6.0, 6.1, 6.4, 6.4, 6.9, 7.1, 7.7, 7.8, 7.4, 7.5])
        self.name = "3C48"

    def p3c48_2019(self):
        self.nu = np.array(np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1]))
        self.polangle = np.array([4.3, -34.0, 23.0, 67.1, -84.0, -72.0, -66.0, -63.0, -62.0, -63.0, -64.0, -68.0, -72.0, -76.0, -77.0, -84.0, -84.0])
        self.polfrac = np.array([0.3, 0.5, 0.9, 1.6, 2.9, 4.3, 5.4, 5.4, 5.7, 6.1, 6.3, 6.5, 7.2, 6.4, 6.7, 5.6, 6.8])
        self.name = "3C48"

    def p3c138(self):
        self.nu = np.array([1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5])
        self.polangle = np.array([-14.0, -11.0, -10.0, -10.0, -9.0, -10.0, -10.0, 0.0, -11.0, -11.0, -12.0, -12.0, -10.0, -8.0, -7.0, -7.0, -8.0, -9.0, -12.0, -13.0, -16.0, -17.0, -24.0, -27.0])
        self.polfrac = np.array([5.6, 7.5, 8.4, 9.0, 10.4, 10.7, 10.0, 0.0, 10.0, 10.4, 9.8, 10.0, 10.4, 10.1, 8.4, 7.9, 7.7, 7.4, 6.7, 6.5, 6.7, 6.6, 6.6, 6.5])
        self.name = "3C138"

    def p3c138_2019(self):
        self.nu = np.array(np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1]))
        self.polangle = np.array([-13, -9.6, -9.3, -10., -9.5, -10.5, -11.5, -9.4, -7.9, -11, -13, -16, -18, -19, -20, -23, -24])
        self.polfrac = np.array([5.5, 7.8, 9.0, 9.9, 10.3, 10.5, 10.2, 10.9, 9.1, 8.2, 8.2, 8.4, 8.4, 8.5, 8.7, 8.8, 9.2])
        self.name = "3C138"

    def p3c286(self):
        self.nu = np.array([1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5])
        self.polangle = np.array([33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 35.0, 35.0, 35.0, 36.0, 36.0])
        self.polfrac = np.array([8.6, 9.5, 9.9, 10.1, 10.5, 10.8, 10.9, 11.1, 11.3, 11.4, 11.6, 11.7, 11.9, 11.9, 11.9, 11.9, 12.1, 12.2, 12.5, 12.5, 12.6, 12.6, 13.1, 13.2])
        self.name = "3C286"

    def p3c286_2019(self):
        self.nu = np.array(np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1]))
        self.polangle = np.array([33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 34.0, 34.0, 35.0, 35.0, 36.0, 36.0, 36.0, 37.0, 36.0])
        self.polfrac = np.array([8.6, 9.8, 10.1, 10.6, 11.2, 11.5, 11.9, 12.1, 12.3, 12.3, 12.5, 12.6, 12.7, 13.1, 13.5, 13.4, 14.6])
        self.name = "3C286"

    def p3c147(self):
        self.nu = np.array([4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5])
        self.polangle = np.array([-100.0, 0.0, -65.0, -39.0, -24.0, -11.0, 43.0, 48.0, 53.0, 59.0, 67.0, 68.0, 75.0, 76.0, 85.0, 86.0])
        self.polfrac = np.array([0.1, 0.3, 0.3, 0.6, 0.7, 0.8, 2.2, 2.4, 2.7, 2.9, 3.4, 3.5, 3.8, 3.8, 4.4, 5.2])
        self.name = "3C147"

    def p3c147_2019(self):
        self.nu = np.array([4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1])
        self.polangle = np.array([-13.0, -57.0, -19.0, 27.0, 53.0, 60.0, 66.0, 79.0, 83.0, 87.0, 87.0, 85.0])
        self.polfrac = np.array([0.16, 0.51, 0.48, 0.85, 1.8, 2.4, 2.9, 3.4, 4.0, 4.5, 4.9, 6.0])
        self.name = "3C147"

    def init_empty(self):
        self.nu = np.array([])
        self.polangle = np.array([])
        self.polfrac = np.array([])
        self.name = self.knownSource

    def getName(self):
        return self.name

    def setName(self, name=""):
        self.name = name

    def getNu(self):
        return self.nu

    def setNu(self, nu=np.array([])):
        self.nu = nu

    def getPolAngle(self):
        return self.polangle

    def setPolDeg(self, polangle=np.array([])):
        self.polangle = polangle

    def getPolFrac(self):
        return self.polfrac

    def setPolFrac(self, polfrac=np.array([])):
        self.polfrac = polfrac

    def isCorrect(self):
        if(len(self.nu) == len(self.polangle) and len(self.nu) == len(self.polfrac)):
            return True
        else:
            return False

    # Returns values in an array lower than a certain frequency
    def filter(self, nu, data, nu_min=0.0, nu_max=np.inf):
        indexes = np.where((nu > nu_min) & (nu < nu_max))
        nu = nu[indexes]
        data = data[indexes]
        return nu, data

    def flux_scalar(self, nu):
        flux_at_nu = 0.0

        for i in range(len(self.spidx_coeffs)):
            flux_at_nu += self.spidx_coeffs[i][0] * (np.log10(nu)**i)

        return 10**flux_at_nu

    def flux(self, nu):
        flux_at_nu = np.zeros(len(nu))

        for i in range(len(self.spidx_coeffs)):
            flux_at_nu += self.spidx_coeffs[i][0] * (np.log10(nu)**i)

        return 10**flux_at_nu

    def getCoeffs(self, standard="Perley-Butler 2013", epoch="2012"):
        coeff_table = os.getenv('CASAPATH').split(' ')[0] + '/data/nrao/VLA/standards/' + self.spix_dict[standard]
        tb_obj = tb()
        tb_obj.open(coeff_table)
        query_table = tb_obj.taql("select * from "+coeff_table+" where Epoch="+epoch)
        coeffs = query_table.getcol(self.name+"_coeffs")
        self.spidx_coeffs = coeffs
        tb_obj.close()
        return coeffs

    def fitAlphaBeta(self, nu, nu_0=0.0):
        if(not nu_0):
            nu_0 = np.median(nu)
        flux_0 = self.flux_scalar(nu_0/1e9)
        fluxes = self.flux(nu/1e9)

        i_coeffs = np.random.rand(2)
        source_func_frac = FluxFunction(flux_0=flux_0, xdata=nu, x_0=nu_0)
        source_func_frac.fit(nu, fluxes, i_coeffs)
        return source_func_frac.getCoeffs()

    def getPolFracCoeffs(self, nu_0=0.0, nterms=3, nu_min=0.0, nu_max=np.inf):
        nu, polfrac = self.filter(self.getNu(), self.getPolFrac(), nu_min, nu_max)
        if(not nu_0):
            nu_0 = np.median(nu)
        ifrac_coeffs = np.random.uniform(0.0, 1.0, nterms)
        source_func_frac = PolFunction(x_0=nu_0, nterms=nterms)
        source_func_frac.fit(nu, polfrac, ifrac_coeffs)
        return source_func_frac.getCoeffs()

    # Returns pol angle coeffs in radians
    def getPolAngleCoeffs(self, nu_0=0.0, nterms=3, nu_min=0.0, nu_max=np.inf):
        nu, polangle = self.filter(self.getNu(), self.getPolAngle(), nu_min, nu_max)
        if(not nu_0):
            nu_0 = np.median(nu)
        iangle_coeffs = np.random.uniform(-np.pi, np.pi, nterms)
        source_func_angle = PolFunction(x_0=nu_0, nterms=nterms)
        source_func_angle.fit(nu, polangle, iangle_coeffs)
        return source_func_angle.getCoeffs()
