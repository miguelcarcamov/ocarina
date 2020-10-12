from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True

def degreesToRadians(self, degrees):
    return degrees * np.pi / 180.0

def radiansToDegrees(self, radians):
    return radians * 180.0 / np.pi

class PolarizedSource(object):
    def __init__(self, nu=np.array([]), poldeg=np.array([]), polfrac=np.array([]), knownSource="", **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)

        switcher = {"3c48" : "p3c48", "3c48_2019": "p3c48_2019", "3c138" : "p3c138", "3c138_2019" : "p3c138_2019", "3c286": "p3c286", "3c286_2019": "p3c286_2019", "3c147" : "p3c147", "3c147_2019": "p3c147_2019"}
        if self.knownSource != "":
            method = getattr(self, switcher.get(self.knownSource, "init_empty"))
            method()


    def p3c48(self):
        self.nu = np.array([1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5]) * 1e9
        self.poldeg = np.array([25.0, 140.0, -5.0, -150.0, -120.0, -100.0, -92.0, -84.0, -75.0, -72.0, -68.0, -67.0, -64.0, -62.0, -62.0, -62.0, -63.0, -64.0, -66.0, -67.0, -70.0, -70.0, -77.0, -85.0])
        self.polfrac = np.array([0.3, 0.5, 0.7, 0.9, 1.4, 2.0, 2.5, 3.2, 3.8, 4.2, 5.2, 5.2, 5.3, 5.4, 6.0, 6.1, 6.4, 6.4, 6.9, 7.1, 7.7, 7.8, 7.4, 7.5])

    def p3c48_2019(self):
        self.nu = np.array(np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1])) * 1e9
        self.poldeg = np.array([4.3, -34, 23, 67.1, -84.0, -72, -66, -63, -62, -63, -64, -68, -72, -76, -77, -84, -84])
        self.polfrac = np.array([0.3, 0.5, 0.9, 1.6, 2.9, 4.3, 5.4, 5.4, 5.7, 6.1, 6.3, 6.5, 7.2, 6.4, 6.7, 5.6, 6.8])

    def p3c138(self):
        self.nu = np.array([1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5]) * 1e9
        self.poldeg = np.array([-14.0, -11.0, -10.0, -10.0, -9.0, -10.0, -10.0, 0.0, -11.0, -11.0, -12.0, -12.0, -10.0, -8.0, -7.0, -7.0, -8.0, -9.0, -12.0, -13.0, -16.0, -17.0, -24.0, -27.0])
        self.polfrac = np.array([5.6, 7.5, 8.4, 9.0, 10.4, 10.7, 10.0, 0.0, 10.0, 10.4, 9.8, 10.0, 10.4, 10.1, 8.4, 7.9, 7.7, 7.4, 6.7, 6.5, 6.7, 6.6, 6.6, 6.5])

    def p3c138_2019(self):
        self.nu = np.array(np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1])) * 1e9
        self.poldeg = np.array([-13, -9.6, -9.3, -10., -9.5, -10.5, -11.5, -9.4, -7.9, -11, -13, -16, -18, -19, -20, -23, -24])
        self.polfrac = np.array([5.5, 7.8, 9.0, 9.9, 10.3, 10.5, 10.2, 10.9, 9.1, 8.2, 8.2, 8.4, 8.4, 8.5, 8.7, 8.8, 9.2])

    def p3c286(self):
        self.nu = np.array([1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5]) * 1e9
        self.poldeg = np.array([33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 35.0, 35.0, 35.0, 36.0, 36.0])
        self.polfrac = np.array([8.6, 9.5, 9.9, 10.1, 10.5, 10.8, 10.9, 11.1, 11.3, 11.4, 11.6, 11.7, 11.9, 11.9, 11.9, 11.9, 12.1, 12.2, 12.5, 12.5, 12.6, 12.6, 13.1, 13.2])

    def p3c286_2019(self):
        self.nu = np.array(np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1])) * 1e9
        self.poldeg = np.array([33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 34, 34, 35, 35, 36, 36, 36, 37, 36])
        self.polfrac = np.array([8.6, 9.8, 10.1, 10.6, 11.2, 11.5, 11.9, 12.1, 12.3, 12.3, 12.5, 12.6, 12.7, 13.1, 13.5, 13.4, 14.6])

    def p3c147(self):
        self.nu = np.array([4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5]) * 1e9
        self.poldeg = np.array([-100.0, 0.0, -65.0, -39.0, -24.0, -11.0, 43.0, 48.0, 53.0, 59.0, 67.0, 68.0, 75.0, 76.0, 85.0, 86.0])
        self.polfrac = np.array([0.1, 0.3, 0.3, 0.6, 0.7, 0.8, 2.2, 2.4, 2.7, 2.9, 3.4, 3.5, 3.8, 3.8, 4.4, 5.2])

    def p3c147_2019(self):
        self.nu = np.array([4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1]) * 1e9
        self.poldeg = np.array([-13, -57, -19, 27, 53, 60, 66, 79, 83, 87, 87, 85])
        self.polfrac = np.array([0.16, 0.51, 0.48, 0.85, 1.8, 2.4, 2.9, 3.4, 4.0, 4.5, 4.9, 6.0])

    def init_empty(self):
        self.nu = np.array([])
        self.poldeg = np.array([])
        self.polfrac = np.array([])

    def getNu(self):
        return self.nu

    def setNu(self, nu=np.array([])):
        self.nu = nu

    def getPolDeg(self):
        return self.poldeg

    def setPolDeg(self, poldeg=np.array([])):
        self.poldeg = poldeg

    def getPolFrac(self):
        return self.polfrac

    def setPolFrac(self, polfrac=np.array([])):
        self.polfrac = polfrac

    def isCorrect(self):
        if(len(self.nu) == len(self.poldeg) and len(self.nu) == len(self.polfrac)):
            return True
        else:
            return False

class PolFunction(object):
    def __init__(self, nu_0=0.0, terms=0, **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)
        self.coeff = []

    def f(self, nu, *args):
        y = np.zeros(len(nu))
        for i in range(0, self.terms + 1):
            y += args[i] * np.power((nu - self.nu_0) / self.nu_0, i)
        return y

    def f_eval(self, nu, coeff):
        y = np.zeros(len(nu))
        for i in range(0, self.terms + 1):
            y += coeff[i] * np.power((nu - self.nu_0) / self.nu_0, i)
        return y

    def fit(self, nu, data, initial_coeffs):
        popt, pcov = curve_fit(self.f, nu, data, p0=initial_coeffs, check_finite=True)
        return popt, pcov


if __name__ == '__main__':
    degree_deg = 11
    degree_frac = 11
    nu_0 = 3.0 * 1e9
    #degree_deg = 2
    #degree_frac = 2
    #nu_0 = 20.0 * 1e9

    # Create a known polarized source
    # Options 3c48, 3c48_2019, 3c138, 3c138_2019, 3c286, 3c286_2019, 3c147, 3c147_2019
    knownSource = '3c286'
    pol_source = PolarizedSource(knownSource=knownSource)
    print pol_source.isCorrect()

    poldeg_func = PolFunction(nu_0=nu_0, terms=degree_deg)
    polfrac_func = PolFunction(nu_0=nu_0, terms=degree_frac)

    deg_coeff = np.random.uniform(-180.0, 180.0, degree_deg + 1)
    frac_coeff =np.random.rand(degree_frac + 1)

    # Do the fitting for polarization degree and fraction
    deg_coeff_popt, deg_coeff_pcov = poldeg_func.fit(pol_source.getNu(), pol_source.getPolDeg(), initial_coeffs=deg_coeff)
    frac_coeff_popt, frac_coeff_pcov = polfrac_func.fit(pol_source.getNu(), pol_source.getPolFrac(), initial_coeffs=frac_coeff)

    # Create a polarized source using the optimized coefficients values of polarization fraction and degree
    opt_polsource = PolarizedSource(nu=pol_source.getNu(), poldeg=poldeg_func.f_eval(pol_source.getNu(), deg_coeff_popt), polfrac=polfrac_func.f_eval(pol_source.getNu(), frac_coeff_popt))

    # Plot polarization fraction and angle both data and fit
    fig1, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(pol_source.getNu() / 1e9, pol_source.getPolFrac(), 'b-', label='data')
    ax1.plot(opt_polsource.getNu() / 1e9, opt_polsource.getPolFrac(), 'r--', label='fit')
    ax1.set_ylabel('Polarization fraction')
    ax1.grid()
    ax1.set_title(knownSource)
    ax1.legend()

    ax2.plot(pol_source.getNu() / 1e9, pol_source.getPolDeg(), 'g-', label='data')
    ax2.plot(opt_polsource.getNu() / 1e9, opt_polsource.getPolDeg(), 'r--', label='fit')
    ax2.set_ylabel('Polarization angle [deg]')
    ax2.grid()
    ax2.set_xlabel('Frequency [GHz]')
    ax2.legend()


    plt.tight_layout()
    plt.show()
