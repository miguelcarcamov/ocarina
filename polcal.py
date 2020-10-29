from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import numpy as np
import matplotlib
import os
import sys
from __casac__.table import table as tb
from applycal import applycal
from plotcal import plotcal
from polcal import polcal
from applycal import applycal
from gaincal import gaincal

matplotlib.rcParams['text.usetex'] = True

def degreesToRadians(degrees):
    return degrees * np.pi / 180.0

def radiansToDegrees(radians):
    return radians * 180.0 / np.pi

def queryTable(query=""):
    return tb.taql(query+" AS "+tablename)

class PolFunction(object):
    def __init__(self, nu_0=0.0, terms=0, **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.coeff = []

    def getCoeffs(self):
        return self.coeff

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

    def fit(self, nu, data, initial_coeffs=None, lowerbound=-np.inf, upperbound=np.inf):
        lowerbounds = np.ones(len(initial_coeffs))*lowerbound
        upperbounds = np.ones(len(initial_coeffs))*upperbound
        popt, pcov = curve_fit(self.f, nu, data, p0=initial_coeffs, check_finite=True, bounds=(lowerbounds,upperbounds))
        self.coeff = popt
        return popt, pcov

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

    def p3c48(self):
        self.nu = np.array([1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80, 12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5])
        self.polangle = np.array([25.0, 140.0, -5.0, -150.0, -120.0, -100.0, -92.0, -84.0, -75.0, -72.0, -68.0, -67.0, -64.0, -62.0, -62.0, -62.0, -63.0, -64.0, -66.0, -67.0, -70.0, -70.0, -77.0, -85.0])
        self.polfrac = np.array([0.3, 0.5, 0.7, 0.9, 1.4, 2.0, 2.5, 3.2, 3.8, 4.2, 5.2, 5.2, 5.3, 5.4, 6.0, 6.1, 6.4, 6.4, 6.9, 7.1, 7.7, 7.8, 7.4, 7.5])
        self.name = "3C48"

    def p3c48_2019(self):
        self.nu = np.array(np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3, 14.1, 16.6, 19.1, 25.6, 32.1, 37.1, 42.1, 48.1]))
        self.polangle = np.array([4.3, -34, 23, 67.1, -84.0, -72, -66, -63, -62, -63, -64, -68, -72, -76, -77, -84, -84])
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
    def filter(self, nu=np.inf):
        indexes = np.where(self.nu < nu)
        self.nu = self.nu[indexes]
        self.polfrac = self.polfrac[indexes]
        self.polangle = self.polangle[indexes]

    def querySpecIdx(self, standard="Perley-Butler 2013", epoch="2012"):
        coeff_table = os.getenv('CASAPATH').split(' ')[0] + '/data/nrao/VLA/standards/' + self.spix_dict[standard]
        print coeff_table
        tb_obj = tb()
        tb_obj.open(coeff_table)
        query_table = tb_obj.taql("select * from "+coeff_table+" where Epoch="+epoch)
        coeffs = query_table.getcol(self.name+"_coeffs")
        # get alpha
        alpha = coeffs[1][0]/coeffs[0][0]
        # get beta
        beta = coeffs[2][0]/coeffs[0][0] - alpha*(alpha-1)/2
        tb_obj.close()
        return [alpha, beta]

    def getPolFracCoeffs(self, nu_0=0.0, nterms=3):
        if(not nu_0):
            nu_0 = np.median(self.nu)
        ifrac_coeffs = np.random.uniform(0.0, 1.0, nterms)
        source_func_frac = PolFunction(nu_0=nu_0, nterms=nterms)
        source_func_frac.fit(self.getNu(), self.getPolFrac(), ifrac_coeffs, 0.0, 1.0)
        return source_func_frac.getCoeffs()

    # Returns pol angle coeffs in radians
    def getPolAngleCoeffs(self, nu_0=0.0, nterms=3):
        if(not nu_0):
            nu_0 = np.median(self.nu)
        iangle_coeffs = np.random.uniform(-np.pi, np.pi, nterms)
        source_func_angle = PolFunction(nu_0=nu_0, nterms=nterms)
        source_func_angle.fit(self.getNu(), self.getPolAngle(), iangle_coeffs, -np.pi, np.pi)
        return source_func_angle.getCoeffs()

class PolCalibration(object):
    def __init__(self, vis="", nspw=0, spw_ids=np.array([]), spw_reffreq=np.array([]), bcal="", lcal="", pcal="", refant="", calibration_tables=[], **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)


    def setModel(self, pol_source_object = None, standard="Perley-Butler 2013", field="", epoch="2012", nu_0=0.0, nterms_angle=3, nterms_frac=3):

        # get spectral idx coeffs from VLA tables
        spec_idx = pol_source_object.querySpecIdx(standard=standard, epoch=epoch)
        pol_frac_coeffs = pol_source_object.getPolFracCoeffs(nu_0=nu_0, nterms=nterms_frac)
        pol_angle_coeffs = pol_source_object.getPolAngleCoeffs(nu_0=nu_0, nterms=nterms_angle)
        # get intensities in frequency
        calflux = setjy(vis=self.vis, standard=standard, field=field, spw="*")
        fid = calflux.keys()[0]

        # loop spws
        for i in range(0, self.nspw):
            spw_id = self.spw_ids[i]
            intensity = calflux[fid][spw_id]['fluxd'][0]
            reffreq_GHz = self.spw_reffreq[i]/1e9
            setjy(vis=self.vis, field=field, standard='manual', spw=str(spw_id), fluxdensity=[intensity,0,0,0], spix=spec_idx, reffreq=str(reffreq_GHz)+"GHz", polindex=pol_frac_coeffs, polangle=pol_angle_coeffs, scalebychan=True, usescratch=False)

    def calibrateLeakage(self, minsnr=3, poltype="D", caltable="D0", gaintable=[], gainfield=[]):
        if os.path.exists(caltable): rmtables(caltable)
        polcal(vis=self.vis, caltable=caltable, field=self.lcal, spw='', refant=self.refant, poltype=poltype, solint='inf', combine='scan', minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)
        return caltable

    def calibratePolAngle(self, minsnr=3, caltable="X0", poltype="Xf", gaintable=[], gainfield=[]):
        if os.path.exists(caltable): rmtables(caltable)
        polcal(vis=self.vis, caltable=caltable, field=self.bcal, spw='', refant=self.refant, poltype=poltype, solint='inf', combine='scan', minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)
        return caltable

    def plotLeakage(self, leakagecaltable="DO", plotdir=""):
        plotcal(caltable=leakagecaltable, xaxis='antenna', yaxis='amp', figfile=plotdir+self.vis[:-3]+'.D0.amp.png', showgui=False)
        plotcal(caltable=leakagecaltable, xaxis='antenna', yaxis='phase', iteration='antenna', figfile=plotdir+self.vis[:-3]+'.D0.phs.png', showgui=False)
        plotcal(caltable=leakagecaltable, xaxis='antenna', yaxis='snr', showgui=False, figfile=plotdir+self.vis[:-3]+'.D0.snr.png')
        plotcal(caltable=leakagecaltable, xaxis='real', yaxis='imag', showgui=False, figfile=plotdir+self.vis[:-3]+'.D0.cmplx.png')

    def applySolutions(self, gaintables=[], gainfields=[]):
        spw = np.array2string(self.spw_ids, separator=',').replace("[", "").replace("]","")
        applycal(vis=self.vis, field=self.bcal, spw=spw, gaintable=[self.curvetable, self.gaincaltable,'D0','X0'], gainfield=['',self.bcal,'',''], calwt=[False], parang=True)
        applycal(vis=self.vis, field=self.lcal, spw=spw, gaintable=[self.curvetable, self.gaincaltable, self.fluxtable,'D0','X0'], gainfield=['',self.lcal,self.lcal,'',''], calwt=[False], parang=True)
        applycal(vis=self.vis, field=self.pcal, spw=spw, gaintable=[self.curvetable, self.gaincaltable, self.fluxtable,'D0','X0'], gainfield=['',self.pcal,self.pcal,'',''], calwt=[False], parang=True)

    #def calibrate(self, known_pol_object=None, standard="Perley-Butler 2013", epoch="2012", plotdir="", nu_0=0.0, nterms_angle=3, nterms_frac=3):
    #    self.setModel(pol_source_object=known_pol_object, standard=standard, field=self.bcal, epoch=epoch, nu_0=nu_0, nterms_angle=nterms_angle, nterms_frac=nterms_frac)
    #    self.setModel(pol_source_object=known_pol_object, standard=standard, field=self.lcal, epoch=epoch, nu_0=nu_0, nterms_angle=nterms_angle, nterms_frac=nterms_frac)
    #    leakage_table = self.calibrateLeakage()
    #    self.plotLeakage(plotdir=plotdir)
    #    polangle_table = self.calibratePolAngle()
    #    self.applySolutions()

    def generate_gain_curve(self, standard='Perley-Butler 2013'):
        #set the flux scale:
        setjy(vis=self.vis, standard=standard, field=self.bcal, usescratch=True)
        if os.path.exists('GC'): rmtables('GC')
        # generate the gain curve calibration table:
        gencal(vis=self.vis, caltype='gc', caltable='GC')

    def complex_gain_calibration(self, minsnr=3, plotdir=""):
        if os.path.exists('G0'): rmtables('G0')
        # do the complex gain calibration:
        callist = ','.join([self.bcal,self.lcal,self.pcal])
        gaincal(vis=self.vis, caltable='G0', gaintype='G', combine='', calmode='ap', field=callist, solint='inf', refant=self.refant, minsnr=minsnr, gaintable=['GC'], parang=False)

        plotcal(caltable='G0',xaxis='time', yaxis='amp', showgui=False, figfile=plotdir+self.vis[:-3]+'.G0.amp.png')
        plotcal(caltable='G0',xaxis='time', yaxis='phase', showgui=False, figfile=plotdir+self.vis[:-3]+'.G0.phs.png')

    def transfer_flux_scale(self):
        # transfer flux scale to phase cal:
        if os.path.exists('F0'): rmtables('F0')
        callist = ','.join([self.bcal,self.lcal,self.pcal])
        myFluxes = fluxscale(vis=self.vis, caltable='G0', reference=self.bcal, transfer=callist, fluxtable='F0', incremental=True, display=False, append=False)

if __name__ == '__main__':
    nu_0 = 1.5 * 1e9

    #print(sys.argv)
    #inputvis = sys.argv[1]
    #refant = sys.argv[2]
    #bcal     = sys.argv[3]
    #lcal     = sys.argv[4]  # Leakage
    #pcal    = sys.argv[5]
    #output  = sys.argv[6]
    p3c286 = PolarizedSource(knownSource="3c286_2019")
    spc_idx = p3c286.querySpecIdx(standard="Perley-Butler 2013", epoch="2012")
    p3c286.filter(nu=2.0*1e9)

    nterms_frac = 2
    nterms_angle = 2
    print(spc_idx)
    print(p3c286.getPolFracCoeffs(nu_0=nu_0, nterms=nterms_frac))
    print(p3c286.getPolAngleCoeffs(nu_0=nu_0, nterms=nterms_angle))


    #spw_table = queryTable(query="SELECT REF_FREQUENCY FROM "+inputvis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
    #spw_ids = spw_table.rownumbers()
    #nspw = len(spw_ids)
    #spw_reffreq = spw_table.getcol("REF_FREQUENCY")[0]
    #polcal = PolCalibration(vis=inputvis, nspw=nspw, spw_ids=spw_ids, spw_reffreq=spw_reffreq, bcal=bcal, lcal=lcal, pcal=pcal, refant=refant)
