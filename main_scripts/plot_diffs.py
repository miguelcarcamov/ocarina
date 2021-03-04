import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '.'))
from polarizedsource import PolarizedSource
from function import PolFunction
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

if __name__ == '__main__':
    nu_0 = 1.5 * 1e9

    #print(sys.argv)
    #inputvis = sys.argv[1]
    #refant = sys.argv[2]
    #bcal     = sys.argv[3]
    #lcal     = sys.argv[4]  # Leakage
    #pcal    = sys.argv[5]
    #output  = sys.argv[6]

    p3c286 = PolarizedSource(knownSource="3c48")
    p3c286.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    #flux_0 = p3c286.flux_scalar(nu_0/1e9)

    print(p3c286.flux_scalar(3.0))


    nu = np.linspace(0.3275*1e9, 50.0*1e9, 40)

    alpha = p3c286.fitAlphaBeta(nu, nu_0=3.0*1e9)
    flux = p3c286.flux_scalar(3.0)*(nu/(3.0*1e9))**(alpha[0]+alpha[1]*np.log(nu/(3.0*1e9)))
    flux_VLA = p3c286.flux_scalar(3.0)*(nu/(3.0*1e9))**(-0.90366565-0.14262821*np.log(nu/(3.0*1e9)))

    fig, axs = plt.subplots(3)
    axs[0].plot(nu/1e9, p3c286.flux(nu/1e9), label='Data')
    axs[0].set_title('Flux Jy')
    axs[0].plot(nu/1e9, flux, label='Fit')
    #axs[0].plot(nu/1e9, flux_VLA, label='Fit VLA Example')
    axs[1].plot(p3c286.getNu()/1e9, p3c286.getPolFrac(), label='Data')
    axs[1].set_title('Polarization fraction')
    axs[2].plot(p3c286.getNu()/1e9, p3c286.getPolAngle(), label='Data')
    axs[2].set_title('Polarization Angle (radians)')

    nterms_frac = 4
    nterms_angle = 5

    pol_frac = PolFunction(x_0=3.0*1e9, nterms=nterms_frac)
    pol_angle = PolFunction(x_0=3.0*1e9, nterms=nterms_angle)

    frac_coeffs = p3c286.getPolFracCoeffs(nu_0=3.0*1e9, nterms=nterms_frac)
    angle_coeffs = p3c286.getPolAngleCoeffs(nu_0=3.0*1e9, nterms=nterms_angle, nu_min=2.5*1e9)
    angle_coeffs_VLA = [1.4215,1.36672,-2.12678,3.48384,-2.71914]
    frac_coeffs_VLA = [0.021429,0.0391826,0.00234878,-0.0230125]
    print(frac_coeffs)
    print(angle_coeffs)
    axs[1].plot(p3c286.getNu()/1e9, pol_frac.f_eval(p3c286.getNu(), frac_coeffs), label='Fit')
    #axs[1].plot(p3c286.getNu()/1e9, pol_frac.f_eval(p3c286.getNu(), frac_coeffs_VLA), label='Fit VLA Example')
    axs[2].plot(p3c286.getNu()/1e9, pol_angle.f_eval(p3c286.getNu(), angle_coeffs), label='Fit')
    #axs[2].plot(p3c286.getNu()/1e9, pol_angle.f_eval(p3c286.getNu(), angle_coeffs_VLA), label='Fit VLA Example')
    print(alpha)
    axs[0].legend(loc='upper right')
    axs[1].legend(loc='upper right')
    axs[2].legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('test.png')
    #spw_table = queryTable(query="SELECT REF_FREQUENCY FROM "+inputvis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
    #spw_ids = spw_table.rownumbers()
    #nspw = len(spw_ids)
    #spw_reffreq = spw_table.getcol("REF_FREQUENCY")[0]
    #polcal = PolCalibration(vis=inputvis, nspw=nspw, spw_ids=spw_ids, spw_reffreq=spw_reffreq, bcal=bcal, lcal=lcal, pcal=pcal, refant=refant)
