import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '.'))
from polcalibration.polarizedsource import PolarizedSource
from polcalibration.function import PolFunction, FluxFunction
from matplotlib.ticker import ScalarFormatter, NullFormatter
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

def formatAxes(ax):
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(NullFormatter())

if __name__ == '__main__':
    nu = np.linspace(0.3275*1e9, 50.0*1e9, 40)
    nu_0 = np.median(nu)
    nu_min = np.min(nu)
    nu_max = np.max(nu)

    p3c286 = PolarizedSource(source="3c286")
    p3c147 = PolarizedSource(source="3c147")

    plot_pol_function_p3c286 = PolFunction(x_0=np.median(p3c286.getNu()))
    plot_pol_function_p3c147 = PolFunction(x_0=np.median(p3c147.getNu()))

    # Get flux density coefficients
    p3c286_flux_coeffs = p3c286.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    p3c147_flux_coeffs = p3c147.getCoeffs(standard="Perley-Butler 2013", epoch="2012")

    flux_p3c286 = p3c286.flux(nu/1e9)
    flux_nu_0_p3c286 = p3c286.flux_scalar(nu_0/1e9)
    flux_p3c147 = p3c147.flux(nu/1e9)
    flux_nu_0_p3c147 = p3c147.flux_scalar(nu_0/1e9)

    p3c286_spidx = p3c286.fitAlphaBeta(nu=nu, nu_0=nu_0)
    p3c147_spidx = p3c147.fitAlphaBeta(nu=nu, nu_0=nu_0)

    fluxFunction_p3c286 = FluxFunction(flux_0 = flux_nu_0_p3c286, x_0 = nu_0)
    fluxFunction_p3c147 = FluxFunction(flux_0 = flux_nu_0_p3c147, x_0 = nu_0)

    fitted_flux_p3c286 = fluxFunction_p3c286.f_eval(nu, p3c286_spidx)
    fitted_flux_p3c147 = fluxFunction_p3c147.f_eval(nu, p3c147_spidx)

    # Get polarization fraction coeffs
    p3c286_polfrac_coeffs = p3c286.getPolFracCoeffs(nu_0=np.median(p3c286.getNu()), nterms=6, nu_min=np.min(p3c286.getNu()), nu_max=np.max(p3c286.getNu()))
    p3c147_polfrac_coeffs = p3c147.getPolFracCoeffs(nu_0=np.median(p3c147.getNu()), nterms=6, nu_min=np.min(p3c147.getNu()), nu_max=np.max(p3c147.getNu()))

    # Get polarization angle coeffs
    p3c286_polangle_coeffs = p3c286.getPolAngleCoeffs(nu_0=np.median(p3c286.getNu()), nterms=6, nu_min=np.min(p3c286.getNu()), nu_max=np.max(p3c286.getNu()))
    p3c147_polangle_coeffs = p3c147.getPolAngleCoeffs(nu_0=np.median(p3c147.getNu()), nterms=6, nu_min=np.min(p3c147.getNu()), nu_max=np.max(p3c147.getNu()))

    fig, axs = plt.subplots(1,1)
    nu_GHz = nu/1e9
    data = flux_p3c286
    #print(data)
    #print(p3c286_polangle_coeffs)
    data_fit = fitted_flux_p3c286
    #print(data_fit)
    axs.loglog(nu_GHz, data, label="Data")
    axs.loglog(nu_GHz, data_fit, label="Fit")
    axs.xaxis.set_major_locator(plt.MaxNLocator(8))
    axs.yaxis.set_major_locator(plt.MaxNLocator(6))
    axs.set_xticks(nu_GHz)
    axs.set_yticks(data)
    axs.set_xlim([np.min(nu_GHz), np.max(nu_GHz)])
    axs.set_ylim([np.min(data), np.max(data)])
    formatAxes(axs)
    axs.xaxis.set_major_locator(plt.MaxNLocator(10))
    axs.yaxis.set_major_locator(plt.MaxNLocator(6))
    axs.set_ylabel("Flux Density (Jy)")
    axs.set_xlabel("Frequency (GHz)")
    axs.grid()
    axs.legend()
    posx = np.max(nu_GHz)- 6
    posy = np.max(data)- 8
    axs.text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    """
    fig, axs = plt.subplots(1,1)

    nu_GHz = p3c286.getNu()/1e9
    data = p3c286.getPolAngleDegrees()
    #print(p3c286_polangle_coeffs)
    data_fit = plot_pol_function_p3c286.f_eval(p3c286.getNu(), p3c286_polangle_coeffs) * 180.0 / np.pi
    print(data_fit)
    axs.plot(nu_GHz, data)
    axs.plot(nu_GHz, data_fit)
    axs.xaxis.set_major_locator(plt.MaxNLocator(8))
    axs.yaxis.set_major_locator(plt.MaxNLocator(6))
    axs.grid()
    posx = np.max(nu_GHz)- 3
    posy = np.max(data)- 0.025
    axs.text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    axs[0,0].plot(nu_GHz, data )
    axs[0,0].set_xlim([np.min(nu_GHz),np.max(nu_GHz)])
    axs[0,0].set_ylim([np.min(data),np.max(data)])
    axs[0,0].xaxis.set_major_locator(plt.MaxNLocator(8))
    axs[0,0].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0,0].grid()
    posx = np.max(nu_GHz)- 3
    posy = np.max(data)- 0.025
    axs[0,0].text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    nu_GHz = p3c147.getNu()/1e9
    data = p3c147.getPolAngleDegrees()
    axs[0,1].plot(nu_GHz, data)
    axs[0,1].set_xlim([np.min(nu_GHz),np.max(nu_GHz)])
    axs[0,1].set_ylim([np.min(data),np.max(data)])
    axs[0,1].xaxis.set_major_locator(plt.MaxNLocator(8))
    axs[0,1].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0,1].grid()
    axs[0,1].set_xlabel("Frequency (GHz)")
    posx = np.max(nu_GHz)- 3
    posy = np.max(data)- 0.1
    axs[0,1].text(posx, posy, p3c147.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    """

    plt.tight_layout()
    plt.savefig('fluxes_fit.png')
    os.system('rm -rf *.log')
    os.system('rm -rf *.last')
