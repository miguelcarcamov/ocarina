import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from scipy.constants import speed_of_light as c
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '.'))
from polcalibration.polarizedsource import PolarizedSource
from polcalibration.function import PolFunction, FluxFunction
from matplotlib.ticker import ScalarFormatter, NullFormatter
GENERAL_SIZE = 21
plt.rcParams['text.latex.preamble']=[r"\usepackage{mathpazo}"]
plt.rcParams.update({
    #'figure.figsize': (W, W/(4/3)),     # 4:3 aspect ratio
    'font.size' : GENERAL_SIZE,                   # Set font size to 11pt
    'axes.titlesize' : GENERAL_SIZE,
    'axes.labelsize': GENERAL_SIZE,               # -> axis labels
    'legend.fontsize': GENERAL_SIZE,              # -> legends
    'xtick.labelsize': GENERAL_SIZE,    # fontsize of the tick labels
    'ytick.labelsize': GENERAL_SIZE,    # fontsize of the tick labels
    'font.family': 'serif',
    'font.serif': 'Palatino',
    'text.usetex': True
})

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

    # Get flux density coefficients
    p3c286_flux_coeffs = p3c286.getCoeffs(standard="Perley-Butler 2017", epoch="2017")
    p3c147_flux_coeffs = p3c147.getCoeffs(standard="Perley-Butler 2017", epoch="2017")

    p3c286_flux_coeff_errs = p3c286.getCoeffErrs()
    p3c147_flux_coeff_errs = p3c147.getCoeffErrs()

    print("Perley and Butler 3C286")
    print(p3c286_flux_coeffs)
    print(p3c286_flux_coeff_errs)

    print("Perley and Butler 3C147")
    print(p3c147_flux_coeffs)
    print(p3c147_flux_coeff_errs)

    p3c286_upper_coeffs = np.array(p3c286_flux_coeffs) + np.array(p3c286_flux_coeff_errs)
    p3c286_lower_coeffs = np.array(p3c286_flux_coeffs) - np.array(p3c286_flux_coeff_errs)

    p3c147_upper_coeffs = np.array(p3c147_flux_coeffs) + np.array(p3c147_flux_coeff_errs)
    p3c147_lower_coeffs = np.array(p3c147_flux_coeffs) - np.array(p3c147_flux_coeff_errs)

    flux_p3c286 = p3c286.flux(nu/1e9)
    flux_nu_0_p3c286 = p3c286.flux_scalar(nu_0/1e9)
    flux_p3c147 = p3c147.flux(nu/1e9)
    flux_nu_0_p3c147 = p3c147.flux_scalar(nu_0/1e9)

    flux_p3c286_upper = p3c286.flux_giving_coeffs(nu/1e9, p3c286_upper_coeffs.tolist())
    flux_nu_0_p3c286_upper = p3c286.flux_scalar_giving_coeffs(nu_0/1e9, p3c286_upper_coeffs.tolist())

    flux_p3c147_upper = p3c147.flux_giving_coeffs(nu/1e9, p3c147_upper_coeffs.tolist())
    flux_nu_0_p3c147_upper = p3c147.flux_scalar_giving_coeffs(nu_0/1e9, p3c147_upper_coeffs.tolist())

    flux_p3c286_lower = p3c286.flux_giving_coeffs(nu/1e9, p3c286_lower_coeffs.tolist())
    flux_nu_0_p3c286_lower = p3c286.flux_scalar_giving_coeffs(nu_0/1e9, p3c286_lower_coeffs.tolist())

    flux_p3c147_lower = p3c147.flux_giving_coeffs(nu/1e9, p3c147_lower_coeffs.tolist())
    flux_nu_0_p3c147_lower = p3c147.flux_scalar_giving_coeffs(nu_0/1e9, p3c147_lower_coeffs.tolist())

    error_flux_3c286 = 0.5*(np.array(flux_p3c286_upper) - np.array(flux_p3c286_lower))
    error_flux_3c147= 0.5*(np.array(flux_p3c147_upper) - np.array(flux_p3c147_lower))

    p3c286_spidx, p3c286_spidx_errs = p3c286.fitAlphaBeta(nu=nu, nu_0=nu_0)
    p3c147_spidx, p3c147_spidx_errs = p3c147.fitAlphaBeta(nu=nu, nu_0=nu_0)

    print("3C286 Flux")
    print(p3c286_spidx)
    print(p3c286_spidx_errs)
    print("3C147 Flux")
    print(p3c147_spidx)
    print(p3c147_spidx_errs)

    p3c286_spidx_upper = np.array(p3c286_spidx) + np.array(p3c286_spidx_errs)
    p3c147_spidx_upper = np.array(p3c147_spidx) + np.array(p3c147_spidx_errs)

    p3c286_spidx_lower= np.array(p3c286_spidx) - np.array(p3c286_spidx_errs)
    p3c147_spidx_lower = np.array(p3c147_spidx) + np.array(p3c147_spidx_errs)

    fluxFunction_p3c286 = FluxFunction(flux_0 = flux_nu_0_p3c286, x_0 = nu_0)
    fluxFunction_p3c147 = FluxFunction(flux_0 = flux_nu_0_p3c147, x_0 = nu_0)

    fluxFunction_p3c286_upper = FluxFunction(flux_0 = flux_nu_0_p3c286_upper, x_0 = nu_0)
    fluxFunction_p3c147_upper = FluxFunction(flux_0 = flux_nu_0_p3c147_upper, x_0 = nu_0)

    fluxFunction_p3c286_lower = FluxFunction(flux_0 = flux_nu_0_p3c286_lower, x_0 = nu_0)
    fluxFunction_p3c147_lower = FluxFunction(flux_0 = flux_nu_0_p3c147_lower, x_0 = nu_0)

    fitted_flux_p3c286_upper = fluxFunction_p3c286.f_eval(nu, p3c286_spidx_upper.tolist())
    fitted_flux_p3c147_upper = fluxFunction_p3c147.f_eval(nu, p3c147_spidx_upper.tolist())

    fitted_flux_p3c286_lower = fluxFunction_p3c286.f_eval(nu, p3c286_spidx_lower.tolist())
    fitted_flux_p3c147_lower = fluxFunction_p3c147.f_eval(nu, p3c147_spidx_lower.tolist())

    fitted_flux_p3c286 = fluxFunction_p3c286.f_eval(nu, p3c286_spidx)
    fitted_flux_p3c147 = fluxFunction_p3c147.f_eval(nu, p3c147_spidx)

    # Get polarization fraction coeffs
    p3c286_polfrac_coeffs, p3c286_polfrac_coeffs_errs = p3c286.getPolFracCoeffs(nterms=2, nu_min=1.008*1e9, nu_max=2.031*1e9)
    p3c147_polfrac_coeffs, p3c147_polfrac_coeffs_errs = p3c147.getPolFracCoeffs(nterms=2, nu_min=1.008*1e9, nu_max=2.031*1e9)

    # Get polarization angle coeffs
    p3c286_polangle_coeffs, p3c286_polangle_coeffs_errs = p3c286.getPolAngleCoeffs(nterms=2, nu_min=1.008*1e9, nu_max=2.031*1e9)
    p3c147_polangle_coeffs, p3c147_polangle_coeffs_errs = p3c147.getPolAngleCoeffs(nterms=2, nu_min=1.008*1e9, nu_max=2.031*1e9)

    print("3C286 Pol angle")
    print(p3c286_polangle_coeffs)
    print(p3c286_polangle_coeffs_errs)
    print("3C286 Pol fraction")
    print(p3c286_polfrac_coeffs)
    print(p3c286_polfrac_coeffs_errs)

    print("3C147 Pol angle")
    print(p3c147_polangle_coeffs)
    print(p3c147_polangle_coeffs_errs)
    print("3C147 Pol fraction")
    print(p3c147_polfrac_coeffs)
    print(p3c147_polfrac_coeffs_errs)

    fig, axs = plt.subplots(1,1)
    nu_GHz = nu/1e9
    data = flux_p3c286
    #print(data)
    #print(p3c286_polangle_coeffs)
    data_fit = fitted_flux_p3c286
    #axs.plot(nu_GHz, data, "o", label="Data", color='black', ms=5.0)
    axs.plot(nu_GHz, data_fit, label="Fit", color='red', linewidth=2)
    axs.errorbar(nu_GHz, data, yerr=error_flux_3c286, fmt="o", label="Data", color='black', ms=5.0, capsize=7.0, barsabove=True)
    #axs.fill_between(nu_GHz, flux_p3c286_lower, flux_p3c286_upper, color='black', alpha=0.2)
    axs.fill_between(nu_GHz, fitted_flux_p3c286_lower, fitted_flux_p3c286_upper, color='red', alpha=0.2)
    axs.set_xscale("log")
    axs.set_yscale("log")
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
    plt.savefig('3c286_fluxfit.pdf', bbox_inches='tight')
    #posx = np.max(nu_GHz)- 6
    #posy = np.max(data)- 8
    #axs.text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    fig, axs = plt.subplots(1,1)
    nu_GHz = nu/1e9
    data = flux_p3c147
    #print(data)
    #print(p3c286_polangle_coeffs)
    data_fit = fitted_flux_p3c147
    #print(data_fit)
    #axs.plot(nu_GHz, data, "o", label="Data", color='black', ms=5.0)
    axs.errorbar(nu_GHz, data, yerr=error_flux_3c286, fmt="o", label="Data", color='black', ms=5.0, capsize=7.0, barsabove=True)
    axs.plot(nu_GHz, data_fit, label="Fit", color='red', linewidth=2)
    axs.fill_between(nu_GHz, fitted_flux_p3c147_lower, fitted_flux_p3c147_upper, color='red', alpha=0.2)
    axs.set_xscale("log")
    axs.set_yscale("log")
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
    plt.savefig('3c147_fluxfit.pdf', bbox_inches='tight')


    fig, axs = plt.subplots(1,1)
    nu_Hz = p3c286.getNu()
    nu_GHz = p3c286.getNu()/1e9
    lambda2_m = (c/nu_Hz)**2
    lambda2_m = lambda2_m[::-1]
    l_band_idx = np.where((nu_GHz>=1.008) & (nu_GHz<=2.031))
    l_band_l2 = (c/nu_Hz[l_band_idx])**2
    l_band_l2 = l_band_l2[::-1]
    data = p3c286.getPolAngleDegrees()
    data_l2 = data[::-1]

    plot_pol_function_p3c286 = PolFunction(x_0=np.median(nu_GHz[l_band_idx]))
    plot_pol_function_p3c147 = PolFunction(x_0=np.median(nu_GHz[l_band_idx]))

    data_fit = plot_pol_function_p3c286.f_eval(nu_GHz[l_band_idx], p3c286_polangle_coeffs) * 180.0 / np.pi
    data_fit_l2 = data_fit[::-1]
    #print(p3c286_polangle_coeffs)
    #print(data_fit)
    #axs.plot(JVLA_nu, data_fit, label="Fit", color='red', linewidth=2)
    axs.plot(lambda2_m, data_l2, 'o', label="Data", color='black', ms=7.0)
    axs.plot(l_band_l2, data_fit_l2, label="Fit", color='red', linewidth=2.5)
    axs.fill_between(l_band_l2, 32, np.max(data_l2), color = 'cyan', alpha = 0.2)
    #axs.xaxis.set_major_locator(plt.MaxNLocator(8))
    #axs.yaxis.set_major_locator(plt.MaxNLocator(6))
    axs.set_ylabel("Polarization angle (degrees)")
    axs.set_xlabel(r'$\lambda^2$ (m$^2$)')
    axs.grid()
    axs.legend()
    axs.set_xlim(np.min(lambda2_m), np.max(lambda2_m))
    axs.set_ylim(32, np.max(data_l2))
    #posx = np.max(nu_GHz)- 3
    #posy = np.max(data)- 0.025
    #axs.text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    plt.savefig('3c286_polanglefit.pdf', bbox_inches='tight')

    range = 1e-2
    fig, axs = plt.subplots(1,1)
    #print(p3c286_polangle_coeffs)
    axs.plot(l_band_l2, data_fit_l2, label="Fit", color='red', linewidth=2.5)
    axs.plot(lambda2_m, data_l2, 'o', label="Data", color='black', ms=7.0)
    axs.fill_between(l_band_l2, 32, np.max(data_l2), color = 'cyan', alpha = 0.2)
    #axs.set_markersize()
    #axs.xaxis.set_major_locator(plt.MaxNLocator(8))
    #axs.yaxis.set_major_locator(plt.MaxNLocator(6))
    axs.set_ylabel("Polarization angle (degrees)")
    axs.set_xlabel(r'$\lambda^2$ (m$^2$)')
    axs.grid()
    axs.legend()
    axs.set_xlim(np.min(l_band_l2), np.max(l_band_l2))
    axs.set_ylim(np.min(data_fit_l2)-range, np.max(data_fit_l2)+range)
    #posx = np.max(nu_GHz)- 3
    #posy = np.max(data)- 0.025
    #axs.text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    plt.savefig('3c286_polanglefit_2.pdf', bbox_inches='tight')

    data_fit = plot_pol_function_p3c286.f_eval(nu_GHz[l_band_idx], p3c286_polfrac_coeffs)
    data_fit_l2 = data_fit[::-1]
    fig, axs = plt.subplots(1,1)
    data = p3c286.getPolFrac()
    data_l2 = data[::-1]
    #print(p3c286_polangle_coeffs)
    #axs.plot(JVLA_nu, data_fit, label="Fit", color='red', linewidth=2)
    axs.plot(lambda2_m, data_l2, 'o', label="Data", color='black', ms=7.0)
    axs.plot(l_band_l2, data_fit_l2, label="Fit", color='red', linewidth=2.5)
    axs.fill_between(l_band_l2, np.min(data_l2), np.max(data_l2), color = 'cyan', alpha = 0.2)
    #axs.xaxis.set_major_locator(plt.MaxNLocator(8))
    #axs.yaxis.set_major_locator(plt.MaxNLocator(6))
    axs.set_ylabel("Polarization fraction")
    axs.set_xlabel(r'$\lambda^2$ (m$^2$)')
    axs.grid()
    axs.legend()
    axs.set_xlim(np.min(lambda2_m), np.max(lambda2_m))
    axs.set_ylim(np.min(data_l2), np.max(data_l2))
    #posx = np.max(nu_GHz)- 3
    #posy = np.max(data)- 0.025
    #axs.text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    plt.savefig('3c286_polfracfit.pdf', bbox_inches='tight')

    range = 1e-2
    fig, axs = plt.subplots(1,1)
    #print(p3c286_polangle_coeffs)
    #print(data_fit_l2)
    axs.plot(l_band_l2, data_fit_l2, label="Fit", color='red', linewidth=2.5)
    axs.plot(lambda2_m, data_l2, 'o', label="Data", color='black', ms=7.0)
    axs.fill_between(l_band_l2, 0.08, np.max(data_l2), color = 'cyan', alpha = 0.2)
    #axs.xaxis.set_major_locator(plt.MaxNLocator(8))
    #axs.yaxis.set_major_locator(plt.MaxNLocator(6))
    axs.set_ylabel("Polarization fraction")
    axs.set_xlabel(r'$\lambda^2$ (m$^2$)')
    axs.grid()
    axs.legend()
    axs.set_xlim(np.min(l_band_l2), np.max(l_band_l2))
    axs.set_ylim(0.08, np.max(data_fit_l2)+range)
    #posx = np.max(nu_GHz)- 3
    #posy = np.max(data)- 0.025
    #axs.text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    plt.savefig('3c286_polfracfit_2.pdf', bbox_inches='tight')

    os.system('rm -rf *.log')
    os.system('rm -rf *.last')
