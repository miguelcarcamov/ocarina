import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '.'))

from polcalibration.polarizedsource import PolarizedSource
from polcalibration.function import PolFunction
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter, NullFormatter



rc('text', usetex=True)
rc('font', family='serif')


def formatAxes(ax):
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(NullFormatter())


if __name__ == '__main__':
    allnupoints = np.array([0.3275, 1.015, 1.275, 1.465, 1.865, 2.565, 3.565, 4.535, 4.835, 4.885,
                   6.135, 6.885, 7.465, 8.435, 8.485, 8.735, 11.06, 12.890, 14.635, 14.715,
                   14.915, 14.965, 17.422, 18.230, 18.485, 18.585, 20.485, 22.460, 22.835, 24.450,
                   25.836, 26.485, 28.450, 29.735, 36.435, 43.065, 43.340, 48.350, 48.565])
    nu_points_3c123 = np.array([0.3275, 1.015, 1.275, 1.465, 1.865, 2.565, 3.565, 4.535, 4.835, 4.885,
                       6.135, 6.885, 7.465, 8.435, 8.485, 8.735, 11.06, 14.635, 14.715,
                       14.915, 14.965, 17.422, 18.485, 18.585, 20.485, 22.460, 22.835, 24.450,
                       25.836, 26.485, 28.450, 29.735, 36.435, 43.340, 48.350])

    datapoints_3c123 = np.array([145.0, 66.2, 46.6, 47.8, 38.7, 28.9, 21.4, 16.9, 16.0, 15.88,
                        12.81, 11.20, 11.01, 9.20, 9.10, 8.86, 6.73, 5.34, 5.02, 5.132,
                        5.092, 4.272, 4.090, 3.934, 3.586, 3.297, 3.334, 2.867, 2.697, 2.716,
                        2.436, 2.453, 1.841, 1.421, 1.269])

    sigmapoints_3c123 = np.array([4.3, 4.3, 3.2, 0.5, 0.6, 0.3, 0.8, 0.2, 0.2, 0.1, 0.15, 0.14, 0.2,
                                   0.04, 0.15, 0.05, 0.15, 0.05, 0.05, 0.025, 0.028, 0.07, 0.055, 0.055, 0.055,
                                   0.022, 0.06, 0.03, 0.06, 0.05, 0.06, 0.05, 0.17, 0.055, 0.12])

    datapoints_3c196 = np.array([46.8, 20.1, 13.3, 14.1, 11.3, 8.16, 6.22, 4.55, 4.22, 4.189, 3.318,
                   2.85, 2.79, 2.294, 2.275, 2.202, 1.64, 1.388, 1.255, 1.206, 1.207, 1.198, 0.988,
                   0.932, 0.947, 0.926, 0.820, 0.745, 0.760, 0.657, 0.620, 0.607, 0.568, 0.529, 0.408,
                   0.367, 0.342, 0.289, 0.272])

    sigmapoints_3c196 = np.array([1.4, 4.8, 2.0, 0.2, 0.2, 0.1, 0.2, 0.06, 0.1, 0.025, 0.05, 0.05, 0.05, 0.010,
                                  0.03, 0.011, 0.03, 0.025, 0.020, 0.020, 0.004, 0.007, 0.02, 0.020, 0.015, 0.015,
                                  0.010, 0.003, 0.010, 0.017, 0.017, 0.017, 0.015, 0.015, 0.005, 0.015, 0.005, 0.005, 0.015])
    datapoints_3c286 = np.array([26.1, 18.4, 13.8, 15.0, 13.2, 10.9, 9.5, 7.68, 7.33, 7.297, 6.49, 5.75,
                                5.70, 5.059, 5.045, 4.930, 4.053, 3.662, 3.509, 3.375, 3.399, 3.387, 2.980,
                                2.860, 2.925, 2.880, 2.731, 2.505, 2.562, 2.387, 2.181, 2.247, 2.079, 2.011,
                                1.684, 1.658, 1.543, 1.449, 1.465])

    sigmapoints_3c286 = np.array([0.8, 4.3, 2.0, 0.2, 0.2, 0.2, 0.1, 0.1, 0.2, 0.046, 0.15, 0.05, 0.10, 0.021,
                                  0.07, 0.024, 0.08, 0.07, 0.04, 0.04, 0.016, 0.015, 0.04, 0.045, 0.045, 0.04,
                                  0.05, 0.016, 0.05, 0.03, 0.06, 0.05, 0.05, 0.05, 0.02, 0.08, 0.024, 0.04, 0.1])
    datapoints_3c295 = np.array([60.8, 30.8, 21.5, 22.2, 17.9, 12.8, 9.62, 6.96, 6.45, 6.37, 4.99, 4.21, 4.13,
                                3.319, 3.295, 3.173, 2.204, 1.904, 1.694, 1.630, 1.626, 1.617, 1.311, 1.222, 1.256,
                                1.221, 1.089, 0.952, 0.967, 0.861, 0.770, 0.779, 0.689, 0.653, 0.484, 0.442, 0.398, 0.359,
                                0.325])

    sigmapoints_3c295 = np.array([1.8, 7.3, 3.0, 0.5, 0.3, 0.2, 0.2, 0.09, 0.15, 0.04, 0.05, 0.05, 0.07, 0.014, 0.05, 0.016,
                                   0.05, 0.04, 0.04, 0.03, 0.008, 0.007, 0.025, 0.05, 0.02, 0.015, 0.015, 0.005, 0.015, 0.02,
                                   0.02, 0.020, 0.020, 0.020, 0.015, 0.020, 0.006, 0.013, 0.025])
    print("len nu: ", len(nu_points_3c123))
    print("3c295: ", len(sigmapoints_3c123))
    nu = np.linspace(0.3275 * 1e9, 50.0 * 1e9, 40)
    p3c123 = PolarizedSource(source="3C123")
    p3c196 = PolarizedSource(source="3C196")
    p3c286 = PolarizedSource(source="3C286")
    p3c295 = PolarizedSource(source="3C295")

    p3c123.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    p3c196.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    p3c286.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    p3c295.getCoeffs(standard="Perley-Butler 2013", epoch="2012")

    fig, axs = plt.subplots(2, 2)

    markersize = 10
    fmt = 'k.'
    nu_GHz = nu / 1e9
    flux0 = p3c123.flux(nu_GHz)
    axs[0, 0].loglog(nu_GHz, flux0, 'r-')
    #axs[0, 0].loglog(nu_points_3c123, datapoints_3c123, 'k.')
    axs[0, 0].errorbar(nu_points_3c123, datapoints_3c123, yerr=sigmapoints_3c123, fmt=fmt, markersize=markersize)
    axs[0, 0].set_xlim([np.min(nu_GHz), np.max(nu_GHz)])
    axs[0, 0].set_ylim([np.min(flux0), np.max(flux0)])
    formatAxes(axs[0, 0])
    axs[0, 0].set_xticks(nu_GHz)
    axs[0, 0].set_yticks(flux0)
    axs[0, 0].xaxis.set_major_locator(plt.MaxNLocator(10))
    axs[0, 0].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0, 0].grid()
    axs[0, 0].set_ylabel("Flux Density (Jy)")
    posx = np.max(nu_GHz) - 10
    posy = np.max(flux0) - 8
    axs[0, 0].text(posx, posy, p3c123.getName(), style='italic', bbox={
                   'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    flux1 = p3c196.flux(nu_GHz)
    axs[0, 1].loglog(nu_GHz, flux1, 'r-')
    #axs[0, 1].loglog(allnupoints, datapoints_3c196, 'k.')
    axs[0, 1].errorbar(allnupoints, datapoints_3c196, yerr=sigmapoints_3c196, fmt=fmt, markersize=markersize)
    axs[0, 1].set_xlim([np.min(nu_GHz), np.max(nu_GHz)])
    axs[0, 1].set_ylim([np.min(flux1), np.max(flux1)])
    formatAxes(axs[0, 1])
    axs[0, 1].set_xticks(nu_GHz)
    axs[0, 1].set_yticks(flux1)
    axs[0, 1].xaxis.set_major_locator(plt.MaxNLocator(10))
    axs[0, 1].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0, 1].grid()
    posx = np.max(nu_GHz) - 12
    posy = np.max(flux1) - 4
    axs[0, 1].text(posx, posy, p3c196.getName(), style='italic', bbox={
                   'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    flux2 = p3c286.flux(nu_GHz)
    axs[1, 0].loglog(nu_GHz, flux2, 'r-')
    axs[1, 0].errorbar(allnupoints, datapoints_3c286, yerr=sigmapoints_3c286, fmt=fmt, markersize=markersize)
    axs[1, 0].set_xlim([np.min(nu_GHz), np.max(nu_GHz)])
    axs[1, 0].set_ylim([np.min(flux2), np.max(flux2)])
    formatAxes(axs[1, 0])
    axs[1, 0].set_xticks(nu_GHz)
    axs[1, 0].set_yticks(flux2)
    axs[1, 0].xaxis.set_major_locator(plt.MaxNLocator(10))
    axs[1, 0].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[1, 0].grid()
    axs[1, 0].set_xlabel("Frequency (GHz)")
    axs[1, 0].set_ylabel("Flux Density (Jy)")
    posx = np.max(nu_GHz) - 12
    posy = np.max(flux2) - 3
    axs[1, 0].text(posx, posy, p3c286.getName(), style='italic', bbox={
                   'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    flux3 = p3c295.flux(nu_GHz)
    axs[1, 1].loglog(nu_GHz, flux3, 'r-')
    axs[1, 1].errorbar(allnupoints, datapoints_3c295, yerr=sigmapoints_3c295, fmt=fmt, markersize=markersize)
    axs[1, 1].set_xlim([np.min(nu_GHz), np.max(nu_GHz)])
    axs[1, 1].set_ylim([np.min(flux3), np.max(flux3)])
    formatAxes(axs[1, 1])
    axs[1, 1].set_xticks(nu_GHz)
    axs[1, 1].set_yticks(flux3)
    axs[1, 1].xaxis.set_major_locator(plt.MaxNLocator(10))
    axs[1, 1].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[1, 1].grid()
    axs[1, 1].set_xlabel("Frequency (GHz)")
    posx = np.max(nu_GHz) - 12
    posy = np.max(flux3) - 10
    axs[1, 1].text(posx, posy, p3c295.getName(), style='italic', bbox={
                   'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    plt.tight_layout()
    plt.savefig('polsources.png')
    print("Done")
