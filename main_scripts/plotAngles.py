import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '.'))
from polcalibration.polarizedsource import PolarizedSource
from polcalibration.function import PolFunction
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

    p3c48 = PolarizedSource(source="3c48")
    p3c286 = PolarizedSource(source="3c286")
    p3c138 = PolarizedSource(source="3c138")
    p3c147 = PolarizedSource(source="3c147")

    p3c48.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    p3c286.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    p3c138.getCoeffs(standard="Perley-Butler 2013", epoch="2012")
    p3c147.getCoeffs(standard="Perley-Butler 2013", epoch="2012")

    fig, axs = plt.subplots(2,2)

    nu_GHz = p3c48.getNu()/1e9
    data = p3c48.getPolAngleDegrees()
    
    axs[0,0].plot(nu_GHz, data)
    axs[0,0].set_xlim([np.min(nu_GHz),np.max(nu_GHz)])
    axs[0,0].set_ylim([np.min(data),np.max(data)])
    axs[0,0].xaxis.set_major_locator(plt.MaxNLocator(8))
    axs[0,0].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0,0].grid()
    axs[0,0].set_ylabel("Polarization angle (degrees)")
    posx = np.max(nu_GHz)- 3
    posy = np.max(data)- 0.1
    axs[0,0].text(posx, posy, p3c48.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    nu_GHz = p3c286.getNu()/1e9
    data = p3c286.getPolAngleDegrees()
    axs[0,1].plot(nu_GHz, data )
    axs[0,1].set_xlim([np.min(nu_GHz),np.max(nu_GHz)])
    axs[0,1].set_ylim([np.min(data),np.max(data)])
    axs[0,1].xaxis.set_major_locator(plt.MaxNLocator(8))
    axs[0,1].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[0,1].grid()
    posx = np.max(nu_GHz)- 3
    posy = np.max(data)- 0.025
    axs[0,1].text(posx, posy, p3c286.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    nu_GHz = p3c138.getNu()/1e9
    data = p3c138.getPolAngleDegrees()
    axs[1,0].plot(nu_GHz, data)
    axs[1,0].set_xlim([np.min(nu_GHz),np.max(nu_GHz)])
    axs[1,0].set_ylim([np.min(data),np.max(data)])
    axs[1,0].xaxis.set_major_locator(plt.MaxNLocator(8))
    axs[1,0].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[1,0].grid()
    axs[1,0].set_xlabel("Frequency (GHz)")
    axs[1,0].set_ylabel("Polarization angle (degrees)")
    posx = np.max(nu_GHz)- 3
    posy = np.max(data)- 0.025
    axs[1,0].text(posx, posy, p3c138.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    nu_GHz = p3c147.getNu()/1e9
    data = p3c147.getPolAngleDegrees()
    axs[1,1].plot(nu_GHz, data)
    axs[1,1].set_xlim([np.min(nu_GHz),np.max(nu_GHz)])
    axs[1,1].set_ylim([np.min(data),np.max(data)])
    axs[1,1].xaxis.set_major_locator(plt.MaxNLocator(8))
    axs[1,1].yaxis.set_major_locator(plt.MaxNLocator(6))
    axs[1,1].grid()
    axs[1,1].set_xlabel("Frequency (GHz)")
    posx = np.max(nu_GHz)- 3
    posy = np.max(data)- 0.1
    axs[1,1].text(posx, posy, p3c147.getName(), style='italic', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})


    plt.tight_layout()
    plt.savefig('polangles.png')
    print("Done")
