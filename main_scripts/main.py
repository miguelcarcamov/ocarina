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
    print sys.argv
    inputvis = sys.argv[3]
    refant = sys.argv[4]
    complexgaintable = sys.argv[5]
    bandpasstable = sys.argv[6]
    phasecaltable = sys.argv[7]
    leakagefield = sys.argv[8]  # Leakage
    polanglefield = sys.argv[9]
    target = sys.argv[10]
    polsource = sys.argv[11]

    # Initialize polarization calibration object
    spw_table = queryTable(query="SELECT REF_FREQUENCY FROM "+inputvis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
    spw_ids = spw_table.rownumbers()
    spw_reffreqs = spw_table.getcol("REF_FREQUENCY")[0]
    nu_0 = np.median(spw_reffreqs)
    caltables = [complexgaintable, bandpasstable, phasecaltable]
    polcal = PolCalibration(vis=inputvis, spw_ids=spw_ids, polanglefield=polanglefield, leakagefield=leakagefield, target=target, refant=refant, calibration_tables=caltables)

    # Create your model source
    p_source = PolarizedSource(knownSource=polsource)
    polcal.setModel(pol_source_object = p_source, standard="Perley-Butler 2013", field=polanglefield, epoch="2012", nu_0=nu_0, nterms_angle=5, nterms_frac=4)
    leakagetable = polcal.calibrateLeakage(gaintable=caltables)
    polcal.plotLeakage(plotdir=".")
    caltablesall = caltables.append(leakagetable)
    polangletable = polcal.calibratePolAngle(gaintable=caltablesall)
    caltablesall2 = caltablesall.append(polangletable)
    polcal.applySolutions(gaintables=caltablesall2)
