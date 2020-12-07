import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from polcalibration.polcalibration import PolCalibration
from polcalibration.utils import *
from polcalibration.polarizedsource import *
from polcalibration.function import *
from flagdata import flagdata
from mstransform import mstransform
from statwt import statwt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

if __name__ == '__main__':
    print sys.argv
    inputvis = sys.argv[3]
    krefant = sys.argv[4]
    refant = sys.argv[5]
    complexgaintable = sys.argv[6]
    leakagefield = sys.argv[7]  # Leakage
    polanglefield = sys.argv[8]
    target = sys.argv[9]
    polsource = sys.argv[10]

    fields=[leakagefield, polanglefield, target]
    fields_str = ",".join(fields)

    # Plot before flagging
    plotms(vis=inputvis, field=fields_str, xaxis='frequency',yaxis='amplitude', ydatacolumn='corrected', correlation='RR,LL', showgui=False, plotfile='beforeflagging_parhands.png', overwrite=True)
    # for cross-hands
    plotms(vis=inputvis, field=fields_str, xaxis='frequency',yaxis='amplitude', ydatacolumn='corrected', correlation='RL,LR', showgui=False, plotfile='beforeflagging_crosshands.png', overwrite=True)
    # Check percentage of flagged data!!!


    mode='rflag'
    flagdata(vis=inputvis, field=fields_str, spw='', mode='tfcrop', correlation='ABS_LL,RR', datacolumn='corrected', freqfit='line', extendflags=False, flagbackup=False)
    flagdata(vis=inputvis, field=fields_str, spw='', mode='rflag', correlation='RL,LR', datacolumn='corrected', extendflags=True, flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="0", correlation='LL,RR', mode=mode, action="apply", flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="0", mode='extend', extendpols=True, action="apply", flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="1~2", correlation='LL,RR', mode=mode, action="apply", flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="1~2", mode='extend', extendpols=True, action="apply", flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="3~5", correlation='LL,RR', mode=mode, action="apply", flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="3~5", mode='extend', extendpols=True, action="apply", flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="5~7", correlation='LL,RR', mode=mode, action="apply", flagbackup=False)
    #flagdata(vis=inputvis, field=f, datacolumn="corrected", spw="5~7", mode='extend', extendpols=True, action="apply", flagbackup=False)

    plotms(vis=inputvis, field=fields_str, xaxis='frequency',yaxis='amplitude', ydatacolumn='corrected', correlation='RR,LL', showgui=False, plotfile='afterflagging_parhands.png', overwrite=True)
    # for cross-hands
    plotms(vis=inputvis, field=fields_str, xaxis='frequency', ydatacolumn='corrected', yaxis='amplitude', correlation='RL,LR', showgui=False, plotfile='afterflagging_crosshands.png', overwrite=True)
    # Initialize polarization calibration object
    spw_table = queryTable(table=inputvis, query="SELECT REF_FREQUENCY FROM "+inputvis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
    min_freq = queryTable(table=inputvis, query="SELECT GMIN(CHAN_FREQ) AS FREQ_MIN FROM "+inputvis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
    max_freq = queryTable(table=inputvis, query="SELECT GMAX(CHAN_FREQ) AS FREQ_MAX FROM "+inputvis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")

    spw_ids = spw_table.rownumbers()
    spw_reffreqs = spw_table.getcol("REF_FREQUENCY")
    min_freqs = min_freq.getcol("FREQ_MIN")[0]
    max_freqs = max_freq.getcol("FREQ_MAX")[0]
    nu_0 = np.median(spw_reffreqs)

    print("Ref freqs: ", spw_reffreqs)
    print("Ref freq: ", nu_0)
    print("Min freq: ", min_freqs)
    print("Max freq: ", max_freqs)

    polcal_obj = PolCalibration(vis=inputvis, spw_ids=spw_ids, polanglefield=polanglefield, leakagefield=leakagefield, target=target, refant=refant, kcross_refant=krefant)

    # Create your model source
    p_source = PolarizedSource(knownSource=polsource)
    p_source2 = PolarizedSource(knownSource="3C84")
    polcal_obj.setModel(pol_source_object=p_source, field=polanglefield, nu_0=nu_0, nterms_angle=1, nterms_frac=1, nu_min=min_freqs, nu_max=max_freqs)

    fluxtable = polcal_obj.setModelFluxScale(pol_source_object=p_source2, field=leakagefield, gaintable=complexgaintable, referencefield=polanglefield, transferfield=leakagefield, nu_0=nu_0)

    kcrosstable = polcal_obj.solveCrossHandDelays(refantmode="strict", channels="13~115")

    leakagetable = polcal_obj.calibrateLeakage(solint='inf,0.5MHz', minsnr=3.0)
    polcal_obj.plotLeakage(plotdir="./")

    polangletable = polcal_obj.calibratePolAngle(solint='inf,0.5MHz', minsnr=3.0)

    polcal_obj.applySolutions()
    polcal_obj.finalPlots()
    mstransform(vis=inputvis, outputvis='A1314.ms', datacolumn='corrected', field=target)
    #mstransform(vis=inputvis, outputvis='A1314_statwt.ms', datacolumn='corrected', field=target)
    #statwt(vis='A1314_statwt.ms', datacolumn='data', timebin='0.001s')
