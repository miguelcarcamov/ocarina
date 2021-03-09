import os
import sys
import numpy as np
import logging
from applycal import applycal
from polcal import polcal
from applycal import applycal
from gaincal import gaincal
from setjy import setjy
from fluxscale import fluxscale
from rmtables import rmtables
from plotms import plotms
from flagdata import flagdata
from utils import queryTable
from __casac__.logsink import logsink

class PolCalibration(object):
    def __init__(self, vis="", spw_ids=np.array([]), polanglefield="", leakagefield="", target="", refant="", kcross_refant="", nu_0=None, nu_min=None, nu_max=None, old_VLA=False, level=logging.INFO, **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.kcrosstable=''
        self.leakagetable=''
        self.polangletable=''
        self.casalog = logsink()
        self.logger = logging.getLogger(self.__class__.__name__)
        self.casalog.origin(self.__class__.__name__)
        self.logger.info("Creating "+self.__class__.__name__)
        self.casalog.post("Creating "+self.__class__.__name__, "INFO")

        if(self.nu_0 == None):
            spw_table = queryTable(table=self.vis, query="SELECT REF_FREQUENCY FROM "+self.vis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
            spw_reffreqs = spw_table.getcol("REF_FREQUENCY")
            self.nu_0 = np.median(spw_reffreqs)

        if(self.nu_min == None):
            min_freq = queryTable(table=self.vis, query="SELECT GMIN(CHAN_FREQ) AS FREQ_MIN FROM "+self.vis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
            self.nu_min = min_freq.getcol("FREQ_MIN")[0]

        if(self.nu_max == None):
            max_freq = queryTable(table=self.vis, query="SELECT GMAX(CHAN_FREQ) AS FREQ_MAX FROM "+self.vis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
            self.nu_max = max_freq.getcol("FREQ_MAX")[0]

        if(len(self.spw_ids)==0):
            spw_table = queryTable(table=self.vis, query="SELECT REF_FREQUENCY FROM "+self.vis+"/SPECTRAL_WINDOW"+" WHERE !FLAG_ROW")
            self.spw_ids = spw_table.rownumbers()

        self.nspw = len(self.spw_ids)
        self.casalog.post("nspw: "+ str(self.nspw), "INFO")
        self.casalog.post("Ref freq: "+ str(self.nu_0), "INFO")
        self.casalog.post("Min freq: "+ str(self.nu_min), "INFO")
        self.casalog.post("Max freq: "+ str(self.nu_max), "INFO")

    def getNu_0(self):
        return self.nu_0

    def getNu_min(self):
        return self.nu_min

    def getNu_max(self):
        return self.nu_max

    def getSpwIds(self):
        return self.spw_ids

    def setUnknownModel(self, pol_source_object=None, field="", gaintable="", referencefield="", transferfield="", fitorder=1, usescratch=False):
        fluxtable = self.vis[:-3]+".F0"
        if os.path.exists(fluxtable): rmtables(fluxtable)
        # From fluxscale documentation we know that the coefficients are return from the natural log nu/nu_0 Taylor expansion
        fluxdict = fluxscale(vis=self.vis, fluxtable=fluxtable, caltable=gaintable, reference=referencefield, transfer=transferfield, fitorder=fitorder)
        print(fluxdict)
        coeffs = fluxdict['2']['spidx'].tolist()
        flux_at_nu_0 = fluxdict['2']['fitFluxd']
        ref_freq = fluxdict['2']['fitRefFreq']
        print("Ref_freq: ", ref_freq)
        print("Coeffs: ", coeffs)
        print("Flux at nu_0: ", flux_at_nu_0)
        print("Flux at nu_0: ", np.exp(coeffs[0]))

        coeffs.pop(0)
        pol_source_object.setCoeffs(coeffs)
        spec_idx = coeffs
    
        self.logger.info("Setting model of: "+pol_source_object.getName())
        self.logger.info("Field: "+ field)
        self.logger.info("Reference freq (GHz): "+ str(ref_freq))
        self.logger.info("I = "+str(flux_at_nu_0))

        self.casalog.post("Setting model of: "+pol_source_object.getName(), "INFO")
        self.casalog.post("Field: "+ field, "INFO")
        self.casalog.post("Reference freq (GHz): "+ str(ref_freq), "INFO")
        self.casalog.post("I = "+ str(flux_at_nu_0), "INFO")

        print("Alpha & Beta: ", spec_idx)
        source_dict = setjy(vis=self.vis, field=field, standard='manual', spw='', fluxdensity=[flux_at_nu_0,0,0,0], spix=spec_idx, reffreq=ref_freq, interpolation="nearest", scalebychan=True, usescratch=usescratch)
        print(source_dict)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RRamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RLamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RRphase_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RLphase_model.png', overwrite=True)
        return fluxtable

    def setKnownModel(self, pol_source_object = None, standard="Perley-Butler 2017", field="", epoch="2017", nterms_angle=3, nterms_frac=3, usescratch=False):

        # get spectral idx coeffs from VLA tables
        intensity, spec_idx = pol_source_object.getKnownSourceInformation(nu_0=self.nu_0, standard=standard, epoch=epoch)
        pol_angle_coeffs, pol_frac_coeffs = pol_source_object.getSourcePolInformation(nu_0=self.nu_0, nterms_angle=nterms_angle, nterms_frac=nterms_frac, nu_min=self.nu_min, nu_max=self.nu_max)
        # get intensity in reference frequency
        self.logger.info("Setting model of: "+pol_source_object.getName())
        self.logger.info("Field: "+ field)
        self.logger.info("Reference freq (GHz): "+ str(self.nu_0/1e9))
        self.logger.info("I = "+ str(intensity))

        self.casalog.post("Setting model of: "+pol_source_object.getName(), "INFO")
        self.casalog.post("Field: "+ field, "INFO")
        self.casalog.post("Reference freq (GHz): "+ str(self.nu_0/1e9), "INFO")
        self.casalog.post("I = "+ str(intensity), "INFO")

        print("Alpha & Beta: ", spec_idx)
        print("Pol fraction coeffs: ", pol_frac_coeffs)
        print("Pol angle coeffs: ", pol_angle_coeffs)
        source_dict = setjy(vis=self.vis, field=field, standard='manual', spw='', fluxdensity=[intensity,0,0,0], spix=spec_idx, reffreq=str(self.nu_0/1e9)+"GHz", polindex=pol_frac_coeffs, polangle=pol_angle_coeffs, interpolation="nearest", scalebychan=True, usescratch=usescratch)
        print(source_dict)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RRamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RLamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RRphase_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RLphase_model.png', overwrite=True)


    def solveCrossHandDelays(self, minsnr=3.0, solint='inf', combine='scan,spw', channels="", refantmode="flex"):
        self.logger.info("Solving Cross-hand Delays")
        self.logger.info("Vis: "+ self.vis)
        self.logger.info("Field: "+ self.polanglefield)
        self.logger.info("Refant: "+ self.refant)

        self.casalog.post("Solving Cross-hand Delays", "INFO")
        self.casalog.post("Vis: "+ self.vis, "INFO")
        self.casalog.post("Field: "+ self.polanglefield, "INFO")
        self.casalog.post("Refant: "+ self.refant, "INFO")
        caltable = self.vis[:-3]+".Kcross"
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]

        if(channels==""):
            spw =str(firstspw)+'~'+str(lastspw)
        else:
            spw =str(firstspw)+'~'+str(lastspw)+':'+channels

        self.logger.info("Spw: " + spw)
        self.casalog.post("Spw: " + spw, "INFO")
        gaincal(vis=self.vis, caltable=caltable, field=self.polanglefield, spw=spw, refant=self.kcross_refant, refantmode=refantmode, minsnr=minsnr, gaintype="KCROSS", solint=solint, combine=combine, calmode="ap", append=False, gaintable=[''], gainfield=[''], interp=[''], spwmap=[[]], parang=True)
        if not os.path.exists(caltable): sys.exit("Caltable was not created and cannot continue. Exiting...")
        plotms(vis=caltable, xaxis='frequency', yaxis='delay', antenna=self.kcross_refant, coloraxis='corr', showgui=False, plotfile=self.vis[:-3]+'.freqvsdelayKcross.png', overwrite=True)
        self.kcrosstable = caltable
        return caltable


    def calibrateLeakage(self, solint='inf', minsnr=3.0, poltype="Df", gaintable=[], gainfield=[], clipmin=0.0, clipmax=0.25, flagclip=True, interpmode='linear'):
        if(gaintable == []):
            if(self.kcrosstable == ""):
                gaintable=[]
            else:
                gaintable=[self.kcrosstable]
        self.logger.info("Leakage calibration")
        self.logger.info("Vis: "+ self.vis)
        self.logger.info("Field: "+ self.leakagefield)

        self.casalog.post("Leakage calibration", "INFO")
        self.casalog.post("Vis: "+ self.vis, "INFO")
        self.casalog.post("Field: "+ self.leakagefield, "INFO")

        print("Gain tables: ", gaintable)
        self.logger.info("Refant: "+ self.refant)
        self.casalog.post("Refant: "+ self.refant, "INFO")
        caltable = self.vis[:-3]+".D0"
        if(gainfield == []): gainfield=[''] * len(gaintable)
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]

        spw = str(firstspw)+'~'+str(lastspw)
        spwmap0 = [0] * self.nspw
        spwmap=[spwmap0, []]
        interp = [interpmode] * len(gaintable)
        if(self.old_VLA):
            spw = ''
            spwmap = []
            interp='nearest'

        self.logger.info("Spw: ", spw)
        self.casalog.post("Spw: ", spw, "INFO")
        print("Spwmap: ", spwmap)
        polcal(vis=self.vis, caltable=caltable, field=self.leakagefield, spw=spw, refant=self.refant, poltype=poltype, solint=solint, spwmap=spwmap, combine='scan', interp=interp, minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)

        if not os.path.exists(caltable): sys.exit("Caltable was not created and cannot continue. Exiting...")

        plotms(vis=caltable,xaxis='freq',yaxis='amp',
       iteraxis='antenna',coloraxis='corr', showgui=False, plotfile=self.vis[:-3]+'.D0.ampvsfreq.png', overwrite=True)

        plotms(vis=caltable,xaxis='chan',yaxis='phase',
               iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,-180,180], showgui=False, plotfile=self.vis[:-3]+'.D0.phasevschan.png', overwrite=True)

        plotms(vis=caltable,xaxis='chan',yaxis='phase',
               iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,-180,180], showgui=False, plotfile=self.vis[:-3]+'.D0.ampvsantenna.png', overwrite=True)

        plotms(vis=caltable,xaxis='antenna1',yaxis='amp',coloraxis='corr', showgui=False, plotfile=self.vis[:-3]+'.D0.ampvsantenna1.png', overwrite=True)

        if(flagclip):
            flagdata(vis=caltable, mode='clip', correlation='ABS_ALL', clipminmax=[clipmin, clipmax], datacolumn='CPARAM', clipoutside=True, action='apply', flagbackup=False, savepars=False)

        self.leakagetable = caltable
        return caltable

    def calibratePolAngle(self, solint='inf', minsnr=3.0, poltype="Xf", gaintable=[], gainfield=[], interpmode='linear'):
        if(gaintable == []):
            if(self.kcrosstable == ""):
                gaintable=[self.leakagetable]
            else:
                gaintable=[self.kcrosstable, self.leakagetable]
        self.logger.info("Polarization angle calibration")
        self.logger.info("Vis: "+ self.vis)
        self.logger.info("Field: "+ self.polanglefield)

        self.casalog.post("Polarization angle calibration", "INFO")
        self.casalog.post("Vis: "+ self.vis, "INFO")
        self.casalog.post("Field: "+ self.polanglefield, "INFO")

        print("Gain tables: ", gaintable)

        self.logger.info("Refant: "+ self.refant)
        self.casalog.post("Refant: "+ self.refant, "INFO")

        caltable = self.vis[:-3]+".X0"
        if(gainfield == []): gainfield=[''] * len(gaintable)
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]

        spw = str(firstspw)+'~'+str(lastspw)
        spwmap0 = [0] * self.nspw
        spwmap=[spwmap0, []]
        interp = [interpmode] * len(gaintable)
        if(self.old_VLA):
            spw = ''
            spwmap = []
            interp = 'nearest'

        self.logger.info("Spw: ", spw)
        self.casalog.post("Spw: ", spw, "INFO")
        print("Spwmap: ", spwmap)
        polcal(vis=self.vis, caltable=caltable, field=self.polanglefield, spw=spw, refant=self.refant, poltype=poltype, solint=solint, combine='scan', spwmap=spwmap, interp=interp, minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)

        if not os.path.exists(caltable): sys.exit("Caltable was not created and cannot continue. Exiting...")
        plotms(vis=caltable,xaxis='frequency',yaxis='phase',coloraxis='spw', showgui=False, plotfile=self.vis[:-3]+'.X0.phasevsfreq.png', overwrite=True)

        self.polangletable = caltable
        return caltable

    def plotLeakage(self, plotdir=""):
        plotms(vis=self.leakagetable, xaxis='antenna', yaxis='amp', plotfile=plotdir+self.vis[:-3]+'.D0.amp.png', showgui=False, overwrite=True)
        plotms(vis=self.leakagetable, xaxis='antenna', yaxis='phase', iteraxis='antenna', plotfile=plotdir+self.vis[:-3]+'.D0.phs.png', showgui=False, overwrite=True)
        plotms(vis=self.leakagetable, xaxis='antenna', yaxis='snr', showgui=False, plotfile=plotdir+self.vis[:-3]+'.D0.snr.png', overwrite=True)
        plotms(vis=self.leakagetable, xaxis='real', yaxis='imag', showgui=False, plotfile=plotdir+self.vis[:-3]+'.D0.cmplx.png', overwrite=True)

    def applySolutions2(self, gainfield=[], applymode="calflagstrict"):
        gaintable=[self.kcrosstable]
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        spw = str(firstspw)+'~'+str(lastspw)
        self.logger.info("Spw: "+ spw)
        self.casalog.post("Spw: "+ spw, "INFO")
        spwmap0 = [0] * self.nspw
        interp = [''] * len(gaintable)
        calwt = [False] * len(gaintable)
        if(gainfield == []): gainfield = ['']
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintable, spwmap=[spwmap0], calwt=calwt, applymode=applymode, interp=interp, gainfield=gainfield, antenna='*&*', parang=True, flagbackup=True)

    def applySolutions3(self, gainfield=[], applymode="calflagstrict"):
        gaintable=[self.kcrosstable, self.leakagetable]
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        spw = str(firstspw)+'~'+str(lastspw)
        self.logger.info("Spw: "+ spw)
        self.casalog.post("Spw: "+ spw, "INFO")
        spwmap0 = [0] * self.nspw
        interp = [''] * len(gaintable)
        calwt = [False] * len(gaintable)
        if(gainfield == []): gainfield = ['']
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintable, spwmap=[spwmap0], calwt=calwt, applymode=applymode, interp=interp, gainfield=gainfield, antenna='*&*', parang=True, flagbackup=True)

    def applySingleSolution(self, field='', spw='', gaintable=[], gainfield=[], selectdata=True, spwmap=[], calwt=[False], applymode="calflagstrict", interp='linear', antenna='', flagbackup=True):
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintable, gainfield=gainfield, selectdata=selectdata, spwmap=spwmap, calwt=calwt, applymode=applymode, interp=interp, antenna=antenna, parang=True, flagbackup=flagbackup)

    def applySolutions(self, gaintable=[], gainfield=[], applymode="calflagstrict", antenna='*&*', flagbackup=True):
        #leakagegain.append(fluxtable)
        if(gaintable == []):
            if(self.kcrosstable == ""):
                gaintable=[self.leakagetable, self.polangletable]
            else:
                gaintable=[self.kcrosstable, self.leakagetable, self.polangletable]
        self.logger.info("Applying solutions")
        self.casalog.post("Applying solutions", "INFO")
        print("Gain tables: ", gaintable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]

        interp = [''] * len(gaintable)
        spw = str(firstspw)+'~'+str(lastspw)
        spwmap0 = [0] * self.nspw
        spwmap = [spwmap0, [], []]
        calwt = [False] * len(gaintable)
        selectdata=True
        if(self.old_VLA):
            interp = 'nearest'
            spwmap = []
            spw = ''
            calwt = [False]
            selectdata=False
            antenna=''

        self.logger.info("Spw: "+ spw)
        self.casalog.post("Spw: "+ spw, "INFO")
        print("Spwmap: ", spwmap)
        if(gainfield == []): gainfield = [''] * len(gaintable)
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintable, selectdata=selectdata, spwmap=spwmap, calwt=calwt, applymode=applymode, interp=interp, gainfield=gainfield, antenna=antenna, parang=True, flagbackup=flagbackup)

    def finalPlots(self):
        plotms(vis=self.vis, field=self.polanglefield, correlation='',
               timerange='',antenna='',avgtime='60',
               xaxis='frequency',yaxis='amp',ydatacolumn='corrected',
               coloraxis='corr',
               plotfile=self.polanglefield+'.corrected-amp.png', showgui=False, overwrite=True)

        plotms(vis=self.vis, field=self.polanglefield, correlation='',
               timerange='',antenna='',avgtime='60',
               xaxis='frequency',yaxis='phase',ydatacolumn='corrected',
               plotrange=[-1,-1,-180,180],coloraxis='corr',
               plotfile=self.polanglefield+'.corrected-phase.png', showgui=False, overwrite=True)

        plotms(vis=self.vis, field=self.leakagefield, correlation='',
               timerange='',antenna='',avgtime='60',
               xaxis='frequency',yaxis='amp',ydatacolumn='corrected', coloraxis='corr',
               plotfile=self.leakagefield+'.corrected-amp.png', showgui=False, overwrite=True)

        plotms(vis=self.vis, field=self.leakagefield, correlation='RR,LL',
               timerange='',antenna='',avgtime='60',
               xaxis='frequency',yaxis='phase',ydatacolumn='corrected',
               plotrange=[-1,-1,-180,180],coloraxis='corr',
               plotfile=self.leakagefield+'.corrected-phase.png', showgui=False, overwrite=True)
