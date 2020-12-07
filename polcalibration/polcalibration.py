import os
import sys
import numpy as np
import logging
from applycal import applycal
from plotcal import plotcal
from polcal import polcal
from applycal import applycal
from gaincal import gaincal
from setjy import setjy
from fluxscale import fluxscale
from rmtables import rmtables
from plotms import plotms
from flagdata import flagdata

class PolCalibration(object):
    def __init__(self, vis="", spw_ids=np.array([]), polanglefield="", leakagefield="", target="", refant="", kcross_refant="", old_VLA=False, level=logging.INFO, **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.nspw = len(spw_ids)
        self.kcrosstable=''
        self.leakagetable=''
        self.polangletable=''
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.info("Creating "+self.__class__.__name__)

    def setModelFluxScale(self, pol_source_object=None, field="", gaintable="", referencefield="", transferfield="", nu_0=0.0):
        fluxtable = self.vis[:-3]+".F0"
        if os.path.exists(fluxtable): rmtables(fluxtable)
        # From fluxscale documentation we know that the returned coefficients come from the log10 Taylor expansion
        fluxdict = fluxscale(vis=self.vis, fluxtable=fluxtable, caltable=gaintable, reference=referencefield, transfer=transferfield, fitorder=4)
        print(fluxdict)
        coeffs = fluxdict['2']['spidx'].tolist()
        pol_source_object.setCoeffs(coeffs)
        nu_fit = np.linspace(0.3275*1e9, 50.0*1e9, 40)
        spec_idx = pol_source_object.fitAlphaBeta(nu_fit, nu_0=nu_0)
        intensity = pol_source_object.flux_scalar(nu_0/1e9)
        self.logger.info("Setting model of: "+pol_source_object.getName())
        self.logger.info("Field: "+ field)
        self.logger.info("Reference freq (GHz): "+ str(nu_0/1e9))
        self.logger.info("I = "+ str(intensity))
        print("Alpha & Beta: ", spec_idx)
        source_dict = setjy(vis=self.vis, field=field, standard='manual', spw='', fluxdensity=[intensity,0,0,0], spix=spec_idx, reffreq=str(nu_0/1e9)+"GHz", interpolation="nearest", scalebychan=True, usescratch=True)
        print(source_dict)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RRamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RLamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RRphase_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RLphase_model.png', overwrite=True)
        return fluxtable

    def setModel(self, pol_source_object = None, standard="Perley-Butler 2017", field="", epoch="2017", nu_0=0.0, nterms_angle=3, nterms_frac=3, nu_min=0.0, nu_max=np.inf):

        # get spectral idx coeffs from VLA tables
        pol_source_object.getCoeffs(standard=standard, epoch=epoch)
        nu_fit = np.linspace(0.3275*1e9, 50.0*1e9, 40)
        spec_idx = pol_source_object.fitAlphaBeta(nu_fit, nu_0=nu_0)
        pol_frac_coeffs = pol_source_object.getPolFracCoeffs(nu_0=nu_0, nterms=nterms_frac, nu_min=nu_min, nu_max=nu_max)
        pol_angle_coeffs = pol_source_object.getPolAngleCoeffs(nu_0=nu_0, nterms=nterms_angle, nu_min=nu_min, nu_max=nu_max)
        # get intensity in reference frequency
        intensity = pol_source_object.flux_scalar(nu_0/1e9)
        self.logger.info("Setting model of: "+pol_source_object.getName())
        self.logger.info("Field: "+ field)
        self.logger.info("Reference freq (GHz): "+ str(nu_0/1e9))
        self.logger.info("I = "+ str(intensity))
        print("Alpha & Beta: ", spec_idx)
        print("Pol fraction coeffs: ", pol_frac_coeffs)
        print("Pol angle coeffs: ", pol_angle_coeffs)
        source_dict = setjy(vis=self.vis, field=field, standard='manual', spw='', fluxdensity=[intensity,0,0,0], spix=spec_idx, reffreq=str(nu_0/1e9)+"GHz", polindex=pol_frac_coeffs, polangle=pol_angle_coeffs, interpolation="nearest", scalebychan=True, usescratch=True)
        print(source_dict)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RRamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RLamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RRphase_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RLphase_model.png', overwrite=True)


    def solveCrossHandDelays(self, solint='inf', combine='scan,spw', channels="", refantmode="flex"):
        self.logger.info("Solving Cross-hand Delays")
        self.logger.info("Vis: "+ self.vis)
        self.logger.info("Field: "+ self.polanglefield)
        self.logger.info("Refant: "+ self.refant)
        caltable = self.vis[:-3]+".Kcross"
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]

        if(channels==""):
            spw =str(firstspw)+'~'+str(lastspw)
        else:
            spw =str(firstspw)+'~'+str(lastspw)+':'+channels

        self.logger.info("Spw: " + spw)
        gaincal(vis=self.vis, caltable=caltable, field=self.polanglefield, spw=spw, refant=self.kcross_refant, refantmode=refantmode, gaintype="KCROSS", solint=solint, combine=combine, calmode="ap", append=False, gaintable=[''], gainfield=[''], interp=[''], spwmap=[[]], parang=True)
        if not os.path.exists(caltable): sys.exit("Caltable was not created and cannot continue. Exiting...")
        plotcal(caltable=caltable, xaxis='freq', yaxis='delay', antenna=self.refant, showgui=False, figfile=self.vis[:-3]+'.freqvsdelayKcross.png')
        self.kcrosstable = caltable
        return caltable


    def calibrateLeakage(self, solint='inf', minsnr=3.0, poltype="Df", gainfield=[], clipmin=0.0, clipmax=0.25, flagclip=True):
        if(self.kcrosstable == ""):
            gaintable=[]
        else:
            gaintable=[self.kcrosstable]
        self.logger.info("Leakage calibration")
        self.logger.info("Vis: "+ self.vis)
        self.logger.info("Field: "+ self.leakagefield)
        print("Gain tables: ", gaintable)
        self.logger.info("Refant: "+ self.refant)
        caltable = self.vis[:-3]+".D0"
        if(gainfield == []): gainfield=[''] * len(gaintable)
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        self.logger.info("Spw: ", str(firstspw)+'~'+str(lastspw))
        spwmap = [0] * self.nspw
        polcal(vis=self.vis, caltable=caltable, field=self.leakagefield, spw=str(firstspw)+'~'+str(lastspw), refant=self.refant, poltype=poltype, solint=solint, spwmap=spwmap, combine='scan', minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)

        if not os.path.exists(caltable): sys.exit("Caltable was not created and cannot continue. Exiting...")

        if(flagclip):
            flagdata(vis=caltable, mode='clip', correlation='ABS_ALL', clipminmax=[clipmin, clipmax], datacolumn='CPARAM', clipoutside=True, action='apply', flagbackup=False, savepars=False)

        plotcal(caltable=caltable,xaxis='freq',yaxis='amp', iteration='antenna', showgui=False, figfile=self.vis[:-3]+'.D0.ampvsfreq.png')

        plotcal(caltable=caltable,xaxis='chan',yaxis='phase', iteration='antenna',plotrange=[-1,-1,-180,180], showgui=False, figfile=self.vis[:-3]+'.D0.phasevschan.png')

        plotcal(caltable=caltable,xaxis='antenna',yaxis='amp', showgui=False, figfile=self.vis[:-3]+'.D0.ampvsantenna.png')
        self.leakagetable = caltable
        return caltable

    def calibratePolAngle(self, solint='inf', minsnr=3.0, poltype="Xf", gainfield=[]):
        if(self.kcrosstable == ""):
            gaintable=[self.leakagetable]
        else:
            gaintable=[self.kcrosstable, self.leakagetable]
        self.logger.info("Polarization angle calibration")
        self.logger.info("Vis: "+ self.vis)
        self.logger.info("Field: "+ self.polanglefield)
        print("Gain tables: ", gaintable)
        self.logger.info("Refant: "+ self.refant)
        caltable = self.vis[:-3]+".X0"
        if(gainfield == []): gainfield=[''] * len(gaintable)
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        self.logger.info("Spw: ", str(firstspw)+'~'+str(lastspw))
        spwmap0 = [0] * self.nspw
        polcal(vis=self.vis, caltable=caltable, field=self.polanglefield, spw=str(firstspw)+'~'+str(lastspw), refant=self.refant, poltype=poltype, solint=solint, combine='scan', spwmap=[spwmap0, []], minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)

        if not os.path.exists(caltable): sys.exit("Caltable was not created and cannot continue. Exiting...")

        plotcal(caltable=caltable, xaxis='freq', yaxis='phase', showgui=False, figfile=self.vis[:-3]+'.X0.phasevsfreq.png')
        self.polangletable = caltable
        return caltable

    def plotLeakage(self, plotdir=""):
        plotcal(caltable=self.leakagetable, xaxis='antenna', yaxis='amp', figfile=plotdir+self.vis[:-3]+'.D0.amp.png', showgui=False)
        plotcal(caltable=self.leakagetable, xaxis='antenna', yaxis='phase', iteration='antenna', figfile=plotdir+self.vis[:-3]+'.D0.phs.png', showgui=False)
        plotcal(caltable=self.leakagetable, xaxis='antenna', yaxis='snr', showgui=False, figfile=plotdir+self.vis[:-3]+'.D0.snr.png')
        plotcal(caltable=self.leakagetable, xaxis='real', yaxis='imag', showgui=False, figfile=plotdir+self.vis[:-3]+'.D0.cmplx.png')

    def applySolutions2(self, gainfield=[], applymode="calflagstrict"):
        gaintables=[self.kcrosstable]
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        spw = str(firstspw)+'~'+str(lastspw)
        self.logger.info("Spw: "+ spw)
        spwmap0 = [0] * self.nspw
        interp = [''] * len(gaintables)
        calwt = [False] * len(gaintables)
        if(gainfield == []): gainfield = ['']
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintables, spwmap=[spwmap0], calwt=calwt, applymode=applymode, interp=interp, gainfield=gainfield, antenna='*&*', parang=True, flagbackup=True)

    def applySolutions3(self, gainfield=[], applymode="calflagstrict"):
        gaintables=[self.kcrosstable, self.leakagetable]
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        spw = str(firstspw)+'~'+str(lastspw)
        self.logger.info("Spw: "+ spw)
        spwmap0 = [0] * self.nspw
        interp = [''] * len(gaintables)
        calwt = [False] * len(gaintables)
        if(gainfield == []): gainfield = ['']
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintables, spwmap=[spwmap0], calwt=calwt, applymode=applymode, interp=interp, gainfield=gainfield, antenna='*&*', parang=True, flagbackup=True)


    def applySolutions(self, gainfield=[], applymode="calflagstrict"):
        #leakagegain.append(fluxtable)
        if(self.kcrosstable == ""):
            gaintables=[self.leakagetable, self.polangletable]
        else:
            gaintables=[self.kcrosstable, self.leakagetable, self.polangletable]
        self.logger.info("Applying solutions")
        print("Gain tables: ", gaintables)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        if(self.old_VLA):
            interp = ['nearest'] * len(gaintables)
            spwmap = []
            spw = ''
            calwt = [False]
        else:
            interp = [''] * len(gaintables)
            spw = str(firstspw)+'~'+str(lastspw)
            spwmap0 = [0] * self.nspw
            spwmap = [spwmap0, [], []]
            calwt = [False] * len(gaintables)
        self.logger.info("Spw: "+ spw)
        if(gainfield == []): gainfield = ['', '', '']
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintables, spwmap=spwmap, calwt=calwt, applymode=applymode, interp=interp, gainfield=gainfield, antenna='*&*', parang=True, flagbackup=True)

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
