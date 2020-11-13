import os
import numpy as np
from applycal import applycal
from plotcal import plotcal
from polcal import polcal
from applycal import applycal
from gaincal import gaincal
from setjy import setjy
from fluxscale import fluxscale
from rmtables import rmtables
from plotms import plotms

class PolCalibration(object):
    def __init__(self, vis="", spw_ids=np.array([]), polanglefield="", leakagefield="", target="", refant="", calibration_tables=[], **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.nspw = len(spw_ids)

    def setModelFluxScale(self, field="", gaintable="", referencefield="", transferfield="", nu_0=0.0):
        fluxtable = self.vis[:-3]+".F0"
        if os.path.exists(fluxtable): rmtables(fluxtable)
        fluxdict = fluxscale(vis=self.vis, fluxtable=fluxtable, caltable=gaintable, reference=referencefield, transfer=transferfield)
        spec_idx = fluxdict['2']['spidx'].tolist()
        I_nu_0 = fluxdict['2']['fitFluxd']
        reffreq = fluxdict['2']['fitRefFreq']
        intensity_at_nu_0 = I_nu_0 * (nu_0/reffreq)**(spec_idx[0]+spec_idx[1]*np.log(nu_0/reffreq))
        print("Intensity: ", intensity_at_nu_0)
        print("Alpha & Beta: ", spec_idx)
        setjy(vis=self.vis, field=field, standard='manual', spw='', fluxdensity=[intensity_at_nu_0,0,0,0], spix=spec_idx, reffreq=str(nu_0/1e9)+"GHz", polindex=[], polangle=[], rotmeas=0, interpolation="nearest", scalebychan=True, usescratch=False)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RRamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RLamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RLphase_model.png', overwrite=True)
        return fluxtable

    def setModel(self, pol_source_object = None, standard="Perley-Butler 2013", field="", epoch="2012", nu_0=0.0, nterms_angle=3, nterms_frac=3, nu_min=0.0, nu_max=np.inf):

        # get spectral idx coeffs from VLA tables
        pol_source_object.getCoeffs(standard=standard, epoch=epoch)
        nu_fit = np.linspace(0.3275*1e9, 50.0*1e9, 40)
        spec_idx = pol_source_object.fitAlphaBeta(nu_fit, nu_0=nu_0)
        pol_frac_coeffs = pol_source_object.getPolFracCoeffs(nu_0=nu_0, nterms=nterms_frac, nu_min=nu_min, nu_max=nu_max)
        pol_angle_coeffs = pol_source_object.getPolAngleCoeffs(nu_0=nu_0, nterms=nterms_angle, nu_min=nu_min, nu_max=nu_max)
        # get intensity in reference frequency
        intensity = pol_source_object.flux_scalar(nu_0/1e9)
        print("Setting model of: ", pol_source_object.getName())
        print("Reference freq (GHz): ", nu_0/1e9)
        print("I = ", intensity)
        print("Alpha & Beta: ", spec_idx)
        print("Pol fraction coeffs: ", pol_frac_coeffs)
        print("Pol angle coeffs: ", pol_angle_coeffs)
        setjy(vis=self.vis, field=field, standard='manual', spw='', fluxdensity=[intensity,0,0,0], spix=spec_idx, reffreq=str(nu_0/1e9)+"GHz", polindex=pol_frac_coeffs, polangle=pol_angle_coeffs, interpolation="nearest", scalebychan=True, usescratch=False)
        plotms(vis=self.vis, field=field, correlation='RR', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RRamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='amp', ydatacolumn='model', showgui=False, plotfile=field+'_RLamp_model.png', overwrite=True)
        plotms(vis=self.vis, field=field, correlation='RL', timerange='', antenna=self.refant, xaxis='frequency', yaxis='phase', ydatacolumn='model', showgui=False, plotfile=field+'_RLphase_model.png', overwrite=True)

    def solveCrossHandDelays(self):
        print("Solving Cross-hand Delays")
        print("Vis: ", self.vis)
        print("Field: ", self.polanglefield)
        print("Refant: ", self.refant)
        caltable = self.vis[:-3]+".Kcross"
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        print("Spw: ", str(firstspw)+'~'+str(lastspw))
        gaincal(vis=self.vis, caltable=caltable, field=self.polanglefield, spw=str(firstspw)+'~'+str(lastspw), refant=self.refant, gaintype="KCROSS", solint="inf", combine="scan,spw", calmode="ap", append=False, gaintable=[''], gainfield=[''], interp=[''], spwmap=[[]], parang=True)
        plotcal(caltable=caltable, xaxis='freq', yaxis='delay', antenna=self.refant, showgui=False, figfile=self.vis[:-3]+'.freqvsdelayKcross.png')
        return caltable


    def calibrateLeakage(self, minsnr=3.0, poltype="D", gaintable=[], gainfield=[]):
        print("Leakage calibration")
        print("Vis: ", self.vis)
        print("Field: ", self.leakagefield)
        print("Gain tables: ", gaintable)
        print("Refant: ", self.refant)
        caltable = self.vis[:-3]+".D0"
        if(gainfield == []): gainfield=['']
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        print("Spw: ", str(firstspw)+'~'+str(lastspw))
        spwmap = [0] * self.nspw
        polcal(vis=self.vis, caltable=caltable, field=self.leakagefield, spw=str(firstspw)+'~'+str(lastspw), refant=self.refant, poltype=poltype, solint='inf', spwmap=spwmap, combine='scan', minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)
        plotcal(caltable=caltable,xaxis='freq',yaxis='amp', iteration='antenna', showgui=False, figfile=self.vis[:-3]+'.D0.ampvsfreq.png')

        plotcal(caltable=caltable,xaxis='chan',yaxis='phase', iteration='antenna',plotrange=[-1,-1,-180,180], showgui=False, figfile=self.vis[:-3]+'.D0.phasevschan.png')

        plotcal(caltable=caltable,xaxis='antenna',yaxis='amp', showgui=False, figfile=self.vis[:-3]+'.D0.ampvsantenna.png')
        return caltable

    def calibratePolAngle(self, minsnr=3.0, poltype="Xf", gaintable=[], gainfield=[]):
        print("Polarization angle calibration")
        print("Vis: ", self.vis)
        print("Field: ", self.polanglefield)
        print("Gain tables: ", gaintable)
        print("Refant: ", self.refant)
        caltable = self.vis[:-3]+".X0"
        if(gainfield == []): gainfield=['']*len(gaintable)
        if os.path.exists(caltable): rmtables(caltable)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        print("Spw: ", str(firstspw)+'~'+str(lastspw))
        spwmap0 = [0] * self.nspw
        polcal(vis=self.vis, caltable=caltable, field=self.polanglefield, spw=str(firstspw)+'~'+str(lastspw), refant=self.refant, poltype=poltype, solint='inf', combine='scan', spwmap=[spwmap0, []], minsnr=minsnr, gaintable=gaintable, gainfield=gainfield)
        plotcal(caltable=caltable, xaxis='freq', yaxis='phase', showgui=False, figfile=self.vis[:-3]+'.X0.phasevsfreq.png')
        return caltable

    def plotLeakage(self, plotdir=""):
        leakagecaltable = self.vis[:-3]+".D0"
        plotcal(caltable=leakagecaltable, xaxis='antenna', yaxis='amp', figfile=plotdir+self.vis[:-3]+'.D0.amp.png', showgui=False)
        plotcal(caltable=leakagecaltable, xaxis='antenna', yaxis='phase', iteration='antenna', figfile=plotdir+self.vis[:-3]+'.D0.phs.png', showgui=False)
        plotcal(caltable=leakagecaltable, xaxis='antenna', yaxis='snr', showgui=False, figfile=plotdir+self.vis[:-3]+'.D0.snr.png')
        plotcal(caltable=leakagecaltable, xaxis='real', yaxis='imag', showgui=False, figfile=plotdir+self.vis[:-3]+'.D0.cmplx.png')

    def applySolutions(self, gaintables=[]):
        #leakagegain.append(fluxtable)
        print("Applying solutions")
        print("Gain tables: ", gaintables)
        firstspw=self.spw_ids[0]
        lastspw=self.spw_ids[-1]
        spw = str(firstspw)+'~'+str(lastspw)
        print("Spw: ", spw)
        spwmap0 = [0] * self.nspw
        interp = ['linear'] * len(gaintables)
        calwt = [False] * len(gaintables)
        gainfield = [''] *  len(gaintables)
        applycal(vis=self.vis, field='', spw=spw, gaintable=gaintables, spwmap=[spwmap0, [], []], calwt=calwt, applymode='calflagstrict', interp=interp, gainfield=gainfield, antenna='*&*', parang=True, flagbackup=True)

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
