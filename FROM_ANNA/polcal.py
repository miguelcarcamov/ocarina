# -------------------------------------------------------------------------
# History:
# [AMS - 130718] created
# [AMS - 200515] updated to include polarization calibration
# -------------------------------------------------------------------------


print(" ")
print("--------------------------------------------")
print("  This script is for processing JVLA data   ")
print("--------------------------------------------")
print(" ")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# TARGET SPECIFIC:

origvis  = './AC461_A960805.xp1'
tmpname  = 'AC461.ms.0'
visname  = 'AC461.ms'
refa     = 'VA06'
bcal     = '1331+305'  # 3C286
lcal     = '0521+166'  # 3C138
pcal     = '1219+484'
spw      = '*'

plotdir = './calplots/'
if not os.path.exists(plotdir): os.mkdir(plotdir)

def p3c286(x):
    return 11.63 - 5.89*np.exp(-0.68*x)

def p3c138(x):
    f = np.array([1.365,1.465,1.435,1.665,4.5351,4.8851])
    p = np.array([7.1,7.5,7.5,8.4,10.0,10.0])
    return p[np.argmin(np.abs(f-x))]

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# MAIN SCRIPT:

# ------------------------------------------------------------------------------------
# 3C286:

nspw = 6
for rf in range(0,nspw-1):

    # get reference flux densities:
    calflux = setjy(vis=visname, standard='Perley-Butler 2013', field=bcal, spw=str(rf)+","+str(rf+1))
    
    # get flux densities:
    fid = calflux.keys()[0]
    i0 = calflux[fid][str(rf)]['fluxd'][0]
    i1 = calflux[fid][str(rf+1)]['fluxd'][0]

    # get reference frequencies:
    tb.open(visname + '/SPECTRAL_WINDOW')
    f0 = tb.getcol('REF_FREQUENCY')[rf]/1e9
    f1 = tb.getcol('REF_FREQUENCY')[rf+1]/1e9
    tb.close()

    alpha=-1.*log(i0/i1)/log(f1/f0) # Values from our setjy() run on Stokes I earlier
    c0=p3c286(f0)/100.  # Fractional polarization @ f0
    d0=33*pi/180 # Polarisation angle 33deg in radians [true from 1 - 6 GHz]

    setjy(vis=visname, field=bcal, standard='manual', spw=str(rf), fluxdensity=[i0,0,0,0], spix=[alpha,0], reffreq=str(f0)+'GHz', polindex=[c0,0], polangle=[d0,0], scalebychan=True, usescratch=False)

c0=p3c286(f1)/100.  # Fractional polarization @ f0
setjy(vis=visname, field=bcal, standard='manual', spw=str(nspw-1), fluxdensity=[i1,0,0,0], spix=[alpha,0], reffreq=str(f1)+'GHz', polindex=[c0,0], polangle=[d0,0], scalebychan=True, usescratch=False)

# ------------------------------------------------------------------------------------
# 3C138:

nspw = 6
for rf in range(0,nspw-1):
    
    # get reference flux densities:
    calflux = setjy(vis=visname, standard='Perley-Butler 2013', field=lcal, spw=str(rf)+","+str(rf+1))
    
    # get flux densities:
    fid = calflux.keys()[0]
    i0 = calflux[fid][str(rf)]['fluxd'][0]
    i1 = calflux[fid][str(rf+1)]['fluxd'][0]
    
    # get reference frequencies:
    tb.open(visname + '/SPECTRAL_WINDOW')
    f0 = tb.getcol('REF_FREQUENCY')[rf]/1e9
    f1 = tb.getcol('REF_FREQUENCY')[rf+1]/1e9
    tb.close()
    
    alpha=-1.*log(i0/i1)/log(f1/f0) # Values from our setjy() run on Stokes I earlier
    c0=p3c138(f0)/100.  # Fractional polarization @ f0
    d0=-10*pi/180 # Polarisation angle -10deg in radians [approx. true from 1 - 6 GHz]
    
    setjy(vis=visname, field=lcal, standard='manual', spw=str(rf), fluxdensity=[i0,0,0,0], spix=[alpha,0], reffreq=str(f0)+'GHz', polindex=[c0,0], polangle=[d0,0], scalebychan=True, usescratch=False)

c0=p3c138(f1)/100.  # Fractional polarization @ f0
setjy(vis=visname, field=lcal, standard='manual', spw=str(nspw-1), fluxdensity=[i1,0,0,0], spix=[alpha,0], reffreq=str(f1)+'GHz', polindex=[c0,0], polangle=[d0,0], scalebychan=True, usescratch=False)


# ------------------------------------------------------------------------------------
# leakage calibration (lcal):

if os.path.exists('D0'): rmtables('D0')

polcal(vis=visname, caltable='D0', field=lcal, spw='', refant=refa, poltype='D', solint='inf', combine='scan', minsnr=3, gaintable=['GC', 'G0'], gainfield=['',lcal])

plotcal(caltable='D0', xaxis='antenna', yaxis='amp', figfile=plotdir+visname[:-3]+'.D0.amp.png', showgui=False)
plotcal(caltable='D0', xaxis='antenna', yaxis='phase', iteration='antenna', figfile=plotdir+visname[:-3]+'.D0.phs.png', showgui=False)
plotcal(caltable='D0', xaxis='antenna', yaxis='snr', showgui=False, figfile=plotdir+visname[:-3]+'.D0.snr.png')
plotcal(caltable='D0', xaxis='real', yaxis='imag', showgui=False, figfile=plotdir+visname[:-3]+'.D0.cmplx.png')


# ------------------------------------------------------------------------------------
# polarization angle calibration (3C286):
if os.path.exists('X0'): rmtables('X0')

polcal(vis=visname, caltable='X0', field=bcal, spw='', refant=refa, poltype='Xf', solint='inf', combine='scan', minsnr=3, gaintable=['GC', 'G0', 'D0'], gainfield=['',bcal,''])


# ------------------------------------------------------------------------------------
# apply calibration solutions to check calibration so far:

applycal(vis=visname, field=bcal, spw=spw, gaintable=['GC','G0','D0','X0'],gainfield=['',bcal,'',''], calwt=[False], parang=True)

applycal(vis=visname, field=lcal, spw=spw, gaintable=['GC','G0','F0','D0','X0'],gainfield=['',lcal,lcal,'',''], calwt=[False], parang=True)

applycal(vis=visname, field=pcal, spw=spw, gaintable=['GC','G0','F0','D0','X0'],gainfield=['',pcal,pcal,'',''], calwt=[False], parang=True)

