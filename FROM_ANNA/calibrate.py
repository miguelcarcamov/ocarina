# -------------------------------------------------------------------------
# History:
# [AMS - 130718] created
# [AMS - 200515] updated to include polarization calibration
# -------------------------------------------------------------------------

import os,sys

print(" ")
print("----------------------------------------------------")
print("  This script is for processing Historic VLA data   ")
print("----------------------------------------------------")
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

nocal=True

if nocal:
    # ------------------------------------------------------------------------------------
    # import the VLA data file:
    os.system('rm -rf '+visname+'* \n') # this line removes the file if it already exists

    importvla(archivefiles=origvis,vis=tmpname)

    # ------------------------------------------------------------------------------------
    # split out only the useful data:
    split(vis=tmpname,outputvis=visname,datacolumn='data',spw=spw)

    # ------------------------------------------------------------------------------------
    # output the header info:
    listname=visname[:-2]+"list"
    listobs(vis=visname,listfile=listname,overwrite=True)

    # ------------------------------------------------------------------------------------
    # plot the configuration:
    plotants(vis=visname, figfile=plotdir+visname[:-3]+'.ants.png', showgui=False)

    # ------------------------------------------------------------------------------------
    # flagging:
    flagdata(vis=visname, mode='quack', quackinterval=10.0, quackmode='beg')
    #flagdata(vis=visname, mode='manual', antenna='VA09&&VA20')

# ------------------------------------------------------------------------------------
# set the flux scale:
setjy(vis=visname, standard='Perley-Butler 2013', field=bcal, usescratch=True)

# ------------------------------------------------------------------------------------
# generate the gain curve calibration table:
if os.path.exists('GC'): rmtables('GC')
    
gencal(vis=visname, caltype='gc', caltable='GC')
    
# ------------------------------------------------------------------------------------
# do the complex gain calibration:
if os.path.exists('G0'): rmtables('G0')

callist = ','.join([bcal,lcal,pcal])
gaincal(vis=visname, caltable='G0', gaintype='G', combine='', calmode='ap', field=callist, solint='inf', refant=refa, minsnr=3, gaintable=['GC'], parang=False)
        
plotcal(caltable='G0',xaxis='time', yaxis='amp', showgui=False, figfile=plotdir+visname[:-3]+'.G0.amp.png')
plotcal(caltable='G0',xaxis='time', yaxis='phase', showgui=False, figfile=plotdir+visname[:-3]+'.G0.phs.png')

# ------------------------------------------------------------------------------------
# transfer flux scale to phase cal:
if os.path.exists('F0'): rmtables('F0')
            
myFluxes = fluxscale(vis=visname, caltable='G0', reference=bcal, transfer=callist, fluxtable='F0', incremental=True, display=False, append=False)


