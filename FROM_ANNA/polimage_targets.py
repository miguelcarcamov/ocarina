# ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------
# History:
# [AMS - 130718] created
# [AMS - 200515] updated to include polarization calibration
# -------------------------------------------------------------------------


print(" ")
print("------------------------------------------------------")
print("  This script is for processing Historical VLA data   ")
print("------------------------------------------------------")
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
srclist  = ['1331+305','0521+166','1219+484','1133+490']

plotdir = './calplots/'
if not os.path.exists(plotdir): os.mkdir(plotdir)

imgdir = './images/'
if not os.path.exists(imgdir): os.mkdir(imgdir)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# MAIN SCRIPT:

casalog.filter('SEVERE')

# -------------------------------------------------------------------------
# last minute flagging:
#flagdata(vis=visname, field='0919+334,0915+320,0921+368,0926+328', mode="manual", antenna="VA08&&VA18")
#flagdata(vis=visname, field='0919+334,0915+320,0921+368,0926+328', mode="manual", antenna="VA08&&VA13")

for src in srclist:

    targ = src
    print("Imaging: ",targ)
    # -------------------------------------------------------------------------
    # apply the calibration:
    applycal(vis=visname, field=targ, spw=spw, gaintable=['GC','G0','F0','D0','X0'],gainfield=['',pcal,pcal,'',''], calwt=[False], parang=True)

    # ------------------------------------------------------------------------------------
    # make an image of the target:
    imagename = targ+".img0"
    os.system('rm -rf '+imagename+'* \n')

    print('>>> Processing target SPW 0 calibrator')
    tclean(vis=visname,
           stokes='IQUV',
           field=targ,
           spw='0',
           imagename=imagename,
           robust= 0.0,
           imsize=[512],
           cell=['2arcsec'],
           weighting='briggs',
           mask='',
           niter=500,
           interactive=False,
           pblimit=-0.2)


    # Create individual Stokes images:
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.I.image',stokes='I')
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.Q.image',stokes='Q')
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.U.image',stokes='U')
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.V.image',stokes='V')
    
    imgstats = imstat(imagename=imagename+'.I.image', box='216,1,255,72')
    sig_i = imgstats['rms'][0]
    print("Stokes I RMS:", imgstats['rms'][0])
    
    imgstats = imstat(imagename=imagename+'.Q.image', box='216,1,255,72')
    sig_q = imgstats['rms'][0]
    print("Stokes Q RMS:", imgstats['rms'][0])
    
    imgstats = imstat(imagename=imagename+'.U.image', box='216,1,255,72')
    sig_u = imgstats['rms'][0]
    print("Stokes U RMS:", imgstats['rms'][0])
    
    imgstats = imstat(imagename=imagename+'.V.image', box='216,1,255,72')
    sig_v = imgstats['rms'][0]
    print("Stokes V RMS:", imgstats['rms'][0])
    
    #  Create polarized intensity image:
    sig_p = 0.5*(sig_q+sig_u)
    sigma = str(sig_p)+"Jy/beam"
    immath(imagename=[imagename+'.Q.image',imagename+'.U.image'], mode='poli', sigma=sigma, outfile=imagename+'.POLI.image')
    
    polithresh = 8.*sig_p
    print("Map RMS for de-biasing:", sig_p)
    print("8 sigma polarization threshold:", polithresh)
    
    # Create polarization angle image:
    immath(imagename=[imagename+'.Q.image',imagename+'.U.image'], mode='pola', outfile=imagename+'.POLA.image',polithresh=str(polithresh)+'Jy/beam')
    


    # ------------------------------------------------------------------------------------
    # make an image of the target:
    imagename = targ+".img1"
    os.system('rm -rf '+imagename+'* \n')

    print('>>> Processing target SPW 1 calibrator')
    tclean(vis=visname,
           stokes='IQUV',
           field=targ,
           spw='1',
           imagename=imagename,
           robust= 0.0,
           imsize=[512],
           cell=['2arcsec'],
           weighting='briggs',
           mask='',
           niter=500,
           interactive=False,
           pblimit=-0.2)


    # Create individual Stokes images:
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.I.image',stokes='I')
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.Q.image',stokes='Q')
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.U.image',stokes='U')
    imsubimage(imagename=imagename+'.image',outfile=imagename+'.V.image',stokes='V')

    imgstats = imstat(imagename=imagename+'.I.image', box='216,1,255,72')
    sig_i = imgstats['rms'][0]
    print("Stokes I RMS:", imgstats['rms'][0])

    imgstats = imstat(imagename=imagename+'.Q.image', box='216,1,255,72')
    sig_q = imgstats['rms'][0]
    print("Stokes Q RMS:", imgstats['rms'][0])

    imgstats = imstat(imagename=imagename+'.U.image', box='216,1,255,72')
    sig_u = imgstats['rms'][0]
    print("Stokes U RMS:", imgstats['rms'][0])

    imgstats = imstat(imagename=imagename+'.V.image', box='216,1,255,72')
    sig_v = imgstats['rms'][0]
    print("Stokes V RMS:", imgstats['rms'][0])

    #  Create polarized intensity image:
    sig_p = 0.5*(sig_q+sig_u)
    sigma = str(sig_p)+"Jy/beam"
    immath(imagename=[imagename+'.Q.image',imagename+'.U.image'], mode='poli', sigma=sigma, outfile=imagename+'.POLI.image')
    
    polithresh = 8.*sig_p
    print("Map RMS for de-biasing:", sig_p)
    print("8 sigma polarization threshold:", polithresh)
    
    # Create polarization angle image:
    immath(imagename=[imagename+'.Q.image',imagename+'.U.image'], mode='pola', outfile=imagename+'.POLA.image',polithresh=str(polithresh)+'Jy/beam')
    


    # ------------------------------------------------------------------------------------
    # clean up:

    os.system('rm -rf images/'+imagename[:-1]+'* \n')
    os.system('mv '+imagename[:-1]+'* images/ \n')
