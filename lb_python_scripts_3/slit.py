import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import pyfits
import sys
from numpy import *
from fit_continuum import *
from median import *

def onecam(file,snap,cam,max_lum=7.5,dynrng_lum=4.5,
           subtract_continuum=True,**kwargs):

    print "subtract_continuum:",subtract_continuum

    ## hbeta rest wavelength
    #lambda_rest = 4861.36
    lambda_rest = 4861.0

    maindir=''
    nonlr_dir=''
    wholesed_dir=''
    if kwargs is not None:
        #print kwargs.items()
        for name, value in kwargs.items():
            if name=='maindir':
                maindir=value+'/'
            elif name=='nonlr_dir':
                nonlr_dir=value+'/'
                print "WARNING! onecam doesn't have anything coded for adding nonlr runs yet."
                print "the keyword 'nonlr_dir' currently has no effect."
            elif name=='wholesed_dir':
                wholesed_dir=value+'/'
                print "WARNING! onecam doesn't have anything coded for adding wholesed runs yet."
                print "the keyword 'wholesed_dir' currently has no effect."
            else: print "kwarg ",name," is not defined."    

    shortdir = (maindir+file).split('lblecha/')
    if len(shortdir) == 2:
        shortdir = shortdir[1]
    else: shortdir=file

    filebase = ((maindir+file).split('mcrx'))[0]
    if filebase != '': filebase = filebase+'/'
    
    f=pyfits.open(maindir+file)    
    print "\nopened file: ",maindir+file


    init_lam = f['lambda'].data.field('lambda')

    ## clip ends of wavelength range
    #llim = array([25,init_lam.shape[0]-26])
    llim = array([0,init_lam.shape[0]-1])

    if max_lum < 2.0 or dynrng_lum < 1.0:
        print "crazy value entered for max_lum or dynrng_lum:", max_lum, dynrng_lum
        quit()

##   *NOTE*: pos and hfw are NOT being entered as kwargs.
##   need an elegant way to pass them to the 'process' function
##   and it doesn't make much sense to enter them at the command line
##   unless we're using physical units.
##   Thus they are being HARD-CODED for now to correspond to
##   2kpc around the central pixel in the image.
##
##     if kwargs is not None:
##         #print kwargs.items()
##         for name, value in kwargs.items():
##             if name is 'slit_pos':
##                 if int(value) >= ylim[0] and int(value) <= ylim[1]:
##                     pos = int(value)
##                 else:
##                     print "value for slit_pos (",int(value),")is outside of plot range:",ylim
##                     print "   ...using default value:",pos
##             elif name is 'slit_width':
##                 if pos-int(value) >= ylim[0] and pos+int(value) <= ylim[1]:
##                     hfw = int(value)
##                 else:
##                     print "value for slit_width (",int(value),")is outside of plot range:",ylim
##                     print "   ...using default value:",hfw
##             elif name is 'max_lum':
##                 if float(value) >= 2.0:
##                     max_lum = float(value)
##                 else: print "invalid value entered for max_lum:", float(value)
##             elif name is 'dynrng_lum':
##                 if float(value) >= 1.0:
##                     dynrng_lum = float(value)
##                 else: print "invalid value entered for dynrng_lum:", float(value)
##             else: print "kwarg ",name," is not defined."


    ## initial values of normalization factors
    sed_norm = 0.0
    sed_scatter_norm = 0.0
    img_vmax = 0.0
    img_vmin = 1.0e6
    img_scatter_vmax = 0.0
    img_scatter_vmin = 1.0e6
    lam_half_rng = 0.0

    print "\nProcessing camera %d..."%cam

    ## process this camera and return relevant values
    pos, hfw, xlim, ylim, kpc_per_pix, lam, sed_noscat, sed_scat, \
        slit_noscat, slit_scat,img_noscat, img_scat, ifu_noscat, ifu_scat, lowlum_sed_noscat = \
        process_onecam(f, snap, cam, init_lam, llim, lambda_rest, \
                           max_lum=max_lum, dynrng_lum=dynrng_lum)

    print "min/max/shape sed_noscat: ",sed_noscat.min(),sed_noscat.max(),sed_noscat.shape
    print "min/max/shape sed_scat: ",sed_scat.min(),sed_scat.max(),sed_scat.shape
    print "min/max/shape slit_noscat: ",slit_noscat.min(),slit_noscat.max(),slit_noscat.shape
    print "min/max/shape slit_scat: ",slit_scat.min(),slit_scat.max(),slit_scat.shape
    print "min/max/shape img_noscat: ",img_noscat.min(),img_noscat.max(),img_noscat.shape
    print "min/max/shape img_scat: ",img_scat.min(),img_scat.max(),img_scat.shape
    print "min/max/shape ifu_noscat: ",ifu_noscat.min(),ifu_noscat.max(),ifu_noscat.shape
    print "min/max/shape ifu_scat: ",ifu_scat.min(),ifu_scat.max(),ifu_scat.shape

    ## convert luminosities to solar units
    sed_noscat /= 3.83e33
    sed_scat /= 3.83e33
    slit_noscat -= log10(3.83e33)
    slit_scat -= log10(3.83e33)
    img_noscat -= log10(3.83e33)
    img_scat -= log10(3.83e33)

    ## subtract continuum from 1-D SED if necessary
    ## applies ONLY to 1-d SED!
    if subtract_continuum:
        
        f_continuum_ns, vfinal_sigma_ns = fit_continuum(lam,sed_noscat, 5, 1.0, 0.02e-7)
        f_continuum, vfinal_sigma = fit_continuum(lam,sed_scat, 5, 1.0, 0.02e-7)
        
        #print shape(f_continuum),shape(vfinal_sigma)            
        #print "min/max lambda=",min(lam),max(lam)
        #print "min/max sed_noscat=",sed_noscat.min(),sed_noscat.max()
        #print "min/max/mean f_continuum_ns=",min(f_continuum_ns),max(f_continuum_ns),mean(f_continuum_ns)
        #print "min/max sed_scat=",min(sed_scat),max(sed_scat)
        #print "min/max/mean f_continuum=",min(f_continuum),max(f_continuum),mean(f_continuum)
        
        sed_noscat -= f_continuum_ns
        sed_scat -= f_continuum
        
        print "min/max subtracted sed_noscat=",min(sed_noscat),max(sed_noscat)
        print "min/max subtracted sed_scat=",min(sed_scat),max(sed_scat)
        print "difference b/t non-scatter & scatter SED normalization (min/max/mean):"
        print min(f_continuum_ns-f_continuum),\
              max(f_continuum_ns-f_continuum),\
              mean(f_continuum_ns-f_continuum)
        
    
    x_kpc_per_pix = kpc_per_pix[0]
    y_kpc_per_pix = kpc_per_pix[1]

    pltpos = ( pos - ylim[0] - 0.5*(ylim[1]-ylim[0]) )*y_kpc_per_pix
    plthfw = hfw * y_kpc_per_pix

    xran = array([(-0.5*(xlim[1]-xlim[0]))*x_kpc_per_pix, (0.5*(xlim[1]-xlim[0]))*x_kpc_per_pix])
    yran = array([(-0.5*(ylim[1]-ylim[0]))*y_kpc_per_pix, (0.5*(ylim[1]-ylim[0]))*y_kpc_per_pix])

    ## set normalizations for axis and color scalings
    if subtract_continuum:
        sed_norm = max([ f_continuum_ns.max(), 1.05*max([sed_noscat.max(),sed_scat.max()]) ])
        sed_scatter_norm = max( [f_continuum.max(), 1.05*sed_scat.max() ])
    else:
        sed_norm = max([ 1.0e9, 1.05*max([sed_noscat.max(),sed_scat.max()]) ])
        sed_scatter_norm = max([ 1.0e9, 1.05*sed_scat.max() ])

 
    ### START HERE WITH UPDATE TO THIS FUNCTION ###

    img_vmax = min([ max_lum, max([slit_noscat.max(), img_noscat.max(),
                                   slit_scat.max(), img_noscat.max()]) ])
    img_vmin = max([ 1.0, img_vmax - dynrng_lum,
                     min([slit_noscat.min(), img_noscat.min(),
                          slit_scat.min(), img_noscat.min()]) ])

    img_scatter_vmax = min([ max_lum, max([slit_scat.max(), img_noscat.max()]) ])
    img_scatter_vmin = max([ 1.0, img_scatter_vmax - dynrng_lum,
                             min([slit_scat.min(), img_noscat.min()]) ])

    ## convert to Angstrom, find normalizations, & set rest wavelength to color scale midpoint
    tmp_vmax = min([ lam[-1]*1.0e10, max([ifu_noscat[img_noscat>=img_vmin].max(),
                                        ifu_scat[img_scat>=img_vmin].max()]) ])
    tmp_vmin = max([ lam[0]*1.0e10, min([ifu_noscat[img_noscat>=img_vmin].min(),
                                        ifu_scat[img_scat>=img_vmin].min()]) ])
    #tmp_vmax = min([ lam[llim[1]]*1.0e10, max([ifu_noscat[img_noscat>=img_vmin].max(),
    #                                    ifu_scat[img_scat>=img_vmin].max()]) ])
    #tmp_vmin = max([ lam[llim[0]]*1.0e10, min([ifu_noscat[img_noscat>=img_vmin].min(),
    #                                    ifu_scat[img_scat>=img_vmin].min()]) ])
    if tmp_vmax < lambda_rest or tmp_vmin > lambda_rest:
        print "Error: ifu wavelength range does not include rest wavelength!"
        print "range = (",tmp_vmin,":",tmp_vmax,")"
        quit()        

    f.close()

    ## set overall normalizations
    print "sed scale range = (0, ",sed_norm,")"
    print "sed scatter scale range = (0, ",sed_scatter_norm,")"

    #img_norm = colors.Normalize(vmin=img_vmin,vmax=img_vmax)
    img_norm = array([img_vmin,img_vmax])
    print "(log) luminosity scale range = ",img_norm
    print "(log) luminosity scatter scale range [NOT CURRENTLY USED]= (",img_scatter_vmin,":",img_scatter_vmax,")"

    #DEBUG -- doing this manually!!
    #lam_half_rng = max([abs(tmp_vmax-lambda_rest),abs(lambda_rest-tmp_vmin)])
    lam_half_rng = 10.0 ## Angstroms
    ifu_vmax = lambda_rest + lam_half_rng
    ifu_vmin = lambda_rest - lam_half_rng
    ifu_norm = array([ifu_vmin, ifu_vmax])
    print "ifu color scale range = (",ifu_norm[0],":",ifu_norm[1],")"

    ## clip low-lum values
    ifu_noscat[img_noscat<img_vmin] = nan
    img_noscat[img_noscat<img_vmin] = nan
    ifu_scat[img_scat<img_vmin] = nan
    img_scat[img_scat<img_vmin] = nan

    lran = ifu_norm

    ## plot it
    plot_onecam(file=shortdir, snap=snap, cam=cam, lam=lam, lran=lran, lambda_rest=lambda_rest,
                xran=xran, yran=yran,
                sed_norm=sed_norm, sed_scatter_norm=sed_scatter_norm, img_norm=img_norm, 
                pltpos=pltpos, plthfw=plthfw, 
                sed_noscat=sed_noscat, sed_scat=sed_scat,
                slit_noscat=slit_noscat, slit_scat=slit_scat, img_noscat=img_noscat,
                img_scat=img_scat, ifu_noscat=ifu_noscat, ifu_scat=ifu_scat)

    plot_onecam_spec_only(file=shortdir, snap=snap, cam=cam, 
                          lam=lam, lran=lran, lambda_rest=lambda_rest,
                          xran=xran, yran=yran,
                          sed_norm=sed_norm, sed_scatter_norm=sed_scatter_norm, 
                          img_norm=img_norm, pltpos=pltpos, plthfw=plthfw, 
                          sed_noscat=sed_noscat, sed_scat=sed_scat,
                          slit_noscat=slit_noscat, slit_scat=slit_scat)



def process_onecam(fileobj, snap, cam, lam, llim, lambda_rest, **kwargs):

    if kwargs is not None:
        print kwargs.items()
        for name, value in kwargs.items():
            if name is 'max_lum':
                max_lum = value
            elif name is 'dynrng_lum':
                dynrng_lum = value
            else:
                print "value for slit_pos (",int(value),")is outside of plot range:",ylim
                print "   ...using default value:",pos
                    

    im_noscat = float128(fileobj['camera%d-nonscatter'%cam].data)
    hdr_noscat = fileobj['camera%d-nonscatter'%cam].header
    im_scat = float128(fileobj['camera%d'%cam].data)
    hdr_scat = fileobj['camera%d'%cam].header
    hdr_campar = fileobj['camera%d-parameters'%cam].header

    lam = float128(lam)

    #print im_noscat.sum(axis=2).shape
    #print im_noscat.sum(axis=2).sum(axis=1).shape
    #print im_noscat.sum(axis=2).sum(axis=1).sum(axis=0).shape

    print "min/max/median lam=",lam.min(),lam.max(),median(lam)
    cen_offset = lam[lam<lambda_rest*1.0e-10].size - \
        lam[lam>lambda_rest*1.0e-10].size
    if mod(lam.size-cen_offset,2) != 0: cen_offset += 1
    # print "WARNING! Manually setting cen_offset to 0."
    # cen_offset = 0
    extra = 0
    if cen_offset > 0:
        print "lambda has ",cen_offset," more blue than red bins."
        lam = lam[cen_offset:]
        im_noscat = im_noscat[cen_offset:, :, :]
        im_scat = im_scat[cen_offset:, :, :]
    elif cen_offset < 0:
        print "lambda has ",cen_offset," more red than blue bins."        
        lam = lam[:cen_offset]
        im_noscat = im_noscat[:cen_offset, :, :]
        im_scat = im_scat[:cen_offset, :, :]
    if extra > 0:
        print "clipping ",extra," extra bins from ends of wavelength range."
        lam = lam[extra:-extra]
        im_noscat = im_noscat[extra:-extra,:,:]
        im_scat = im_scat[extra:-extra,:,:]
    llim = array([0,lam.size])

    print "num of blue lam bins: ",lam[lam<lambda_rest*1.0e-10].size
    print "num of red lam bins: ",lam[lam>lambda_rest*1.0e-10].size
    print "im_noscat, im_scat shapes: ", im_noscat.shape, im_scat.shape

    cameradist = hdr_campar['cameradist'] ## in kpc
    print "\ncameradist=",cameradist
    #for old code:
    #sr_per_pix = hdr_scat['PIXEL_SR']
    #for new code (sunrise 4):
    xpixels = hdr_campar['XSIZE']
    ypixels = hdr_campar['YSIZE']
    sr_per_pix = hdr_campar['solid_angle']/(xpixels*ypixels)
    x_kpc_per_pix = hdr_scat['CD2_2']
    y_kpc_per_pix = hdr_scat['CD1_1']    
    #print "x_kpc_per_pix, y_kpc_per_pix = ",x_kpc_per_pix,y_kpc_per_pix
    print "sr_per_pix = ",sr_per_pix
    if x_kpc_per_pix != hdr_noscat['CD2_2'] or \
            y_kpc_per_pix != hdr_noscat['CD1_1']:
        print "Error! Pixel scale mismatch between scatter and nonscatter cameras."
        print "nonscatter kpc_per_pix = ",hdr_noscat['CD2_2'],hdr_noscat['CD1_1']
        quit()

    ## code returns surface brightness in W/m/m^2/sr (cameradist is in kpc)
    ## here we convert image unit to lambda*L_lambda in erg/s
    codeSB_to_ergs = sr_per_pix * 1.0e7 * 4*pi * (cameradist*3.09e19)**2
    print "codeSB_to_ergs=",codeSB_to_ergs

    deltalam = float128(zeros(lam.size))
    for l in arange(lam.size-2)+1: deltalam[l] = 0.5*(lam[l+1]-lam[l-1])
    deltalam[0] = deltalam[1]
    deltalam[-1] = deltalam[-2]

    for l in range(lam.size):        
        #im_noscat[l,:,:] *= lam[l]
        #im_scat[l,:,:] *= lam[l]
        im_noscat[l,:,:] *= deltalam[l]
        im_scat[l,:,:] *= deltalam[l]
    im_noscat_ergs = im_noscat * codeSB_to_ergs
    im_scat_ergs = im_scat * codeSB_to_ergs

    print "min/max im_noscat_ergs: ",\
        im_noscat_ergs[im_noscat_ergs>0].min(),\
        im_noscat_ergs[im_noscat_ergs>0].max()
    print "min/max im_scat_ergs: ",\
        im_scat_ergs[im_scat_ergs>0].min(), im_scat_ergs[im_scat_ergs>0].max()

    ## note: imshow plots the 1st dimension of a 2-element array as y.
    #xpixels = im_noscat_ergs.sum(axis=0).sum(axis=0).shape[0]
    #ypixels = im_noscat_ergs.sum(axis=0).sum(axis=1).shape[0]

    if xpixels != ypixels:
        print "Warning! Camera array is not square."
        print "Plot aspect ratios may not be to scale."

    if xpixels != im_scat_ergs.sum(axis=0).sum(axis=0).shape[0] or \
            ypixels != im_scat_ergs.sum(axis=0).sum(axis=1).shape[0]:
        print "Error! Pixel number mismatch between scatter and nonscatter cameras."
        print "nonscatter: ",xpixels,ypixels
        print "scatter: ",im_scat_ergs.sum(axis=0).sum(axis=0).shape[0],\
            im_scat_ergs.sum(axis=0).sum(axis=1).shape[0]
        quit()

    xlim = array([0,xpixels])
    ylim = array([0,ypixels])

    ## hardcoding these for now;
    ## not worth defining as kwargs unless they're in more physical units.
    pos = int(ypixels/2.0)
    hfw = max([2*int(1.0/y_kpc_per_pix),1])

    ## convert luminosity range from Lsun to erg/s:    
    max_lum = 10.0**(max_lum) * 3.83e33
    dynrng_lum = 10.0**(dynrng_lum)
    min_lum = max_lum/dynrng_lum 
    print "Clipping luminosity range:"
    print "max_lum = ",max_lum," erg/s"
    print "min_lum = ",min_lum," erg/s"
    ## **** lum from each pixel. 

    im_tmp = im_noscat_ergs[llim[0]:llim[1],ylim[0]:ylim[1],xlim[0]:xlim[1]]
    arr_tmp = copy(im_tmp).sum(axis=0)
    print "Clipping ",arr_tmp[arr_tmp<min_lum].size,\
        " elements of arr_tmp (for sed_noscat)."
    for i in range(arr_tmp[:,0].size):
        for j in range(arr_tmp[0,:].size):
            if arr_tmp[i,j]<min_lum: im_tmp[:,i,j] = 0.0
    #print "Checking..."
    #print "Shape of arr_tmp and im_tmp:",arr_tmp.shape,im_tmp.shape
    #print "zero elements of im_tmp:",im_tmp[im_tmp==0.0].size
    if arr_tmp[arr_tmp>max_lum].size > 0:
        print "WARNING! ", arr_tmp[arr_tmp>max_lum].size, \
            " elements of arr_tmp (for sed_noscat) exceed max_lum."
    sed_noscat = copy(im_tmp).sum(axis=2).sum(axis=1)
    #sed_noscat = (im_noscat_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]]).sum(axis=2).sum(axis=1)

    slit_tmp = im_tmp[:,pos-hfw:pos+hfw,:]
    slit_noscat = copy(slit_tmp).sum(axis=1)
    del im_tmp, arr_tmp, slit_tmp
    #slit_noscat = sum(im_noscat_ergs[llim[0]:llim[1],pos-hfw:pos+hfw,
    #                                 xlim[0]:xlim[1]],axis=1)
    print "Setting ",slit_noscat[slit_noscat==0.0].size, \
        " elements of slit_noscat to 1.0e-30."
    slit_noscat[slit_noscat==0.0] = 1.0e-30
    slit_noscat = log10(slit_noscat)

    img_noscat = sum(im_noscat_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]],axis=0)
    print img_noscat[img_noscat==0.0].size, \
        " elements of img_noscat are initially 0."
    print "initial min/max/mean/median of (nonzero) img_noscat: ",\
        img_noscat[img_noscat>0.0].min(),img_noscat.max(), \
        img_noscat[img_noscat>0.0].mean(),median(img_noscat[img_noscat>0.0])
    print "Clipping ",img_noscat[img_noscat<min_lum].size, \
        " low-lum elements of img_noscat."
    img_noscat[img_noscat<min_lum] = 0.0
    if img_noscat[img_noscat>max_lum].size > 0:
        print "WARNING! ",img_noscat[img_noscat>max_lum].size,\
            " elements of img_noscat exceed max_lum=",max_lum
    print "Setting ",img_noscat[img_noscat==0.0].size, \
        " elements of img_noscat to 1.0e-30."
    print "final min/max/mean/median of (nonzero) img_noscat: ",\
        img_noscat[img_noscat>0.0].min(),img_noscat.max(), \
        img_noscat[img_noscat>0.0].mean(),median(img_noscat[img_noscat>0.0])
    img_noscat[img_noscat==0.0] = 1.0e-30
    img_noscat = log10(img_noscat)


    srt_lum_noscat = sort(reshape(img_noscat,(img_noscat.size)))
    lowlum_noscat = srt_lum_noscat[floor(0.2*srt_lum_noscat.size)]
    print 0.2*srt_lum_noscat.size
    print "lowlum_noscat = ",lowlum_noscat
    print img_noscat[img_noscat<lowlum_noscat].size
    print "min/max/mean img_noscat = ",img_noscat.min(),img_noscat.max(),img_noscat.mean()

    lowlum_sed_noscat = zeros((sed_noscat.size))
    for i in range(xpixels):
        for j in range(ypixels):
            if img_noscat[j,i] < lowlum_noscat:
                lowlum_sed_noscat += im_noscat_ergs[:,j,i]
    print lowlum_sed_noscat.shape
    Lblue_lowlum = lowlum_sed_noscat[lam<1.0e-10*lambda_rest].sum()
    Lred_lowlum = lowlum_sed_noscat[lam>1.0e-10*lambda_rest].sum()
    Ltot_lowlum = lowlum_sed_noscat.sum()
    #print "Lblue_lowlum, Lred_lowlum, Ltot_lowlum = ",Lblue_lowlum, Lred_lowlum, Ltot_lowlum


    #dlam = lam - 1.0e-10*lambda_rest
    #print "min/max/mean/med dlam=",dlam.min(),dlam.max(),dlam.mean(),median(dlam)
    ##Calculate IFU spectrum in Angstrom
    #allpix_noscat = im_noscat.sum(axis=0)
    allpix_noscat = im_noscat_ergs.sum(axis=0)
    print "initial min/max/mean/median of (nonzero) allpix_noscat: ",\
        allpix_noscat[allpix_noscat>0.0].min(),allpix_noscat.max(), \
        allpix_noscat[allpix_noscat>0.0].mean(),\
        median(allpix_noscat[allpix_noscat>0.0])
    print allpix_noscat[allpix_noscat==0.0].size, \
        " elements of allpix_noscat are initially 0."
    print "Clipping ",allpix_noscat[allpix_noscat<min_lum].size, \
        " low-lum elements of allpix_noscat."
    allpix_noscat[allpix_noscat<min_lum] = 0.0
    if allpix_noscat[allpix_noscat>max_lum].size > 0:
        print "WARNING! ",allpix_noscat[allpix_noscat>max_lum].size,\
            " elements of allpix_noscat exceed max_lum=",max_lum
    print "final min/max/mean/median of (nonzero) allpix_noscat: ",\
        allpix_noscat[allpix_noscat>0.0].min(),allpix_noscat.max(), \
        allpix_noscat[allpix_noscat>0.0].mean(),\
        median(allpix_noscat[allpix_noscat>0.0])


    ################################################################
    ### TESTING ONLY: continuum subtraction for IFU spectra
    test_continuum_ns, test_vfinal_sigma_ns = fit_continuum(lam,sed_noscat, 5, 1.0, 0.02e-7)

    test_im_noscat_ergs_nocont = copy(im_noscat_ergs)
    cont_diff = zeros((lam.size,ypixels,xpixels))
    tot_allpix_noscat = allpix_noscat.sum(axis=1).sum(axis=0)
    for i in range(xpixels):
        #print i
        if mod(i,20) == 0: print i
        for j in range(ypixels):            
            if allpix_noscat[j,i] > 0.0:
                #f_cont, vf_sig = fit_continuum(lam,im_noscat_ergs[:,j,i], 5, 1.0, 0.02e-7)
                #test_im_noscat_ergs_nocont[:,j,i] = im_noscat_ergs[:,j,i] - f_cont
                test_im_noscat_ergs_nocont[:,j,i] = \
                    im_noscat_ergs[:,j,i] - test_continuum_ns*allpix_noscat[j,i]/tot_allpix_noscat 
                #cont_diff[:,j,i] = test_continuum_ns*allpix_noscat[j,i]/tot_allpix_noscat - f_cont
    print "min/max test_im_noscat_ergs_nocont=",\
        test_im_noscat_ergs_nocont.min(),test_im_noscat_ergs_nocont.max()
    print "setting ",test_im_noscat_ergs_nocont[test_im_noscat_ergs_nocont<0.0].size,\
        " negative elementso of test_im_noscat_ergs_nocont to 0.0."
    print "   (",100.0*test_im_noscat_ergs_nocont[test_im_noscat_ergs_nocont<0.0].size/test_im_noscat_ergs_nocont.size,"% of all elements.)" 
    test_im_noscat_ergs_nocont[test_im_noscat_ergs_nocont<0.0] = 0.0
    print "min/max test_im_noscat_ergs_nocont=",\
        test_im_noscat_ergs_nocont.min(),test_im_noscat_ergs_nocont.max()
    #print "min/max/mean continuum normalization = ",\
    #    allpix_noscat[allpix_noscat>0].min()*test_continuum_ns/tot_allpix_noscat,\
    #    allpix_noscat.max()*test_continuum_ns/tot_allpix_noscat,\
    #    allpix_noscat[allpix_noscat>0].mean()*test_continuum_ns/tot_allpix_noscat

    #print "min_cont_diff = ",abs(cont_diff[cont_diff!=0.0]).min()
    #print "max_cont_diff = ",abs(cont_diff).max()
    #print "mean_cont_diff = ",abs(cont_diff[cont_diff!=0.0]).mean()

    ################################################################

    ifu_noscat = tensordot(im_noscat_ergs,lam,axes=(0,0))
    ##ifu_noscat = tensordot(im_noscat,lam,axes=(0,0))
    ifu_noscat[allpix_noscat!=0.0] = ( 1.0e10 * ifu_noscat[allpix_noscat!=0.0] /
                                            allpix_noscat[allpix_noscat!=0.0] ) 
    ifu_noscat[allpix_noscat==0.0] = lambda_rest
    print "min/max/mean/med ifu_noscat=",ifu_noscat.min(),ifu_noscat.max(),\
        ifu_noscat.mean(),median(ifu_noscat)



    ### TESTING IFU array:
    ## next plan: figure out what the difference is between averaging over im_noscat and im_noscat_ergs (the nonzero values?)
    ##also, make some fake data with perfect gaussian and test that. --> NO blue/red difference for fake data, 
    ##but IFU data points are still a tiny bit shifted. small enough to attribute to rounding error? probably. 
    ##so need to take a closer look at the real data.

    print "total Lred =  ",(im_noscat[lam*1.0e10>lambda_rest,:,:]).sum(axis=2).sum(axis=1).sum(axis=0)
    print "total Lblue = ",(im_noscat[lam*1.0e10<lambda_rest,:,:]).sum(axis=2).sum(axis=1).sum(axis=0)


    #bluered_diff=float128(zeros((xpixels,ypixels)))
    #bluered_diff2=float128(zeros((xpixels,ypixels)))
    #for k in range(xpixels):
    #    for l in range(ypixels):
    #        bluered_diff[l,k]=(sum(im_noscat_ergs[lam*1.0e10<lambda_rest,l,k])-
    #                           sum(im_noscat_ergs[lam*1.0e10>lambda_rest,l,k]))
    #        bluered_diff2[l,k]=(sum(im_noscat[lam*1.0e10<lambda_rest,l,k])-
    #                            sum(im_noscat[lam*1.0e10>lambda_rest,l,k]))
    #print "min/max/mean bluered_diff = ",bluered_diff.min(),bluered_diff.max(),bluered_diff.mean()
    #print "min/max/mean bluered_diff2 = ",bluered_diff2.min(),bluered_diff2.max(),bluered_diff2.mean()
    #print bluered_diff[bluered_diff>0.0].shape," elements of bluered_diff are preferentially blue"
    #print bluered_diff2[bluered_diff2>0.0].shape," elements of bluered_diff2 are preferentially blue"


    Lblue = sum(sed_noscat[lam*1.0e10<lambda_rest])
    Lred = sum(sed_noscat[lam*1.0e10>lambda_rest])
    Ltot = sum(sed_noscat)
    print "Lblue=",Lblue
    print "Lred=",Lred
    print "Ltot=",Lblue+Lred, Ltot

    im_tmp = im_scat_ergs[llim[0]:llim[1],ylim[0]:ylim[1],xlim[0]:xlim[1]]
    arr_tmp = copy(im_tmp).sum(axis=0)
    print "Clipping ",arr_tmp[arr_tmp<min_lum].size,\
        " elements of arr_tmp (for sed_scat)."
    for i in range(arr_tmp[:,0].size):
        for j in range(arr_tmp[0,:].size):
            if arr_tmp[i,j]<min_lum: im_tmp[:,i,j] = 0.0
    #print "Checking..."
    #print "Shape of arr_tmp and im_tmp:",arr_tmp.shape,im_tmp.shape
    #print "zero elements of im_tmp:",im_tmp[im_tmp==0.0].size
    if arr_tmp[arr_tmp>max_lum].size > 0:
        print "WARNING! ", arr_tmp[arr_tmp>max_lum].size, \
            " elements of arr_tmp (for sed_scat) exceed max_lum."
    sed_scat = copy(im_tmp).sum(axis=2).sum(axis=1)
    #sed_scat = (im_scat_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]]).sum(axis=2).sum(axis=1)

    slit_tmp = im_tmp[:,pos-hfw:pos+hfw,:]
    slit_scat = copy(slit_tmp).sum(axis=1)
    del im_tmp, arr_tmp, slit_tmp
    #slit_scat = sum(im_scat_ergs[llim[0]:llim[1],pos-hfw:pos+hfw,
    #                             xlim[0]:xlim[1]],axis=1)
    print slit_scat[slit_scat==0.0].size, \
        " elements of slit_scat are initially 0."
    print "Setting ",slit_scat[slit_scat==0.0].size, \
        " elements of slit_scat to 1.0e-30."
    slit_scat[slit_scat==0.0] = 1.0e-30
    slit_scat = log10(slit_scat)

    img_scat = sum(im_scat_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]],axis=0)
    print img_scat[img_scat==0.0].size, \
        " elements of img_scat are initially 0."
    print "initial min/max/mean/median of (nonzero) img_scat: ",\
        img_scat[img_scat>0.0].min(),img_scat.max(), \
        img_scat[img_scat>0.0].mean(),median(img_scat[img_scat>0.0])
    print "Clipping ",img_scat[img_scat<min_lum].size, \
        " low-lum elements of img_scat."
    img_scat[img_scat<min_lum] = 0.0
    if img_scat[img_scat>max_lum].size > 0:
        print "WARNING! ",img_scat[img_scat>max_lum].size,\
            " elements of img_scat exceed max_lum=",max_lum
    print "Setting ",img_scat[img_scat==0.0].size, \
        " elements of img_scat to 1.0e-30."
    img_scat[img_scat==0.0] = 1.0e-30
    img_scat = log10(img_scat)

    
    ##Calculate IFU spectrum in Angstrom
    #allpix_scat = im_scat.sum(axis=0)
    allpix_scat = im_scat_ergs.sum(axis=0)
    print allpix_scat[allpix_scat==0.0].size, \
        " elements of allpix_scat are initially 0."
    print "initial min/max/mean/med allpix_scat=",\
        allpix_scat.min(),allpix_scat.max(),\
        allpix_scat.mean(),median(allpix_scat)
    print "Clipping ",allpix_scat[allpix_scat<min_lum].size, \
        " low-lum elements of allpix_scat."
    allpix_scat[allpix_scat<min_lum] = 0.0
    if allpix_scat[allpix_scat>max_lum].size > 0:
        print "WARNING! ",allpix_scat[allpix_scat>max_lum].size,\
            " elements of allpix_scat exceed max_lum=",max_lum
    print "final min/max/mean/med allpix_scat=",\
        allpix_scat.min(),allpix_scat.max(),\
        allpix_scat.mean(),median(allpix_scat)
    ifu_scat = tensordot(im_scat_ergs,lam,axes=(0,0))
    ifu_scat[allpix_scat!=0.0] = ( 1.0e10 * ifu_scat[allpix_scat!=0.0] /
                                   allpix_scat[allpix_scat!=0.0] ) 
    ifu_scat[allpix_scat==0.0] = lambda_rest
    print "min/max/mean/med ifu_scat=",\
        ifu_scat.min(),ifu_scat.max(),ifu_scat.mean(),median(ifu_scat)


    print "total Lred =  ",(im_noscat[lam*1.0e10>lambda_rest,:,:]).sum(axis=2).sum(axis=1).sum(axis=0)
    print "total Lblue = ",(im_noscat[lam*1.0e10<lambda_rest,:,:]).sum(axis=2).sum(axis=1).sum(axis=0)


    #print "\nimg_noscat[0,0] = ",img_noscat[0,0]
    #print "ifu_noscat[0,0] = ",ifu_noscat[0,0]
    #print "Lblue [0,0] = ",sum(im_noscat[lam*1.0e10<lambda_rest,0,0]),\
    #    " Lred [0,0] = ",sum(im_noscat[lam*1.0e10>lambda_rest,0,0])
 
    return pos, hfw, xlim, ylim, [x_kpc_per_pix,y_kpc_per_pix], lam, \
        sed_noscat, sed_scat, slit_noscat, slit_scat, \
        img_noscat, img_scat, ifu_noscat, ifu_scat, \
        lowlum_sed_noscat
    


def plot_onecam(file, snap, cam, lam, lran, lambda_rest, xran, yran,
                sed_norm, sed_scatter_norm, img_norm, pltpos, plthfw, 
                sed_noscat, sed_scat,
                slit_noscat, slit_scat, img_noscat, img_scat, ifu_noscat, ifu_scat):


    print "checking input to plot onecam:"
    print "snap=",snap,", cam=",cam
    print "lam: ",lam.shape
    print " lambda rest=",lambda_rest
    print "sed_norm:",sed_norm
    print "sed_scatter_norm:",sed_scatter_norm
    print "img_norm:",img_norm
    print "sed_noscat:",sed_noscat.shape,sed_noscat.dtype
    print "sed_scat:",sed_scat.shape,sed_scat.dtype
    print "slit_noscat:",slit_noscat.shape,slit_noscat.dtype
    print "slit_scat:",slit_scat.shape,slit_scat.dtype
    print "img_noscat:",img_noscat.shape,img_noscat.dtype
    print "img_scat:",img_scat.shape,img_scat.dtype

    sed_noscat = float64(sed_noscat)
    sed_scat = float64(sed_scat)
    slit_noscat = float64(slit_noscat)
    slit_scat = float64(slit_scat)
    img_noscat = float64(img_noscat)
    img_scat = float64(img_scat)
    ifu_noscat = float64(ifu_noscat)
    ifu_scat = float64(ifu_scat)

    ## get intrinsic aspect ratios for axes
    slit_aspect = (lran[1]-lran[0])/float(yran[1]-yran[0])
    img_aspect = float(yran[1]-yran[0])/(xran[1]-xran[0])

    #plt.ioff()
    plt.clf()
    fig = plt.figure(figsize=(7,8))
    matplotlib.rcParams.update({'font.size': 12})

    ### 1-D SPECTRA ###
    
    gs0 = gridspec.GridSpec(1,2)
    gs0.update(left=0.12,right=0.79, bottom=0.82, top=0.92, wspace=0.15, hspace=0.3)

    ax01 = plt.subplot(gs0[0])
    plt.plot(array([lambda_rest,lambda_rest]),array([0,sed_norm]),'black',linestyle=':')
    plt.plot(lam*1.0e10, sed_noscat)
    plt.xlim(lran[0]-10,lran[1]+10)
    #plt.ylim(0,sed_norm*1.05)
    plt.ylim(1e2,sed_norm*1.05)
    plt.xlabel(r'$\lambda$ [$\AA$]',fontsize=11)
    #plt.ylabel(r'$\lambda L_{\lambda}$ [erg s$^{-1}$]')
    #plt.ylabel(r'$\lambda L_{\lambda}$ [log L$_{\odot}$]')
    plt.ylabel(r'$\lambda L_{\lambda}$ [arb. units]')
    plt.title('no scatter',fontsize=12)
    #ax01.set_xticks([4840,4860,4880])
    ax01.xaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in ax01.get_xticklabels(): item.set_fontsize(10)
    #ax01.yaxis.set_major_locator(ticker.MaxNLocator(4))
    #ax01.ticklabel_format(axis='y',style='sci',scilimits=(-3,4))
    ax01.set_yticklabels([])

    ax02 = plt.subplot(gs0[1])
    plt.plot(array([lambda_rest,lambda_rest]),array([0,sed_scatter_norm]),'black',linestyle=':')
    plt.plot(lam*1.0e10, sed_scat)
    plt.xlim(lran[0]-10,lran[1]+10)
    plt.ylim(1e2,sed_scatter_norm*1.05)
    #plt.ylim(0,sed_scatter_norm*1.05)
    #plt.ylim(0,sed_norm)
    plt.xlabel(r'$\lambda$ [$\AA$]',fontsize=11)
    plt.title('scatter',fontsize=12)
    ax02.xaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in ax02.get_xticklabels(): item.set_fontsize(10)
    #ax02.yaxis.set_major_locator(ticker.MaxNLocator(4))
    #ax02.ticklabel_format(axis='y',style='sci',scilimits=(-3,4))
    ax02.set_yticklabels([])

     ### SLIT SPECTRA ###
    xtrim = 1.1

    gs1 = gridspec.GridSpec(1,2)
    gs1.update(left=0.12,right=0.79, bottom=0.6, top=0.755, wspace=0.15, hspace=0.3)

    ax11 = plt.subplot(gs1[0])
    plt.imshow(slit_noscat, 
               extent=array([[xran[0]+xtrim,xran[1]-xtrim],lran]).flatten(),
               cmap='hot',norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),origin='lower',
               interpolation='nearest',aspect='auto')
    #interpolation='nearest',aspect=0.5/slit_aspect)
    plt.ylabel(r'$\lambda$ [$\AA$]')
    ax11.yaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in (ax11.get_yticklabels()): item.set_fontsize(10)
    ax11.set_xticklabels([])

    ax12 = plt.subplot(gs1[1])
    plt.imshow(slit_scat,
               extent=array([[xran[0]+xtrim,xran[1]-xtrim],lran]).flatten(),
               origin='lower',cmap='hot',
               norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),
               interpolation='nearest',aspect='auto')
    ax12.set_xticklabels([])
    ax12.set_yticklabels([])

    ### IMAGES ###

    gs2 = gridspec.GridSpec(1,2)
    gs2.update(left=0.12,right=0.79, bottom=0.33, top=0.59, wspace=0.15, hspace=0.3)
    
    ax21 = plt.subplot(gs2[0])
    plt.plot([xran[0]+xtrim,xran[1]-xtrim],
             array([pltpos-plthfw,pltpos-plthfw]),'w--',linewidth=1.0)
    plt.plot([xran[0]+xtrim,xran[1]-xtrim],
             array([pltpos+plthfw,pltpos+plthfw]),'w--',linewidth=1.0)
    plt.imshow(img_noscat,interpolation='nearest',
               norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),
               cmap='hot', extent=array([[xran[0]+xtrim,xran[1]-xtrim],
                                         [yran[0]+xtrim,yran[1]-xtrim]]).flatten(),
               origin='lower',aspect='auto')
    ax21.patch.set_facecolor('black')
    #ax21.set_xticklabels([])
    plt.xlabel("x [kpc]")
    plt.ylabel("y [kpc]")

    ax22 = plt.subplot(gs2[1])
    plt.plot([xran[0]+xtrim,xran[1]-xtrim],
             array([pltpos-plthfw,pltpos-plthfw]),'w--',linewidth=1.0)
    plt.plot([xran[0]+xtrim,xran[1]-xtrim],
             array([pltpos+plthfw,pltpos+plthfw]),'w--',linewidth=1.0)
    plt.imshow(img_scat,interpolation='nearest',
               norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),
               cmap='hot',extent=array([[xran[0]+xtrim,xran[1]-xtrim],
                                        [yran[0]+1,yran[1]-1]]).flatten(),
               origin='lower',aspect='auto')
    ax22.patch.set_facecolor('black')
    plt.xlabel("x [kpc]")
    #ax22.set_xticklabels([])
    ax22.set_yticklabels([])

    ## colorbar for luminosity scale
    cbaxes_lum = fig.add_axes([0.82,0.34,0.03,0.41])
    cb = plt.colorbar(cax = cbaxes_lum)
    #cb.set_label(r'$\lambda L_{\lambda}$ [erg s$^{-1}$]')
    cb.set_label(r'$\lambda L_{\lambda}$ [log L$_{\odot}$]')

    ### IFU SPECTRA ###

    #gs3 = gridspec.GridSpec(1,2)
    #gs3.update(left=0.12,right=0.79, bottom=0.06, top=0.32, wspace=0.15, hspace=0.3)

    #ax31 = plt.subplot(gs3[0])
    ##plt.imshow(ifu_noscat,interpolation='nearest',cmap='seismic',
    ##           norm=colors.Normalize(vmin=lran[0]+0.1*(lran[1]-lran[0]),vmax=lran[1]-0.1*(lran[1]-lran[0])),
    ##           extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    #plt.imshow(ifu_noscat,interpolation='nearest',cmap='seismic',
    #           extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    #plt.xlabel("x [kpc]")
    #plt.ylabel("y [kpc]")

    #ax32 = plt.subplot(gs3[1])
    #plt.imshow(ifu_scat,interpolation='nearest',cmap='seismic',
    #           extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    #ax32.set_yticklabels([])
    #plt.xlabel("x [kpc]")

    ## colorbar for wavelength scale
    #cbaxes_ifu = fig.add_axes([0.82,0.07,0.03,0.24])
    #cb2 = plt.colorbar(cax = cbaxes_ifu)
    #cb2.set_label(r'$\lambda$ [$\AA$]')

    #plt.suptitle('%s: camera %d'%(file,cam))

    #DefaultSize = fig.get_size_inches()
    #print "default size in inches = ",DefaultSize
    #fig.set_size_inches(8.0,10.0)

    #filebase = ((maindir+file).split('mcrx'))[0]
    #if filebase != '': filebase = filebase+'/'
 
    #filename=filebase+"%d-c%d-hbetaspec.eps"%(snap,cam)
    filename="%d-c%d-hbetaspec.eps"%(snap,cam)
    #if subtract_continuum: filename=filebase+"%d-c%d-hbetaspec.eps"%(snap,cam)
    #else: filename=filebase+"%d-c%d-hbetaspec_withcontinuum.eps"%(snap,cam)
    plt.savefig(filename)
    plt.clf()
    plt.cla()
    plt.close('all')

def plot_onecam_spec_only(file, snap, cam, lam, lran, lambda_rest, xran, yran,
                          sed_norm, sed_scatter_norm, img_norm, pltpos, plthfw, 
                          sed_noscat, sed_scat,slit_noscat, slit_scat):


    print "checking input to plot onecam:"
    print "snap=",snap,", cam=",cam
    print "lam: ",lam.shape
    print " lambda rest=",lambda_rest
    print "sed_norm:",sed_norm
    print "img_norm:",img_norm
    print "sed_scatter_norm:",sed_scatter_norm
    print "sed_noscat:",sed_noscat.shape,sed_noscat.dtype
    print "sed_scat:",sed_scat.shape, sed_scat.dtype
    print "slit_noscat:",slit_noscat.shape,slit_noscat.dtype
    print "slit_scat:",slit_scat.shape,slit_scat.dtype

    sed_noscat = float64(sed_noscat)
    sed_scat = float64(sed_scat)
    slit_noscat = float64(slit_noscat)
    slit_scat = float64(slit_scat)

    ## get intrinsic aspect ratios for axes
    slit_aspect = (lran[1]-lran[0])/float(yran[1]-yran[0])

   #plt.ioff()
    plt.clf()
    fig = plt.figure(figsize=(7,4))
    matplotlib.rcParams.update({'font.size': 12})

    ### 1-D SPECTRA ###
    gs0 = gridspec.GridSpec(1,2)
    gs0.update(left=0.12,right=0.79, bottom=0.6, top=0.9, wspace=0.15, hspace=0.3)

    ax01 = plt.subplot(gs0[0])
    plt.plot(array([lambda_rest,lambda_rest]),array([0,sed_norm]),'black',linestyle=':')
    plt.plot(lam*1.0e10, sed_noscat)
    plt.xlim(lran[0]-10,lran[1]+10)
    #plt.ylim(0,sed_norm*1.05)
    plt.ylim(1e2,sed_norm*1.05)
    plt.xlabel(r'$\lambda$ [$\AA$]',fontsize=11)
    #plt.ylabel(r'$\lambda L_{\lambda}$ [erg s$^{-1}$]')
    #plt.ylabel(r'$\lambda L_{\lambda}$ [log L$_{\odot}$]')
    plt.ylabel(r'$\lambda L_{\lambda}$ [arb. units]',fontsize=10)
    plt.title('no scatter',fontsize=12)
    #ax01.set_xticks([4840,4860,4880])
    ax01.xaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in ax01.get_xticklabels(): item.set_fontsize(10)
    #ax01.yaxis.set_major_locator(ticker.MaxNLocator(4))
    #ax01.ticklabel_format(axis='y',style='sci',scilimits=(-3,4))
    ax01.set_yticklabels([])

    ax02 = plt.subplot(gs0[1])
    plt.plot(array([lambda_rest,lambda_rest]),array([0,sed_scatter_norm]),'black',linestyle=':')
    plt.plot(lam*1.0e10, sed_scat)
    plt.xlim(lran[0]-10,lran[1]+10)
    plt.ylim(1e2,sed_scatter_norm*1.05)
    #plt.ylim(0,sed_scatter_norm*1.05)
    #plt.ylim(0,sed_norm)
    plt.xlabel(r'$\lambda$ [$\AA$]',fontsize=11)
    plt.title('scatter',fontsize=12)
    ax02.xaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in ax02.get_xticklabels(): item.set_fontsize(10)
    #ax02.yaxis.set_major_locator(ticker.MaxNLocator(4))
    #ax02.ticklabel_format(axis='y',style='sci',scilimits=(-3,4))
    ax02.set_yticklabels([])

     ### SLIT SPECTRA ###
    xtrim = 1.1

    gs1 = gridspec.GridSpec(1,2)
    gs1.update(left=0.12,right=0.79, bottom=0.15, top=0.45, wspace=0.15, hspace=0.3)

    ax11 = plt.subplot(gs1[0])
    plt.imshow(transpose(slit_noscat), 
               extent=array([[lran[0]-10,lran[1]+10],
                             [xran[0]+xtrim,xran[1]-xtrim]]).flatten(),
               cmap='hot',norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),origin='lower',
               interpolation='nearest',aspect='auto')
    #interpolation='nearest',aspect=0.5/slit_aspect)
    plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.ylabel("position [kpc]",fontsize=10)
    ax11.yaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in (ax11.get_yticklabels()): item.set_fontsize(10)
    ax11.xaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in (ax11.get_xticklabels()): item.set_fontsize(10)
    #ax11.set_xticklabels([])

    ax12 = plt.subplot(gs1[1])
    plt.imshow(transpose(slit_scat),
               extent=array([[lran[0]-10,lran[1]+10],
                             [xran[0]+xtrim,xran[1]-xtrim]]).flatten(),
               origin='lower',cmap='hot',
               norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),
               interpolation='nearest',aspect='auto')
    plt.xlabel(r'$\lambda$ [$\AA$]')
    #ax12.set_xticklabels([])
    ax12.xaxis.set_major_locator(ticker.MaxNLocator(4))
    for item in (ax12.get_xticklabels()): item.set_fontsize(10)
    ax12.set_yticklabels([])

    ## colorbar for luminosity scale
    cbaxes_lum = fig.add_axes([0.82,0.15,0.03,0.3])
    cbaxes_lum.yaxis.set_major_locator(ticker.MaxNLocator(4))
    cb = plt.colorbar(cax = cbaxes_lum)
    #cb.set_label(r'$\lambda L_{\lambda}$ [erg s$^{-1}$]')
    cb.set_label(r'$\lambda L_{\lambda}$ [log L$_{\odot}$]',fontsize=10)


    filename="%d-c%d-hbeta-speconly.eps"%(snap,cam)
    #if subtract_continuum: filename=filebase+"%d-c%d-hbetaspec.eps"%(snap,cam)
    #else: filename=filebase+"%d-c%d-hbetaspec_withcontinuum.eps"%(snap,cam)
    plt.savefig(filename)
    #plt.clf()
    #plt.cla()
    #plt.close('all')

def multicam(file,snap,cam_arr=arange(7),max_lum=7.5,dynrng_lum=4.5,
             subtract_continuum=True,**kwargs):

    ## hbeta rest wavelength
    #lambda_rest = 4861.36
    lambda_rest = 4861.0
    #lambda_rest = 4860.0743711

    maindir=''
    nonlr_dir=''
    wholesed_dir=''
    if kwargs is not None:
        #print kwargs.items()
        for name, value in kwargs.items():
            if name=='maindir' and value != '':
                maindir=value+'/'
            elif name=='nonlr_dir' and value != '':
                nonlr_dir=value+'/'
            elif name=='wholesed_dir' and value != '':
                wholesed_dir=value+'/'
            else: print "kwarg ",name," is not defined."    

    shortdir = (maindir+file).split('lblecha/')
    if len(shortdir) == 2:
        shortdir = shortdir[1]
    else: shortdir=file

    filebase = ((maindir+file).split('mcrx'))[0]
    if filebase != '': filebase = filebase+'/'
    
    f=pyfits.open(maindir+file)    
    print "\nopened file: ",maindir+file


    init_lam = f['lambda'].data.field('lambda')

    ## clip ends of wavelength range
    #llim = array([25,init_lam.shape[0]-26])
    llim = array([0,init_lam.shape[0]-1])
    
    if max_lum < 2.0 or dynrng_lum < 1.0:
        print "crazy value entered for max_lum or dynrng_lum:", max_lum, dynrng_lum
        quit()

    ## initial values of normalization factors
    sed_norm = 0.0
    sed_scatter_norm = 0.0
    img_vmax = 0.0
    img_vmin = 1.0e6
    img_scatter_vmax = 0.0
    img_scatter_vmin = 1.0e6
    lam_half_rng = 0.0

    ## iterate thru cameras to get data
    for icam in cam_arr:

        print "\nProcessing camera %d..."%icam

        ## process this camera and return relevant values
        pos, hfw, xlim, ylim, kpc_per_pix, lam, sed_noscat, sed_scat, \
            slit_noscat, slit_scat, img_noscat, img_scat, ifu_noscat, ifu_scat, \
            lowlum_sed_noscat = \
            process_onecam(f, snap, icam, init_lam, llim, lambda_rest, 
                           max_lum=max_lum, dynrng_lum=dynrng_lum)

        print "min/max/shape sed_noscat: ",sed_noscat.min(),sed_noscat.max(),sed_noscat.shape
        print "min/max/shape sed_scat: ",sed_scat.min(),sed_scat.max(),sed_scat.shape
        print "min/max/shape slit_noscat: ",slit_noscat.min(),slit_noscat.max(),slit_noscat.shape
        print "min/max/shape slit_scat: ",slit_scat.min(),slit_scat.max(),slit_scat.shape
        print "min/max/shape img_noscat: ",img_noscat.min(),img_noscat.max(),img_noscat.shape
        print "min/max/shape img_scat: ",img_scat.min(),img_scat.max(),img_scat.shape
        #print "min/max/shape ifu_noscat: ",ifu_noscat.min(),ifu_noscat.max(),ifu_noscat.shape
        #print "min/max/shape ifu_scat: ",ifu_scat.min(),ifu_scat.max(),ifu_scat.shape    
        
        ## convert luminosities to solar units
        sed_noscat /= 3.83e33
        sed_scat /= 3.83e33
        slit_noscat -= log10(3.83e33)
        slit_scat -= log10(3.83e33)
        img_noscat -= log10(3.83e33)
        img_scat -= log10(3.83e33)

        #ifu_noscat[img_noscat<(max_lum-dynrng_lum)] = nan
        #ifu_scat[img_scat<(max_lum-dynrng_lum)] = nan
        print "min/max/shape ifu_noscat: ",ifu_noscat.min(),ifu_noscat.max(),ifu_noscat.mean(),median(ifu_noscat),ifu_noscat.shape
        print "min/max/shape ifu_scat: ",ifu_scat.min(),ifu_scat.max(),ifu_scat.mean(),median(ifu_scat),ifu_scat.shape

        ## subtract continuum from 1-D SED if necessary
        ## applies ONLY to 1-d SED!
        if subtract_continuum:

            f_continuum_ns, vfinal_sigma_ns = fit_continuum(lam,sed_noscat, 5, 1.0, 0.02e-7)
            f_continuum, vfinal_sigma = fit_continuum(lam,sed_scat, 5, 1.0, 0.02e-7)

            print shape(f_continuum),shape(vfinal_sigma)            
            #print "min/max lambda=",min(lam),max(lam)
            #print "min/max sed_noscat=",sed_noscat.min(),sed_noscat.max()
            #print "min/max/mean f_continuum_ns=",min(f_continuum_ns),max(f_continuum_ns),mean(f_continuum_ns)
            #print "min/max sed_scat=",min(sed_scat),max(sed_scat)
            #print "min/max/mean f_continuum=",min(f_continuum),max(f_continuum),mean(f_continuum)

            sed_noscat -= f_continuum_ns
            sed_scat -= f_continuum

            print "min/max subtracted sed_noscat=",min(sed_noscat),max(sed_noscat)
            print "min/max subtracted sed_scat=",min(sed_scat),max(sed_scat)
            print "difference b/t non-scatter & scatter SED normalization (min/max/mean):"
            print min(f_continuum_ns-f_continuum),\
                  max(f_continuum_ns-f_continuum),\
                  mean(f_continuum_ns-f_continuum)


        ## are we plotting any data from a no-NLR run?
        lam_nonlr = array([0.0])
        sed_noscat_nonlr = array([0.0])
        sed_scat_nonlr = array([0.0])
        if nonlr_dir != '':
            fno=pyfits.open(nonlr_dir+file)    
            print "\nopened no-NLR file: ",nonlr_dir+file
            lam_nonlr = fno['lambda'].data.field('lambda')

            im_noscat_nonlr = fno['camera%d-nonscatter'%icam].data
            hdr_noscat_nonlr = fno['camera%d-nonscatter'%icam].header
            im_scat_nonlr = fno['camera%d'%icam].data
            hdr_scat_nonlr = fno['camera%d'%icam].header
            hdr_campar_nonlr = fno['camera%d-parameters'%icam].header
            cameradist_nonlr = hdr_campar_nonlr['cameradist'] ## in kpc
            sr_per_pix_nonlr = hdr_scat_nonlr['PIXEL_SR']
            print "cameradist_nonlr=",cameradist_nonlr
            print "sr_per_pix_nonlr=",sr_per_pix_nonlr
            ## code returns surface brightness in W/m/m^2/sr (cameradist is in kpc)
            ## here we convert image unit to lambda*L_lambda in erg/s
            codeSB_to_ergs_nonlr = sr_per_pix_nonlr * 1.0e7 * 4*pi * (cameradist_nonlr*3.09e19)**2
            print "codeSB_to_ergs_nonlr=",codeSB_to_ergs_nonlr
            for l in range(lam_nonlr.shape[0]):
                im_noscat_nonlr[l,:,:] *= lam_nonlr[l]
                im_scat_nonlr[l,:,:] *= lam_nonlr[l]
            im_noscat_nonlr_ergs = im_noscat_nonlr * codeSB_to_ergs_nonlr
            im_scat_nonlr_ergs = im_scat_nonlr * codeSB_to_ergs_nonlr
            print "min/max im_noscat_nonlr_ergs: ",im_noscat_nonlr_ergs[im_noscat_nonlr_ergs>0].min(),\
                  im_noscat_nonlr_ergs[im_noscat_nonlr_ergs>0].max()
            print "min/max im_scat_nonlr_ergs: ",im_scat_nonlr_ergs[im_scat_nonlr_ergs>0].min(), \
                  im_scat_nonlr_ergs[im_scat_nonlr_ergs>0].max()
            im_noscat_nonlr_ergs[im_noscat_nonlr_ergs==0.0] = 1.0e-20
            im_scat_nonlr_ergs[im_scat_nonlr_ergs==0.0] = 1.0e-20
                
            sed_noscat_nonlr = (im_noscat_nonlr_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]]).sum(axis=2).sum(axis=1)
            sed_scat_nonlr = (im_scat_nonlr_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]]).sum(axis=2).sum(axis=1)
            sed_noscat_nonlr /= 3.83e33
            sed_scat_nonlr /= 3.83e33

            ## subtract continuum from 1-D SED if necessary
            ## applies ONLY to 1-d SED!
            if subtract_continuum:

                f_continuum_ns_nonlr, vfinal_sigma_ns_nonlr = fit_continuum(lam_nonlr,sed_noscat_nonlr, 5, 1.0, 0.02e-7)
                f_continuum_nonlr, vfinal_sigma_nonlr = fit_continuum(lam_nonlr,sed_scat_nonlr, 5, 1.0, 0.02e-7)
                
                sed_noscat_nonlr -= f_continuum_ns_nonlr 
                sed_scat_nonlr  -= f_continuum_nonlr 

                print "min/max subtracted sed_noscat_nonlr=",min(sed_noscat_nonlr),max(sed_noscat_nonlr)
                print "min/max subtracted sed_scat_nonlr=",min(sed_scat_nonlr),max(sed_scat_nonlr) 
                print "difference b/t nonlr noscat & scatter SED normalization (min/max/mean):"
                print min(f_continuum_ns_nonlr-f_continuum_nonlr),\
                      max(f_continuum_ns_nonlr-f_continuum_nonlr),\
                      mean(f_continuum_ns_nonlr-f_continuum_nonlr)                

            fno.close()


        ## deal with whole-SED data, if any
        lam_ws = array([0.0])
        sed_noscat_ws = array([0.0])
        sed_scat_ws = array([0.0])
        if wholesed_dir != '':
            fws=pyfits.open(wholesed_dir+file)    
            print "\nopened whole-SED file: ",wholesed_dir+file
            lam_ws = fws['lambda'].data.field('lambda')

            im_noscat_ws = fws['camera%d-nonscatter'%icam].data
            hdr_noscat_ws = fws['camera%d-nonscatter'%icam].header
            im_scat_ws = fws['camera%d'%icam].data
            hdr_scat_ws = fws['camera%d'%icam].header
            hdr_campar_ws = fws['camera%d-parameters'%icam].header
            cameradist_ws = hdr_campar_ws['cameradist'] ## in kpc
            sr_per_pix_ws = hdr_scat_ws['PIXEL_SR']
            print "cameradist_ws=",cameradist_ws
            print "sr_per_pix_ws=",sr_per_pix_ws
            ## code returns surface brightness in W/m/m^2/sr (cameradist is in kpc)
            ## here we convert image unit to lambda*L_lambda in erg/s
            codeSB_to_ergs_ws = sr_per_pix_ws * 1.0e7 * 4*pi * (cameradist_ws*3.09e19)**2
            print "codeSB_to_ergs_ws=",codeSB_to_ergs_ws
            for l in range(lam_ws.shape[0]):
                im_noscat_ws[l,:,:] *= lam_ws[l]
                im_scat_ws[l,:,:] *= lam_ws[l]
            im_noscat_ws_ergs = im_noscat_ws * codeSB_to_ergs_ws
            im_scat_ws_ergs = im_scat_ws * codeSB_to_ergs_ws
            print "min/max im_noscat_ws_ergs: ",im_noscat_ws_ergs[im_noscat_ws_ergs>0].min(),\
                  im_noscat_ws_ergs[im_noscat_ws_ergs>0].max()
            print "min/max im_scat_ws_ergs: ",im_scat_ws_ergs[im_scat_ws_ergs>0].min(),\
                  im_scat_ws_ergs[im_scat_ws_ergs>0].max()
            im_noscat_ws_ergs[im_noscat_ws_ergs==0.0] = 1.0e-20
            im_scat_ws_ergs[im_scat_ws_ergs==0.0] = 1.0e-20
            
            sed_noscat_ws = (im_noscat_ws_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]]).sum(axis=2).sum(axis=1)
            sed_scat_ws = (im_scat_ws_ergs[:,ylim[0]:ylim[1],xlim[0]:xlim[1]]).sum(axis=2).sum(axis=1)
            sed_noscat_ws /= 3.83e33
            sed_scat_ws /= 3.83e33

            ## subtract continuum from 1-D SED if necessary
            ## applies ONLY to 1-d SED!
            if subtract_continuum:

                f_continuum_ns_ws, vfinal_sigma_ns_ws = fit_continuum(lam_ws,sed_noscat_ws, 5, 1.0, 0.02e-7)
                f_continuum_ws, vfinal_sigma_ws = fit_continuum(lam_ws,sed_scat_ws, 5, 1.0, 0.02e-7)
                
                sed_noscat_ws -= f_continuum_ns_ws 
                sed_scat_ws  -= f_continuum_ws 
                
                print "min/max subtracted sed_noscat_ws=",min(sed_noscat_ws),max(sed_noscat_ws)
                print "min/max subtracted sed_scat_ws=",min(sed_scat_ws),max(sed_scat_ws)
                print "difference b/t whole-sed noscat & scatter SED normalization (min/max/mean):"
                print min(f_continuum_ns_ws-f_continuum_ws),\
                      max(f_continuum_ns_ws-f_continuum_ws),\
                      mean(f_continuum_ns_ws-f_continuum_ws)            

            fws.close()

        
        if icam == 0:
            x_kpc_per_pix = kpc_per_pix[0]
            y_kpc_per_pix = kpc_per_pix[1]
            
            pltpos = ( pos - ylim[0] - 0.5*(ylim[1]-ylim[0]) )*y_kpc_per_pix
            plthfw = hfw * y_kpc_per_pix

            xran = array([(-0.5*(xlim[1]-xlim[0]))*x_kpc_per_pix, (0.5*(xlim[1]-xlim[0]))*x_kpc_per_pix])
            yran = array([(-0.5*(ylim[1]-ylim[0]))*y_kpc_per_pix, (0.5*(ylim[1]-ylim[0]))*y_kpc_per_pix])

            allcam_sed_noscat =  zeros((len(cam_arr), sed_noscat.shape[0]), dtype=float64)
            allcam_sed_scat =    zeros((len(cam_arr), sed_scat.shape[0]), dtype=float64)
            allcam_slit_noscat = zeros((len(cam_arr), slit_noscat[:,0].shape[0], slit_noscat[0,:].shape[0]), dtype=float64)
            allcam_slit_scat =   zeros((len(cam_arr), slit_scat[:,0].shape[0], slit_scat[0,:].shape[0]), dtype=float64)
            allcam_img_noscat =  zeros((len(cam_arr), img_noscat[:,0].shape[0], img_noscat[0,:].shape[0]), dtype=float64)
            allcam_img_scat =    zeros((len(cam_arr), img_scat[:,0].shape[0], img_scat[0,:].shape[0]), dtype=float64)
            allcam_ifu_noscat =  zeros((len(cam_arr), ifu_noscat[:,0].shape[0], ifu_noscat[0,:].shape[0]), dtype=float64)
            allcam_ifu_scat =    zeros((len(cam_arr), ifu_scat[:,0].shape[0], ifu_scat[0,:].shape[0]), dtype=float64)
            print "sed array shapes: ", allcam_sed_noscat.shape, allcam_sed_scat.shape
            print "slit array shapes: ", allcam_slit_noscat.shape, allcam_slit_scat.shape
            print "img array shapes: ", allcam_img_noscat.shape, allcam_img_scat.shape
            print "ifu array shapes: ", allcam_ifu_noscat.shape, allcam_ifu_scat.shape

            allcam_lowlum_sed_noscat = zeros((len(cam_arr), lowlum_sed_noscat.shape[0]), dtype=float64)

            if nonlr_dir != '':
                allcam_sed_noscat_nonlr = zeros((len(cam_arr), sed_noscat_nonlr.shape[0]), dtype=float64)
                allcam_sed_scat_nonlr = zeros((len(cam_arr), sed_scat_nonlr.shape[0]), dtype=float64)
            if wholesed_dir != '':
                allcam_sed_noscat_ws = zeros((len(cam_arr), sed_noscat_ws.shape[0]), dtype=float64)
                allcam_sed_scat_ws = zeros((len(cam_arr), sed_scat_ws.shape[0]), dtype=float64)

                
        ## set normalizations for axis and color scalings
        if subtract_continuum:
            sed_norm = max([sed_norm, f_continuum_ns.max(), 1.05*max([sed_noscat.max(),sed_scat.max()])])
            sed_scatter_norm = max( [sed_scatter_norm, f_continuum.max(), 1.05*sed_scat.max() ])
        else:
            sed_norm = max([sed_norm, 1.0e9, 1.05*max([sed_noscat.max(),sed_scat.max()])])
            sed_scatter_norm = max( [sed_scatter_norm, 1.0e9, 1.05*sed_scat.max() ])

        
        img_vmax = max([ img_vmax, min([ max_lum, max([slit_noscat.max(), img_noscat.max(),
                                                      slit_scat.max(), img_noscat.max()]) ]) ])
        img_vmin = min([ img_vmin, 
                         max([ 1.0, img_vmax - dynrng_lum, 
                               min([slit_noscat.min(), img_noscat.min(),
                                    slit_scat.min(), img_noscat.min()]) ]) ])
        print img_vmax,img_vmin
        
        img_scatter_vmax = max([ img_scatter_vmax, 
                                 min([ max_lum, max([slit_scat.max(), img_noscat.max()]) ]) ])
        img_scatter_vmin = min([ img_scatter_vmin, 
                                 max([ 1.0, img_scatter_vmax - dynrng_lum,
                                       min([slit_scat.min(), img_noscat.min()]) ]) ])
        print img_scatter_vmax,img_scatter_vmin
        
        ## convert to Angstrom, find normalizations, & set rest wavelength to color scale midpoint
        tmp_vmax = min([ lam[-1]*1.0e10, max([ifu_noscat[img_noscat>=img_vmin].max(),
                                                   ifu_scat[img_scat>=img_vmin].max()]) ])
        tmp_vmin = max([ lam[0]*1.0e10, min([ifu_noscat[img_noscat>=img_vmin].min(),
                                                   ifu_scat[img_scat>=img_vmin].min()]) ])
        #tmp_vmax = min([ lam[llim[1]]*1.0e10, max([ifu_noscat[img_noscat>=img_vmin].max(),
        #                                           ifu_scat[img_scat>=img_vmin].max()]) ])
        #tmp_vmin = max([ lam[llim[0]]*1.0e10, min([ifu_noscat[img_noscat>=img_vmin].min(),
        #                                           ifu_scat[img_scat>=img_vmin].min()]) ])
        if tmp_vmax < lambda_rest or tmp_vmin > lambda_rest:
            print "Error: ifu wavelength range does not include rest wavelength!"
            print "range = (",tmp_vmin,":",tmp_vmax,")"
            print "max ifu_noscat=",ifu_noscat.max()
            print "max ifu_noscat cut=",ifu_noscat[img_noscat>=img_vmin].max()
            print "max ifu_scat=",ifu_scat.max()
            print "max ifu_scat cut=",ifu_scat[img_scat>=img_vmin].max()
            quit()

        lam_half_rng = max([ lam_half_rng, 
                             max([abs(tmp_vmax-lambda_rest),
                                  abs(lambda_rest-tmp_vmin)]) ])

            

        allcam_sed_noscat[icam,:] = sed_noscat
        allcam_sed_scat[icam,:] = sed_scat
        allcam_slit_noscat[icam,:,:] = slit_noscat
        allcam_slit_scat[icam,:,:] = slit_scat
        allcam_img_noscat[icam,:,:] = img_noscat
        allcam_img_scat[icam,:,:] = img_scat
        allcam_ifu_noscat[icam,:,:] = ifu_noscat
        allcam_ifu_scat[icam,:,:] = ifu_scat
        del sed_noscat, sed_scat, slit_noscat, slit_scat
        del img_noscat, img_scat, ifu_noscat, ifu_scat

        if nonlr_dir != '':
            allcam_sed_noscat_nonlr[icam,:] = sed_noscat_nonlr
            allcam_sed_scat_nonlr[icam,:] = sed_scat_nonlr
            del sed_noscat_nonlr, sed_scat_nonlr
        if wholesed_dir != '':
            allcam_sed_noscat_ws[icam,:] = sed_noscat_ws
            allcam_sed_scat_ws[icam,:] = sed_scat_ws
            del sed_noscat_ws, sed_scat_ws

        allcam_lowlum_sed_noscat[icam,:] = lowlum_sed_noscat
        del lowlum_sed_noscat

    ## done with fits file
    f.close()

    ### TEMPORARY -- DEBUG ONLY
    #plt.clf()
    ##fig, axarr = plt.subplots(len(cam_arr),sharex=True)
    #for i in cam_arr:
    #    print "plotting camera", i,"..."
    #    plt.subplot(7,1,i)
    #    plt.plot(lam,allcam_lowlum_sed_noscat[i,:])
    #    plt.plot([lambda_rest*1.0e-10,lambda_rest*1.0e-10],[0.0,allcam_lowlum_sed_noscat.max()],'k')

    #    Lblue_lowlum = allcam_lowlum_sed_noscat[i,init_lam<1.0e-10*lambda_rest].sum()
    #    Lred_lowlum = allcam_lowlum_sed_noscat[i,init_lam>1.0e-10*lambda_rest].sum()
    #    Ltot_lowlum = allcam_lowlum_sed_noscat[i,:].sum()
    #    print "Lblue_lowlum, Lred_lowlum, Ltot_lowlum = ",Lblue_lowlum, Lred_lowlum, Ltot_lowlum

    #plt.savefig("lowlum_sed_test.eps")
    ###

    ## set overall normalizations
    print "sed scale range = (0, ",sed_norm,")"
    print "sed scatter scale range = (0, ",sed_scatter_norm,")"

    #img_norm = colors.Normalize(vmin=img_vmin,vmax=img_vmax)
    img_norm = array([img_vmin,img_vmax])
    print "(log) luminosity scale range = ",img_norm
    print "(log) luminosity scatter scale range [NOT CURRENTLY USED]= (",img_scatter_vmin,":",img_scatter_vmax,")"


    #DEBUG -- doing this manually!!
    lam_half_rng = 10.0 ## Angstroms
    ifu_vmax = lambda_rest + lam_half_rng
    ifu_vmin = lambda_rest - lam_half_rng
    ifu_norm = array([ifu_vmin, ifu_vmax])
    print "ifu color scale range = (",ifu_norm[0],":",ifu_norm[1],")"

    lran = ifu_norm

    ## clip low-lum values
    allcam_ifu_noscat[allcam_img_noscat<img_vmin] = nan
    allcam_img_noscat[allcam_img_noscat<img_vmin] = nan
    allcam_ifu_scat[allcam_img_scat<img_vmin] = nan
    allcam_img_scat[allcam_img_scat<img_vmin] = nan

    ## set up the plot
    #plt.ioff()
    plt.clf()
    fig = plt.figure(figsize=(9.6,8))

    row_size = 0.85/len(cam_arr)
    wspace = 0.3
    hspace = 0.2

    cbaxes_lum = fig.add_axes([0.31,0.05,0.42,0.02])
    cbaxes_ifu = fig.add_axes([0.78,0.05,0.19,0.02])
    #cbaxes_lum = [0.3,0.05,0.44,0.02]
    #cbaxes_ifu = [0.77,0.05,0.21,0.02]

    ## iterate thru cameras again to make nonscatter plot
    for icam in cam_arr:

        print "\nPlotting camera %d..."%icam

        if nonlr_dir != '':
            thiscam_sed_noscat_nonlr = allcam_sed_noscat_nonlr[icam,:]
            thiscam_sed_scat_nonlr = allcam_sed_scat_nonlr[icam,:]
        else:
            thiscam_sed_noscat_nonlr = array([0.0])
            thiscam_sed_scat_nonlr = array([0.0])
        if wholesed_dir != '':
            thiscam_sed_noscat_ws = allcam_sed_noscat_ws[icam,:]
            thiscam_sed_scat_ws = allcam_sed_scat_ws[icam,:]
        else:
            thiscam_sed_noscat_ws = array([0.0])
            thiscam_sed_scat_ws = array([0.0])


        ## plot rows for each camera
        plot_onecam_row(icam, len(cam_arr), lam, lran, lambda_rest, xran, yran,
                        row_size, wspace, hspace,
                        sed_norm, sed_scatter_norm, img_norm, pltpos, plthfw,
                        cbaxes_lum, cbaxes_ifu,
                        lam_nonlr, thiscam_sed_noscat_nonlr, thiscam_sed_scat_nonlr,
                        lam_ws, thiscam_sed_noscat_ws, thiscam_sed_scat_ws,
                        allcam_sed_noscat[icam,:], allcam_slit_noscat[icam,:,:],
                        allcam_img_noscat[icam,:,:], allcam_ifu_noscat[icam,:,:],
                        allcam_sed_scat[icam,:], allcam_slit_scat[icam,:,:],
                        allcam_img_scat[icam,:,:], allcam_ifu_scat[icam,:,:])

    
    ## finish and save the plot
    plt.suptitle('%s'%shortdir,fontsize=10)
    print "dpi:", fig.get_dpi()
    #plt.savefig("%d-multicam_hbetaspec.pdf"%snap)
    if subtract_continuum: filename=filebase+"%d-multicam_hbetaspec.eps"%snap
    else: filename=filebase+"%d-multicam_hbetaspec_withcontinuum.eps"%snap
    plt.savefig(filename)

        
def plot_onecam_row(icam, ncam, lam, lran, lambda_rest, xran, yran,
                    row_size, wspace, hspace,
                    sed_norm, sed_scatter_norm, img_norm, pltpos, plthfw,
                    cbaxes_lum, cbaxes_ifu,
                    lam_nonlr, sed_noscat_nonlr, sed_scat_nonlr,
                    lam_ws, sed_noscat_ws, sed_scat_ws,
                    sed_noscat, slit_noscat, img_noscat, ifu_noscat,
                    sed_scat, slit_scat, img_scat, ifu_scat):
    
    
    print "min/max/shape sed_noscat: ",sed_noscat.min(),sed_noscat.max(),sed_noscat.shape
    print "min/max/shape sed_scat: ",sed_scat.min(),sed_scat.max(),sed_scat.shape
    print "min/max/shape slit_noscat: ",slit_noscat.min(),slit_noscat.max(),slit_noscat.shape
    print "min/max/shape slit_scat: ",slit_scat.min(),slit_scat.max(),slit_scat.shape
    print "min/max/shape img_noscat: ",img_noscat.min(),img_noscat.max(),img_noscat.shape
    print "min/max/shape img_scat: ",img_scat.min(),img_scat.max(),img_scat.shape
    print "min/max/shape ifu_noscat: ",ifu_noscat.min(),ifu_noscat.max(),ifu_noscat.shape
    print "min/max/shape ifu_scat: ",ifu_scat.min(),ifu_scat.max(),ifu_scat.shape    
    

    ## get intrinsic aspect ratios for axes
    slit_aspect = (lran[1]-lran[0])/float(yran[1]-yran[0])
    img_aspect = float(yran[1]-yran[0])/(xran[1]-xran[0])
    
    gs0 = gridspec.GridSpec(1,8)
    gs0.update(left=0.06,right=0.98, bottom=1.0-0.04-(icam+1)*(row_size),
               top=1.0-0.04-icam*(row_size)-0.015, wspace=wspace, hspace=hspace)

    ### NOSCATTER SED
    ax0 = plt.subplot(gs0[0])
    plt.plot(array([lambda_rest,lambda_rest]),array([0,sed_norm]),'black',linestyle=':')
    SEDNOSCAT = plt.plot(lam*1.0e10, sed_noscat)
    plt.plot(lam_nonlr*1.0e10, sed_noscat_nonlr,'magenta',linestyle='--')
    plt.plot(lam_ws*1.0e10, sed_noscat_ws,'green',linestyle='--')
    plt.xlim(lran[0],lran[1])
    plt.ylim(0,sed_norm)
    if icam == ncam-1:
        plt.xlabel(r'$\lambda$ [$\AA$]',fontsize=8)
        ax0.xaxis.set_major_locator(ticker.MaxNLocator(4))
        for tick in ax0.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax0.set_xticklabels([])
    plt.ylabel(('CAMERA%d\n'%icam)+(r'$\lambda L_{\lambda}$ [$L_{\odot}$]'),fontsize=8)
    if icam == 0: plt.title('non-scatter',fontsize=9)
    ax0.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax0.ticklabel_format(axis='y',style='sci',scilimits=(-3,4))
    for tick in ax0.yaxis.get_major_ticks():
        tick.label1.set_fontsize(8)
    ax0.yaxis.get_offset_text().set_fontsize(8)

    ### SCATTER SED
    ax4 = plt.subplot(gs0[1])
    plt.plot(array([lambda_rest,lambda_rest]),array([0,sed_scatter_norm]),'black',linestyle=':')
    SEDSCAT = plt.plot(lam*1.0e10, sed_scat)
    plt.plot(lam_nonlr*1.0e10, sed_scat_nonlr,'magenta',linestyle='--')
    plt.plot(lam_ws*1.0e10, sed_scat_ws,'green',linestyle='--')
    plt.xlim(lran[0],lran[1])
    #plt.ylim(0,sed_norm)
    plt.ylim(0,sed_scatter_norm)
    if icam == ncam-1:
        plt.xlabel(r'$\lambda$ [$\AA$]',fontsize=8)
        ax4.xaxis.set_major_locator(ticker.MaxNLocator(4))
        for tick in ax4.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax4.set_xticklabels([])
    #plt.ylabel(r'$\lambda L_{\lambda}$ [erg s$^{-1}$]')
    if icam == 0: plt.title('scatter',fontsize=9)
    #ax4.set_yticklabels([])
    ax4.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax4.ticklabel_format(axis='y',style='sci',scilimits=(-3,4))
    for tick in ax4.yaxis.get_major_ticks():
        tick.label1.set_fontsize(8)
    ax4.yaxis.get_offset_text().set_fontsize(8)



    ### NOSCATTER SLIT SPECTRUM
    ax1 = plt.subplot(gs0[2])
    SLITNOSCAT = plt.imshow(slit_noscat, extent=array([xran,lran]).flatten(),
               cmap='hot',norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),origin='lower',
               interpolation='nearest',aspect='auto')
    #interpolation='nearest',aspect=0.5/slit_aspect)
    if icam == 0: plt.title('non-scatter',fontsize=10)
    plt.ylabel(r'$\lambda$ [$\AA$]',fontsize=8)
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontsize(8)
    if icam == ncam-1:
        plt.xlabel("x [kpc]",fontsize=8)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax1.set_xticklabels([])


    ### SCATTER SLIT SPECTRUM
    ax5 = plt.subplot(gs0[3])
    SLITSCAT = plt.imshow(slit_scat, extent=array([xran,lran]).flatten(),
               cmap='hot',norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),origin='lower',
               interpolation='nearest',aspect='auto')
    #interpolation='nearest',aspect=0.5/slit_aspect)
    if icam == 0: plt.title('scatter',fontsize=10)
    #plt.ylabel(r'$\lambda$ [$\AA$]')
    #ax5.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax5.set_yticklabels([])
    if icam == ncam-1:
        plt.xlabel("x [kpc]",fontsize=8)
        ax5.ticklabel_format(axis='x')
        for tick in ax5.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax5.set_xticklabels([])


    ### NOSCATTER HBETA IMAGE
    ax2 = plt.subplot(gs0[4])
    plt.plot(xran*0.99,array([pltpos-plthfw,pltpos-plthfw]),'white',linewidth=1.5)
    plt.plot(xran*0.99,array([pltpos+plthfw,pltpos+plthfw]),'white',linewidth=1.5)
    IMNOSCAT = plt.imshow(img_noscat,interpolation='nearest',norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),
                          cmap='hot',extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    ax2.patch.set_facecolor('black')
    if icam == 0: plt.title('non-scatter',fontsize=10)
    if icam == ncam-1:
        plt.xlabel("x [kpc]",fontsize=8)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax2.set_xticklabels([])
    plt.ylabel("y [kpc]",fontsize=8)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label1.set_fontsize(8)

    ### SCATTER HBETA IMAGE
    ax6 = plt.subplot(gs0[5])
    plt.plot(xran*0.99,array([pltpos-plthfw,pltpos-plthfw]),'white',linewidth=1.5)
    plt.plot(xran*0.99,array([pltpos+plthfw,pltpos+plthfw]),'white',linewidth=1.5)
    IMSCAT = plt.imshow(img_scat,interpolation='nearest',norm=colors.Normalize(vmin=img_norm[0],vmax=img_norm[1]),
                        cmap='hot',extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    ax6.patch.set_facecolor('black')
    if icam == 0: plt.title('scatter',fontsize=10)
    if icam == ncam-1:
        plt.xlabel("x [kpc]",fontsize=8)
        for tick in ax6.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax6.set_xticklabels([])
    ax6.set_yticklabels([])
    #plt.ylabel("y [kpc]")

    if icam==ncam-1:
        ## colorbar for luminosity scale
        cb = plt.colorbar(cax = cbaxes_lum, orientation='horizontal')
        #cb = plt.colorbar(IMNOSCAT, orientation='horizontal')
        cb.set_label(r'$\lambda L_{\lambda}$ [L$_{\odot}$]',fontsize=8)
        #cb.ax.set_position(cbaxes_lum)
        #cb.ax.xaxis.set_major_locator(ticker.MaxNLocator(8))
        cb.ax.tick_params(labelsize=8)


    ### NOSCATTER IFU
    ## try colormaps bwr, coolwarm, or seismic instead of jet

    #norm=colors.Normalize(vmin=lran[0]+0.1*(lran[1]-lran[0]),
    #                      vmax=lran[1]-0.1*(lran[1]-lran[0])),

    ax3 = plt.subplot(gs0[6])
    ### TEST ###
    ##ifu_noscat = zeros((300,300)) + lambda_rest + 2.0
    ### TEST - DEBUG ###
    lvmin=lran[0]-5.0
    lvmax=lran[1]+5.0
    print "lvmin,lvmax=",lvmin,lvmax
    IFUNOSCAT = plt.imshow(ifu_noscat,interpolation='nearest',cmap='seismic',
                           norm=colors.Normalize(vmin=lvmin,vmax=lvmax),
                           extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    #IFUNOSCAT = plt.imshow(ifu_noscat,interpolation='nearest',cmap='seismic',
    #                       extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    if icam == 0: plt.title('non-scatter',fontsize=10)
    if icam == ncam-1:
        plt.xlabel("x [kpc]",fontsize=8)
        for tick in ax3.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax3.set_xticklabels([])
    plt.ylabel("y [kpc]",fontsize=8)
    for tick in ax3.yaxis.get_major_ticks():
        tick.label1.set_fontsize(8)


    ### SCATTER IFU
    ax7 = plt.subplot(gs0[7])
    ### TEST ###
    ##ifu_scat = zeros((300,300)) + lambda_rest + 1.0
    IFUSCAT = plt.imshow(ifu_scat,interpolation='nearest',cmap='seismic',
                         norm=colors.Normalize(vmin=lran[0],vmax=lran[1]),
                         extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    #IFUSCAT = plt.imshow(ifu_scat,interpolation='nearest',cmap='seismic',
    #                     extent=array([xran,yran]).flatten(),origin='lower',aspect='auto')
    if icam == 0: plt.title('scatter',fontsize=10)
    if icam == ncam-1:
        plt.xlabel("x [kpc]",fontsize=8)
        for tick in ax7.xaxis.get_major_ticks():
            tick.label1.set_fontsize(8)
    else: ax7.set_xticklabels([])
    ax7.set_yticklabels([])
    #plt.ylabel("y [kpc]")

    if icam==ncam-1:
        ## colorbar for wavelength scale
        cb2 = plt.colorbar(cax = cbaxes_ifu, orientation='horizontal',ticks=[lran[0],lambda_rest,lran[1]])
        #cb2 = plt.colorbar(cax = cbaxes_ifu, orientation='horizontal',ticks=[lran[0]+0.1*(lran[1]-lran[0]),
        #                                                                     lambda_rest,lran[1]-0.1*(lran[1]-lran[0])])
        #cb2 = plt.colorbar(IFUNOSCAT, orientation='horizontal',ticks=[lran[0],lambda_rest,lran[1]])
        cb2.set_label(r'$\lambda$ [$\AA$]',fontsize=8)
        #cb2.ax.set_position(cbaxes_ifu)
        #cb2.ax.set_xticklabels([ "%.1f"%((lran[0])+0.1*(lran[1]-lran[0])), "%.1f"%lambda_rest,
        #                         "%.1f"%((lran[1])-0.1*(lran[1]-lran[0])) ])
        cb2.ax.tick_params(labelsize=8)



def mytest(inval, outval):

    return inval*3, arange(24).reshape(3,4,2), outval*-4.7



if __name__ == "__main__":
    nonlr_dir=""
    wholesed_dir=""    
    if len(sys.argv)==3:
        multicam(sys.argv[1],long(sys.argv[2]))
    elif len(sys.argv)>3:
        if len(sys.argv)>=5:
            nonlr_dir=sys.argv[4]
            if len(sys.argv)==6: wholesed_dir=sys.argv[5]        
        
        if sys.argv[3]=="False":
            multicam(sys.argv[1],long(sys.argv[2]),subtract_continuum=False,nonlr_dir=nonlr_dir,wholesed_dir=wholesed_dir)
        else:
            multicam(sys.argv[1],long(sys.argv[2]),nonlr_dir=nonlr_dir,wholesed_dir=wholesed_dir)

