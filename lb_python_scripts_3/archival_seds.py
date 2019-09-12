import matplotlib
#matplotlib.use('Agg')
import pyfits, glob, pylab, numpy, copy as cpy, stats, vecops, types, h5py, re
import sys
import pj_constants as pjc
from pylab import *
import scipy.io
import broadband as bb

## Jarrett et al. 2011                       
W1_Vega_to_AB = 2.699
W2_Vega_to_AB = 3.339
W3_Vega_to_AB = 5.174
W4_Vega_to_AB = 6.620

def lanz_plots(path='/oasis/projects/nsf/hvd115/lblecha/lanz_2014_data',
               AISM=False):

    global W1_Vega_to_AB, W2_Vega_to_AB, W3_Vega_to_AB, W4_Vega_to_AB

    if AISM:
        data = scipy.io.readsav('%s/SIGS_sims_I_AISM_SEDs.sav'%path)
    else:
        data = scipy.io.readsav('%s/SIGS_sims_I_DISM_SEDs.sav'%path)

    lam = data['lambda']
    Llam = data['l_lambda']
    time = data['time']
    sfr = data['sfr']
    d_bh = data['d_bh']
    mgas = data['mgas']
    mdust = data['mdust']
    mstar = data['mstar']
    fgas = mgas/(mgas+mstar)
    fdust = mdust/(mdust+mgas)
    runnames = data['runname'].tolist()
    print runnames
    #print np.where(time[:,-1]>0)[0]
    iso_plotidx = [runnames.index(name) for name in ['M0','M1','M2','M3']]
    equal_plotidx = [runnames.index(name) for name in ['M0M0e','M1M1e','M2M2e','M3M3e']]
    uneq_plotidx = [runnames.index(name) for name in ['M1M0e','M2M0e','M2M1e',
                                                             'M3M0e','M3M1e','M3M2e']]

    IC_mstar = {'M0':0.061, 'M1':0.38, 'M2':1.18, 'M3':4.22}
    IC_mgas =  {'M0':0.035, 'M1':0.14, 'M2':0.33, 'M3':0.80}
    IC_mtot =  {'M0':5.0, 'M1':20.0, 'M2':51.0, 'M3':116.0}
    IC_mstar.update({'M0M0e':(IC_mstar['M0']*2), 'M1M1e':(IC_mstar['M1']*2),
                     'M2M2e':(IC_mstar['M2']*2), 'M3M3e':(IC_mstar['M3']*2),
                     'M1M0e':(IC_mstar['M1']+IC_mstar['M0']), 'M2M0e':(IC_mstar['M2']+IC_mstar['M0']),
                     'M2M1e':(IC_mstar['M2']+IC_mstar['M1']), 'M3M0e':(IC_mstar['M3']+IC_mstar['M0']),
                     'M3M1e':(IC_mstar['M3']+IC_mstar['M1']), 'M3M2e':(IC_mstar['M3']+IC_mstar['M2'])})
    IC_mgas.update({'M0M0e':(IC_mgas['M0']*2), 'M1M1e':(IC_mgas['M1']*2),
                    'M2M2e':(IC_mgas['M2']*2), 'M3M3e':(IC_mgas['M3']*2),
                    'M1M0e':(IC_mgas['M1']+IC_mgas['M0']), 'M2M0e':(IC_mgas['M2']+IC_mgas['M0']),
                    'M2M1e':(IC_mgas['M2']+IC_mgas['M1']), 'M3M0e':(IC_mgas['M3']+IC_mgas['M0']),
                    'M3M1e':(IC_mgas['M3']+IC_mgas['M1']), 'M3M2e':(IC_mgas['M3']+IC_mgas['M2'])})
    IC_mtot.update({'M0M0e':(IC_mtot['M0']*2), 'M1M1e':(IC_mtot['M1']*2),
                    'M2M2e':(IC_mtot['M2']*2), 'M3M3e':(IC_mtot['M3']*2),
                    'M1M0e':(IC_mtot['M1']+IC_mtot['M0']), 'M2M0e':(IC_mtot['M2']+IC_mtot['M0']),
                    'M2M1e':(IC_mtot['M2']+IC_mtot['M1']), 'M3M0e':(IC_mtot['M3']+IC_mtot['M0']),
                    'M3M1e':(IC_mtot['M3']+IC_mtot['M1']), 'M3M2e':(IC_mtot['M3']+IC_mtot['M2'])})
    IC_fgas = {}
    IC_fgas.update({key:IC_mgas[key]/(IC_mgas[key]+IC_mstar[key]) for key in IC_mgas})
    print "initial fgas: "
    for key,fg in IC_fgas.iteritems():
        print '%s: %.3g'%(key,fg)
    IC_q = {'M0':0.0, 'M1':0.0, 'M2':0.0, 'M3':0.0}
    IC_q.update({'M0M0e':1.0, 'M1M1e':1.0, 'M2M2e':1.0, 'M3M3e':1.0,
                 'M1M0e':(IC_mtot['M0']/IC_mtot['M1']), 'M2M0e':(IC_mtot['M0']/IC_mtot['M2']),
                 'M2M1e':(IC_mtot['M1']/IC_mtot['M2']), 'M3M0e':(IC_mtot['M0']/IC_mtot['M3']),
                 'M3M1e':(IC_mtot['M1']/IC_mtot['M3']), 'M3M2e':(IC_mtot['M2']/IC_mtot['M3'])})
    print "q for unequal-mass mergers: "
    for key,q in IC_q.iteritems():
        print '%s %.3g'%(key,q)
    #print [IC_q[runnames[n]] for n in uneq_plotidx]
    plt.clf()
    plt.close('all')
    extra='AISM_' if AISM else ''
    oplotvals_vs_t(path,'iso_%ssfr_bhsep'%extra,runnames,iso_plotidx,
                   time,d_bh,sfr,val3=sfr/mstar, 
                   ylim1=(.1,250),ylim2=(1.0e-3,100),ylim3=(1.0e-12,1.0e-9),
                   ylabel1='BH sep [kpc]',ylabel2='SFR [Msun/yr]',ylabel3='sSFR [yr$^{-1}$]')
    oplotvals_vs_t(path,'equal_%ssfr_bhsep'%extra,runnames,equal_plotidx,
                   time,d_bh,sfr,val3=sfr/mstar, 
                   ylim1=(.1,250),ylim2=(1.0e-3,100),ylim3=(1.0e-12,1.0e-9),
                   ylabel1='BH sep [kpc]',ylabel2='SFR [Msun/yr]',ylabel3='sSFR [yr$^{-1}$]')
    oplotvals_vs_t(path,'uneq_%ssfr_bhsep'%extra,runnames,uneq_plotidx,
                   time,d_bh,sfr,val3=sfr/mstar, 
                   ylim1=(.1,250),ylim2=(1.0e-3,100),ylim3=(1.0e-12,1.0e-9),
                   ylabel1='BH sep [kpc]',ylabel2='SFR [Msun/yr]',ylabel3='sSFR [yr$^{-1}$]')

    oplotvals_vs_t(path,'iso_%sfgas_fdust'%extra,runnames,iso_plotidx,
                   time,fgas,fdust,ylim1=(0,0.4),ylim2=(0,0.02),ylog1=False,ylog2=False,
                   ylabel1='Mgas/(Mgas+M*)',ylabel2='Mdust/(Mdust+Mgas)')
    oplotvals_vs_t(path,'equal_%sfgas_fdust'%extra,runnames,equal_plotidx,
                   time,fgas,fdust,ylim1=(0,0.4),ylim2=(0,0.02),ylog1=False,ylog2=False,
                   ylabel1='Mgas/(Mgas+M*)',ylabel2='Mdust/(Mdust+Mgas)')
    oplotvals_vs_t(path,'uneq_%sfgas_fdust'%extra,runnames,uneq_plotidx,
                   time,fgas,fdust,ylim1=(0,0.4),ylim2=(0,0.02),ylog1=False,ylog2=False,
                   ylabel1='Mgas/(Mgas+M*)',ylabel2='Mdust/(Mdust+Mgas)')


    L_W1 = np.zeros(Llam.shape[1:])*np.nan
    L_W2 = np.zeros(Llam.shape[1:])*np.nan
    L_W3 = np.zeros(Llam.shape[1:])*np.nan
    MAB_W1 = np.zeros(Llam.shape[1:])*np.nan
    MAB_W2 = np.zeros(Llam.shape[1:])*np.nan
    MAB_W3 = np.zeros(Llam.shape[1:])*np.nan
    #W1W2_Vega = np.zeros(Llam.shape[1:]) - 1
    #W1W2_Vega_cammean = np.zeros(Llam.shape[2:]) - 1
    for ksim in range(Llam.shape[3]):
        print "sim: %s"%data['runname'][ksim]
        for jt in range(Llam.shape[2]):
            if time[jt,ksim]==0: continue
            for icam in range(Llam.shape[1]):
                L_W1[icam,jt,ksim],MAB_W1[icam,jt,ksim] = bb.Leff('WISE-W1',lam,Llam[:,icam,jt,ksim],
                                                                  sed_lam_units='m',sed_lum_units='W',
                                                                  return_ABMag=True)
                L_W2[icam,jt,ksim],MAB_W2[icam,jt,ksim] = bb.Leff('WISE-W2',lam,Llam[:,icam,jt,ksim],
                                                                  sed_lam_units='m',sed_lum_units='W',
                                                                  return_ABMag=True)
                L_W3[icam,jt,ksim],MAB_W3[icam,jt,ksim] = bb.Leff('WISE-W3',lam,Llam[:,icam,jt,ksim],
                                                                  sed_lam_units='m',sed_lum_units='W',
                                                                  return_ABMag=True)

    W1W2_Vega = (MAB_W1 - W1_Vega_to_AB) - (MAB_W2 - W2_Vega_to_AB)
    W1W2_Vega_cammean = np.mean(W1W2_Vega,axis=0)
    W1W2_Vega_cammax = np.max(W1W2_Vega,axis=0)
    W1W2_Vega_cammin = np.min(W1W2_Vega,axis=0)

    W2W3_Vega = (MAB_W2 - W2_Vega_to_AB) - (MAB_W3 - W3_Vega_to_AB)
    W2W3_Vega_cammean = np.mean(W2W3_Vega,axis=0)
    W2W3_Vega_cammax = np.max(W2W3_Vega,axis=0)
    W2W3_Vega_cammin = np.min(W2W3_Vega,axis=0)

    plot_LOSerr_vs_t(path,'iso_%sW1W2'%extra,time,W1W2_Vega_cammean,W1W2_Vega_cammin,
                     W1W2_Vega_cammax,runnames,iso_plotidx,oplot_vals=(0.5,0.8),
                     oplot_vals_style=('r:','r'),ylim=(-0.1,1.5),ylabel='W1-W2')
    plot_LOSerr_vs_t(path,'equal_%sW1W2'%extra,time,W1W2_Vega_cammean,W1W2_Vega_cammin,
                     W1W2_Vega_cammax,runnames,equal_plotidx,oplot_vals=(0.5,0.8),
                     oplot_vals_style=('r:','r'),ylim=(-0.1,1.5),ylabel='W1-W2')
    plot_LOSerr_vs_t(path,'uneq_%sW1W2'%extra,time,W1W2_Vega_cammean,W1W2_Vega_cammin,
                     W1W2_Vega_cammax,runnames,uneq_plotidx,oplot_vals=(0.5,0.8),
                     oplot_vals_style=('r:','r'),ylim=(-0.1,1.5),ylabel='W1-W2')

    plot_LOSerr_vs_t(path,'iso_%sW2W3'%extra,time,W2W3_Vega_cammean,W2W3_Vega_cammin,
                     W2W3_Vega_cammax,runnames,iso_plotidx,ylabel='W2-W3')
    plot_LOSerr_vs_t(path,'equal_%sW2W3'%extra,time,W2W3_Vega_cammean,W2W3_Vega_cammin,
                     W2W3_Vega_cammax,runnames,equal_plotidx,ylabel='W2-W3')
    plot_LOSerr_vs_t(path,'uneq_%sW2W3'%extra,time,W2W3_Vega_cammean,W2W3_Vega_cammin,
                     W2W3_Vega_cammax,runnames,uneq_plotidx,ylabel='W2-W3')

    print time.shape,sfr.shape,d_bh.shape
    print Llam.shape
    #for n in range(time.shape[1]):
    #    print len(np.where(time[:,n]>0)[0])
    #    print Llam[:,:,np.where(time[:,n]>0)[0],n].min()
    fig = plt.figure(figsize=(6,8))
    ax1 = fig.add_subplot(5,1,1)
    plt.yscale('log')
    plt.ylim(1.0e-3,100)
    for n in range(time.shape[1]):
        ix = np.where(time[:,n]>0)[0]
        ax1.plot(time[ix,n],sfr[ix,n])
    ax2 = fig.add_subplot(5,1,2)
    plt.yscale('log')
    plt.ylim(0.1,250)
    for n in range(time.shape[1]):
        ix = np.where(time[:,n]>0)[0]
        ax2.plot(time[ix,n],d_bh[ix,n])

    print W1W2_Vega.shape,W1W2_Vega_cammean.shape
    ax3 = fig.add_subplot(5,1,3)
    plt.xlim(0,6)
    for n in range(time.shape[1]):
        ix = np.where(time[:,n]>0)[0]
        #print "len time, W1W2 for %s:"%runnames[n],time[ix,n].shape,W1W2_Vega_cammean[ix,n].shape
        #print time[ix,n],W1W2_Vega_cammean[ix,n]
        ax3.plot(time[ix,n],W1W2_Vega_cammean[ix,n])
    #W1W2_Vega_cammin = np.min(W1W2_Vega,axis=0)
    #W1W2_Vega_cammax = np.max(W1W2_Vega,axis=0)

    ax4 = fig.add_subplot(5,1,4)
    for n in range(time.shape[1]):
        ix = np.where(time[:,n]>0)[0]
        ax4.plot(time[ix,n],fgas[ix,n])
    ax5 = fig.add_subplot(5,1,5)
    for n in range(time.shape[1]):
        ix = np.where(time[:,n]>0)[0]
        ax5.plot(time[ix,n],mdust[ix,n])

def gfs_plots(path='/oasis/projects/nsf/hvd115/lblecha/snyder_2013_data',
              file='agn_midir_Snyder_et_al_2013.fits'):

    global W1_Vega_to_AB, W2_Vega_to_AB, W3_Vega_to_AB, W4_Vega_to_AB

    ## key simulation times as defined in Greg's paper (pre-burst, peak SFR, peak AGN, post-burst)
    hi_times = [0.53, 0.69, 0.74, 0.99]
    low_times = [1.35, 1.63, 1.71, 1.79]

    f=pyfits.open("%s/%s"%(path,file))

    ix_Av = 11 ## index of Av param in data file
    ## DIAGNOSTICS_HIGHLYOBSCURED:
    ## shape (142,11,11,27) for time, simulation, viewing angle, param
    hdu_hiobsc = 1 
    ## simulations (MP OFF='default ISM', MP ON='alt ISM'):
    ##1.  AGNx1, MP OFF   (default ISM, multiphase model off, MW dust)
    ##2.  AGNx1, MP ON
    ##3.  AGNx10, MP OFF
    ##4.  AGNx10, MP ON
    ##5.  AGNx0, MP OFF
    ##6.  AGNx0, MP ON
    ##7.  AGNx1, LMC average dust (see paper), MP OFF
    ##8.  AGNx1, SMCbar dust (see paper), MP OFF
    ##9.  AGNx10, SMCbar dust, MP OFF
    ##10. AGNx10, SMCbar dust, MP ON
    ##11. AGNx1, SMCbar, MP ON
    simname_hiobsc = ['AGN1_MPoff_MWdust','AGN1_MPon_MWdust',
                      'AGN10_MPoff_MWdust','AGN10_MPon_MWdust',
                      'AGN0_MPoff_MWdust','AGN0_MPon_MWdust',
                      'AGN1_MPoff_LMCdust','AGN1_MPoff_SMCdust',
                      'AGN10_MPoff_SMCdust','AGN10_MPon_SMCdust',
                      'AGN1_MPon_SMCdust']
    hi_carr = ['b','c','b','c','gray','gray','g','r','r','m','m']
    hi_lsarr = ['-',':','-',':','','','-','-','-',':',':']
    hi_lwarr = [1,1,2,2,0,0,1,1,2,2,1]

    ## DIAGNOSTICS_LESSOBSCURED:
    ## shape (142,5,11,27) for time, simulation, viewing angle, param
    hdu_lowobsc = 3
    ## simulations:
    ##1.  AGNx1, ALT ISM, MW dust, less obscured merger
    ##2.  AGNx10, ALT ISM, MW dust
    ##3.  AGNx0, ALT ISM, MW dust
    ##4.  AGNx1, ALT ISM, SMCbar dust
    ##5.  AGNx10, ALT ISM, SMCbar dust
    simname_lowobsc = ['AGN1_MPon_MWdust','AGN10_MPon_MWdust',
                       'AGN0_MPon_MWdust','AGN1_MPon_SMCdust',
                       'AGN10_MPon_SMCdust']
    low_carr = ['c','c','gray','m','m']
    low_lwarr = [1,2,0,1,2]

    # for viewing angle: 
    #Description:
    #1.  average (arithmetic mean) of 7 viewing directions (where applicable)
    #2.  maximum of 7 viewing directions (where applicable)
    #3.  minimum of 7 viewing directions (where applicable)
    #4. through 10.  7 viewing angles evenly distributed in solid angle.
    #These are the same at all times for all simulations.
    #*Note*: Some quantities (i.e., SFR), are single-valued
    #11.  Measure of width of distribution with viewing angle. See Paper.
    #Measure defined S = MAD/0.67 = (median absolute deviation) / 0.67

    nsims_hiobsc = f[hdu_hiobsc].data.shape[1]
    nsims_lowobsc = f[hdu_lowobsc].data.shape[1]
    print nsims_hiobsc, nsims_lowobsc
    hiobsc_Av = f[hdu_hiobsc].data[:,:,:,ix_Av]
    lowobsc_Av = f[hdu_lowobsc].data[:,:,:,ix_Av]
    print "min/max hiobsc_Av: ",hiobsc_Av[hiobsc_Av==hiobsc_Av].min(),hiobsc_Av[hiobsc_Av==hiobsc_Av].max()
    print "min/max lowobsc_Av: ",lowobsc_Av[lowobsc_Av==lowobsc_Av].min(),lowobsc_Av[lowobsc_Av==lowobsc_Av].max()
    time_hiobsc = f[2].data
    time_lowobsc = f[4].data
    print time_hiobsc[0],time_hiobsc[30],time_hiobsc[60],time_hiobsc[90]
    print time_lowobsc[50],time_lowobsc[70],time_lowobsc[90],time_lowobsc[110],time_lowobsc[130]
    print time_lowobsc[102]
    #sys.exit()

    ## note: viewing angles are meaningless for these values, just use the first value (the 'average'):
    hiobsc_fAGN = f[hdu_hiobsc].data[:,:,0,0]
    #hiobsc_LAGN = f[hdu_hiobsc].data[:,:,0,26] ## something appears to be wrong with LAGN values. dont use.
    hiobsc_SFR = f[hdu_hiobsc].data[:,:,0,16]
    lowobsc_fAGN = f[hdu_lowobsc].data[:,:,0,0]
    #lowobsc_LAGN = f[hdu_lowobsc].data[:,:,0,26]
    lowobsc_SFR = f[hdu_lowobsc].data[:,:,0,16]

    ## now get the SEDs
    lam = f['LAMBDA_MICRONS'].data ##in microns
    ## shape (120,7,142,nsim) for wavelength, viewing angle, time, simulation
    hiobsc_Llam = f['L_LAMBDA_HIGHLYOBSCURED'].data
    lowobsc_Llam = f['L_LAMBDA_LESSOBSCURED'].data

    hiobsc_L_W1 = np.zeros(hiobsc_Llam.shape[1:])
    lowobsc_L_W1 = np.zeros(lowobsc_Llam.shape[1:])
    hiobsc_MAB_W1 = np.zeros(hiobsc_Llam.shape[1:])
    lowobsc_MAB_W1 = np.zeros(lowobsc_Llam.shape[1:])
    hiobsc_L_W2 = np.zeros(hiobsc_Llam.shape[1:])
    lowobsc_L_W2 = np.zeros(lowobsc_Llam.shape[1:])
    hiobsc_MAB_W2 = np.zeros(hiobsc_Llam.shape[1:])
    lowobsc_MAB_W2 = np.zeros(lowobsc_Llam.shape[1:])
    hiobsc_L_W3 = np.zeros(hiobsc_Llam.shape[1:])
    lowobsc_L_W3 = np.zeros(lowobsc_Llam.shape[1:])
    hiobsc_MAB_W3 = np.zeros(hiobsc_Llam.shape[1:])
    lowobsc_MAB_W3 = np.zeros(lowobsc_Llam.shape[1:])
    #print hiobsc_MAB_W1.shape,hiobsc_MAB_W2.shape,hiobsc_Llam.shape
    for icam in range(hiobsc_Llam.shape[1]):
        print "cam %d"%icam
        for jt in range(hiobsc_Llam.shape[2]):
            #if jt % 10 == 0: print "jt %d"%jt

            for ksim in range(hiobsc_Llam.shape[3]):
                #print "obtaining hiobsc W1 data for [%d,%d,%d]"%(icam,jt,ksim)
                hiobsc_L_W1[icam,jt,ksim],hiobsc_MAB_W1[icam,jt,ksim] = bb.Leff('WISE-W1',lam,hiobsc_Llam[:,icam,jt,ksim],sed_lum_units='W',return_ABMag=True)
                #print "obtaining hiobsc W2 data for [%d,%d,%d]"%(icam,jt,ksim)
                hiobsc_L_W2[icam,jt,ksim],hiobsc_MAB_W2[icam,jt,ksim] = bb.Leff('WISE-W2',lam,hiobsc_Llam[:,icam,jt,ksim],sed_lum_units='W',return_ABMag=True)
                #print "obtaining hiobsc W3 data for [%d,%d,%d]"%(icam,jt,ksim)
                hiobsc_L_W3[icam,jt,ksim],hiobsc_MAB_W3[icam,jt,ksim] = bb.Leff('WISE-W3',lam,hiobsc_Llam[:,icam,jt,ksim],sed_lum_units='W',return_ABMag=True)
                #hiobsc_MAB_W1[icam,jt,ksim] = bb.ABMag('WISE-W1',lam,hiobsc_Llam[:,icam,jt,ksim],sed_lum_units='W')

            for ksim in range(lowobsc_Llam.shape[3]):
                #print "obtaining lowobsc W1 data for [%d,%d,%d]"%(icam,jt,ksim)
                lowobsc_L_W1[icam,jt,ksim],lowobsc_MAB_W1[icam,jt,ksim] = bb.Leff('WISE-W1',lam,lowobsc_Llam[:,icam,jt,ksim],sed_lum_units='W',return_ABMag=True)
                #print "obtaining lowobsc W2 data for [%d,%d,%d]"%(icam,jt,ksim)
                lowobsc_L_W2[icam,jt,ksim],lowobsc_MAB_W2[icam,jt,ksim] = bb.Leff('WISE-W2',lam,lowobsc_Llam[:,icam,jt,ksim],sed_lum_units='W',return_ABMag=True)
                #print "obtaining lowobsc W3 data for [%d,%d,%d]"%(icam,jt,ksim)
                lowobsc_L_W3[icam,jt,ksim],lowobsc_MAB_W3[icam,jt,ksim] = bb.Leff('WISE-W3',lam,lowobsc_Llam[:,icam,jt,ksim],sed_lum_units='W',return_ABMag=True)
                #lowobsc_MAB_W1[icam,jt,ksim] = bb.ABMag('WISE-W1',lam,lowobsc_Llam[:,icam,jt,ksim],sed_lum_units='W')
        
        #print "hiobsc_MAB_W2:",hiobsc_MAB_W2[icam,:,:]
        #print "hiobsc_MAB_W3:",hiobsc_MAB_W2[icam,:,:]
        #tmp_hiobsc_W2W3_Vega = (hiobsc_MAB_W2[icam,:,:] - W2_Vega_to_AB) - (hiobsc_MAB_W3[icam,:,:] - W3_Vega_to_AB)
        #print "tmp_hiobsc_W2W3_Vega:",tmp_hiobsc_W2W3_Vega

    hiobsc_W1W2_Vega = (hiobsc_MAB_W1 - W1_Vega_to_AB) - (hiobsc_MAB_W2 - W2_Vega_to_AB)
    hiobsc_W1W2_Vega_cammean = np.mean(hiobsc_W1W2_Vega,axis=0)
    hiobsc_W1W2_Vega_cammax = np.max(hiobsc_W1W2_Vega,axis=0)
    hiobsc_W1W2_Vega_cammin = np.min(hiobsc_W1W2_Vega,axis=0)
    lowobsc_W1W2_Vega = (lowobsc_MAB_W1 - W1_Vega_to_AB) - (lowobsc_MAB_W2 - W2_Vega_to_AB)
    lowobsc_W1W2_Vega_cammean = np.mean(lowobsc_W1W2_Vega,axis=0)
    lowobsc_W1W2_Vega_cammax = np.max(lowobsc_W1W2_Vega,axis=0)
    lowobsc_W1W2_Vega_cammin = np.min(lowobsc_W1W2_Vega,axis=0)

    hiobsc_W2W3_Vega = (hiobsc_MAB_W2 - W2_Vega_to_AB) - (hiobsc_MAB_W3 - W3_Vega_to_AB)
    hiobsc_W2W3_Vega_cammean = np.mean(hiobsc_W2W3_Vega,axis=0)
    hiobsc_W2W3_Vega_cammax = np.max(hiobsc_W2W3_Vega,axis=0)
    hiobsc_W2W3_Vega_cammin = np.min(hiobsc_W2W3_Vega,axis=0)
    lowobsc_W2W3_Vega = (lowobsc_MAB_W2 - W2_Vega_to_AB) - (lowobsc_MAB_W3 - W3_Vega_to_AB)
    lowobsc_W2W3_Vega_cammean = np.mean(lowobsc_W2W3_Vega,axis=0)
    lowobsc_W2W3_Vega_cammax = np.max(lowobsc_W2W3_Vega,axis=0)
    lowobsc_W2W3_Vega_cammin = np.min(lowobsc_W2W3_Vega,axis=0)

    #print 'min/max hiobsc_L_W2:',hiobsc_L_W2.min(),hiobsc_L_W2.max()
    #print 'min/max hiobsc_L_W3:',hiobsc_L_W3.min(),hiobsc_L_W3.max()
    #print 'hiobsc_MAB_W2:',hiobsc_MAB_W2
    #print 'hiobsc_MAB_W3:',hiobsc_MAB_W3
    print 'min/max hiobsc_W2W3_Vega:',hiobsc_W2W3_Vega[hiobsc_W2W3_Vega==hiobsc_W2W3_Vega].min(),hiobsc_W2W3_Vega[hiobsc_W2W3_Vega==hiobsc_W2W3_Vega].max()
    #print 'min/max lowobsc_L_W2:',lowobsc_L_W2.min(),lowobsc_L_W2.max()
    #print 'min/max lowobsc_L_W3:',lowobsc_L_W3.min(),lowobsc_L_W3.max()
    #print 'min/max lowobsc_MAB_W2:',lowobsc_MAB_W2.min(),lowobsc_MAB_W2.max()
    #print 'min/max lowobsc_MAB_W3:',lowobsc_MAB_W3.min(),lowobsc_MAB_W3.max()
    print 'min/max lowobsc_W2W3_Vega:',lowobsc_W2W3_Vega[lowobsc_W2W3_Vega==lowobsc_W2W3_Vega].min(),lowobsc_W2W3_Vega[lowobsc_W2W3_Vega==lowobsc_W2W3_Vega].max()
    #print 'min/max lowobsc_W2W3_Vega:',lowobsc_W2W3_Vega.min(),lowobsc_W2W3_Vega.max()
    #sys.exit()

    #xlim_lowobsc=(0.6,2.1)
    #xlim_hiobsc=(0,1.5)
    xlim_lowobsc=(0.7,2.1)
    xlim_hiobsc=(0.0,1.4)
    plt.clf()
    plt.close("all")

    low_plotidx = [simname_lowobsc.index(name) for name in ['AGN0_MPon_MWdust','AGN1_MPon_MWdust','AGN1_MPon_SMCdust','AGN10_MPon_MWdust','AGN10_MPon_SMCdust']]

    print "making lowobsc plots."

    matplotlib.rcParams.update({'font.size': 11})

    plot_LOSerr_vs_t(path,'lowobsc_Av',time_lowobsc,lowobsc_Av[:,:,0],
                     lowobsc_Av[:,:,2],lowobsc_Av[:,:,1],simname_lowobsc,
                     low_plotidx, oplot_times=(low_times[1:3]),
                     xlim=xlim_lowobsc,ylim=(-0.03,3.5),ylabel='Av')

    plot_LOSerr_vs_t(path,'lowobsc_W1W2',time_lowobsc,
                     lowobsc_W1W2_Vega_cammean,lowobsc_W1W2_Vega_cammin,
                     lowobsc_W1W2_Vega_cammax,simname_lowobsc,
                     low_plotidx, oplot_times=low_times[1:3],
                     oplot_vals=(0.5,0.8),oplot_vals_style=('r:','r'),
                     xlim=xlim_lowobsc,ylim=(-0.1,2.5),ylabel='W1-W2')

    plot_LOSerr_vs_t(path,'lowobsc_W2W3',time_lowobsc,
                     lowobsc_W2W3_Vega_cammean,lowobsc_W2W3_Vega_cammin,
                     lowobsc_W2W3_Vega_cammax,simname_lowobsc,
                     low_plotidx, oplot_times=low_times[1:3],
                     #oplot_vals=(0.5,0.8),oplot_vals_style=('r:','r'),
                     xlim=xlim_lowobsc,ylim=(),ylabel='W2-W3')

    w1w2_scatter_plot(path,'lowobsc_fAGN_vs_W1W2',lowobsc_fAGN,
                      lowobsc_W1W2_Vega,simname_lowobsc,low_plotidx,
                      xlim=(-2.5,0),xlabel='log(LAGN/Lbol)')

    w1w2_scatter_plot(path,'lowobsc_Av_vs_W1W2',lowobsc_Av[:,:,3:10],
                      lowobsc_W1W2_Vega,simname_lowobsc,low_plotidx,
                      xlim=(-0.03,3.5),xlabel='Av')

    ### plot lowobsc fAGN and SFR for all sims
    fig = plt.figure(figsize=(6,8))
    ax1 = fig.add_subplot(2,1,1)
    plt.xlabel('time [Gyr]')
    plt.ylabel('log(LAGN/Lbol)')
    plt.xlim(xlim_lowobsc)
    ylim = (-2.5,0.1)
    plt.ylim(ylim)
    ax1.plot([low_times[1],low_times[1]],ylim,'k',ls='--')
    ax1.plot([low_times[2],low_times[2]],ylim,'k',ls='--')
    for nsim in range(nsims_lowobsc):
        if simname_lowobsc[nsim]=='AGN10_MPon_SMCdust':
            ax1.plot(time_lowobsc[time_lowobsc<1.95],
                     lowobsc_fAGN[time_lowobsc<1.95,nsim],
                     color=low_carr[nsim],lw=low_lwarr[nsim])  
        else:
            ax1.plot(time_lowobsc,lowobsc_fAGN[:,nsim],
                     color=low_carr[nsim],lw=low_lwarr[nsim])
    ax2 = fig.add_subplot(2,1,2)
    plt.xlabel('time [Gyr]')
    plt.ylabel('log(SFR) [Msun/yr]')
    plt.xlim(xlim_lowobsc)
    ylim = (-2,4)
    plt.ylim(ylim)
    ax2.plot([low_times[1],low_times[1]],ylim,'k',ls='--')
    ax2.plot([low_times[2],low_times[2]],ylim,'k',ls='--')
    for nsim in range(nsims_lowobsc):
        ax2.plot(time_lowobsc,lowobsc_SFR[:,nsim],
                 color=low_carr[nsim],lw=low_lwarr[nsim])
    fig.subplots_adjust(hspace=0.4,left=0.2)
    fig.savefig('%s/lowobsc_fAGN_SFR.eps'%path)


    hi_plotidx = [simname_hiobsc.index(name) for name in ['AGN0_MPon_MWdust','AGN1_MPon_MWdust','AGN1_MPon_SMCdust','AGN10_MPon_MWdust','AGN10_MPon_SMCdust']]

    print "making hiobsc plots."

    plot_LOSerr_vs_t(path,'hiobsc_Av',time_hiobsc, hiobsc_Av[:,:,0],
                     hiobsc_Av[:,:,2], hiobsc_Av[:,:,1],simname_hiobsc,
                     hi_plotidx, oplot_times=hi_times[1:3],
                     xlim=xlim_hiobsc,ylim=(-0.03,3.5),ylabel='Av')

    plot_LOSerr_vs_t(path,'hiobsc_W1W2',time_hiobsc,
                     hiobsc_W1W2_Vega_cammean,hiobsc_W1W2_Vega_cammin,
                     hiobsc_W1W2_Vega_cammax,simname_hiobsc,
                     hi_plotidx, oplot_times=hi_times[1:3],
                     oplot_vals=(0.5,0.8),oplot_vals_style=('r:','r'),
                     xlim=xlim_hiobsc,ylim=(-0.1,2.5),ylabel='W1-W2')

    plot_LOSerr_vs_t(path,'hiobsc_W2W3',time_hiobsc,
                     hiobsc_W2W3_Vega_cammean,hiobsc_W2W3_Vega_cammin,
                     hiobsc_W2W3_Vega_cammax,simname_hiobsc,
                     hi_plotidx, oplot_times=hi_times[1:3],
                     #oplot_vals=(0.5,0.8),oplot_vals_style=('r:','r'),
                     xlim=xlim_hiobsc,ylim=(),ylabel='W2-W3')

    w1w2_scatter_plot(path,'hiobsc_fAGN_vs_W1W2',hiobsc_fAGN,
                      hiobsc_W1W2_Vega,simname_hiobsc,hi_plotidx,
                      xlim=(-2.5,0),xlabel='log(LAGN/Lbol)')

    w1w2_scatter_plot(path,'hiobsc_Av_vs_W1W2',hiobsc_Av[:,:,3:10],
                      hiobsc_W1W2_Vega,simname_hiobsc,hi_plotidx,
                      xlim=(-0.03,3.5),xlabel='Av')


    hi_nomp_plotidx = [simname_hiobsc.index(name) for name in ['AGN0_MPoff_MWdust','AGN1_MPoff_MWdust','AGN1_MPoff_LMCdust','AGN1_MPoff_SMCdust','AGN10_MPoff_MWdust','AGN10_MPoff_SMCdust']]

    plot_LOSerr_vs_t(path,'hiobsc_MPoff_Av',time_hiobsc,hiobsc_Av[:,:,0],
                     hiobsc_Av[:,:,2],hiobsc_Av[:,:,1],simname_hiobsc,
                     hi_nomp_plotidx, oplot_times=hi_times[1:3],
                     xlim=xlim_hiobsc,ylim=(-0.03,3.5),ylabel='Av')

    plot_LOSerr_vs_t(path,'hiobsc_MPoff_W1W2',time_hiobsc,
                     hiobsc_W1W2_Vega_cammean,hiobsc_W1W2_Vega_cammin,
                     hiobsc_W1W2_Vega_cammax,simname_hiobsc,
                     hi_nomp_plotidx, oplot_times=hi_times[1:3],
                     oplot_vals=(0.5,0.8),oplot_vals_style=('r:','r'),
                     xlim=xlim_hiobsc,ylim=(-0.1,2.5),ylabel='W1-W2')

    plot_LOSerr_vs_t(path,'hiobsc_MPoff_W2W3',time_hiobsc,
                     hiobsc_W2W3_Vega_cammean,hiobsc_W2W3_Vega_cammin,
                     hiobsc_W2W3_Vega_cammax,simname_hiobsc,
                     hi_nomp_plotidx, oplot_times=hi_times[1:3],
                     #oplot_vals=(0.5,0.8),oplot_vals_style=('r:','r'),
                     xlim=xlim_hiobsc,ylim=(),ylabel='W2-W3')

    w1w2_scatter_plot(path,'hiobsc_MPoff_fAGN_vs_W1W2',hiobsc_fAGN,
                      hiobsc_W1W2_Vega,simname_hiobsc,hi_nomp_plotidx,
                      xlim=(-2.5,0),xlabel='log(LAGN/Lbol)')

    w1w2_scatter_plot(path,'hiobsc_MPoff_Av_vs_W1W2',hiobsc_Av[:,:,3:10],
                      hiobsc_W1W2_Vega,simname_hiobsc,hi_nomp_plotidx,
                      xlim=(-0.03,3.5),xlabel='Av')


    ### plot hiobsc fAGN and SFR for all sims
    fig = plt.figure(figsize=(6,8))
    ax1 = fig.add_subplot(2,1,1)
    plt.xlabel('time [Gyr]')
    plt.ylabel('log(LAGN/Lbol)')
    plt.xlim(xlim_hiobsc)
    ylim = (-3.5,0.1)
    plt.ylim(ylim)
    ax1.plot([hi_times[1],hi_times[1]],ylim,'k',ls='--')
    ax1.plot([hi_times[2],hi_times[2]],ylim,'k',ls='--')
    for nsim in range(nsims_hiobsc):
        ax1.plot(time_hiobsc,hiobsc_fAGN[:,nsim],color=hi_carr[nsim],
                 lw=hi_lwarr[nsim],ls=hi_lsarr[nsim])
    ax2 = fig.add_subplot(2,1,2)
    plt.xlabel('time [Gyr]')
    plt.ylabel('log(SFR) [Msun/yr]')
    plt.xlim(xlim_hiobsc)
    ylim = (-2,4)
    plt.ylim(ylim)
    ax2.plot([hi_times[1],hi_times[1]],ylim,'k',ls='--')
    ax2.plot([hi_times[2],hi_times[2]],ylim,'k',ls='--')
    for nsim in range(nsims_hiobsc):
        ax2.plot(time_hiobsc,hiobsc_SFR[:,nsim],color=hi_carr[nsim],lw=hi_lwarr[nsim],ls=hi_lsarr[nsim])
    fig.subplots_adjust(hspace=0.4,left=0.2)
    fig.savefig('%s/hiobsc_fAGN_SFR.eps'%path)




def plot_LOSerr_vs_t(path,plotname,time,meanval,minval,maxval,simnames,
                     plotidx,oplot_times=(),oplot_vals=(),oplot_vals_style=(),
                     xlim=(),ylim=(),ylabel=''):
    
    nsims = len(plotidx)
    if nsims<1 or nsims>6:
        print "Error: plot_LOSerr_vs_t() can plot data for 1-6 simulations."
        print simnames[plotidx]
        raise AssertionError

    if len(xlim)!=2: xlim = (time.min(),time.max())
    if len(ylim)!=2: ylim = (0.95*minval[minval==minval].min(),1.05*maxval[maxval==maxval].max())


    fig = plt.figure(figsize=(8,8))
    for n in range(nsims):
        if nsims == 5:
            ax = fig.add_subplot(3,2,1) if n==0 else fig.add_subplot(3,2,n+2)
        else:
            ax = fig.add_subplot(3,2,n+1)
        
        plt.xlabel('time [Gyr]')
        plt.ylabel(ylabel)
        plt.xlim(xlim)
        plt.ylim(ylim)

        if len(oplot_times)>0:
            for t in oplot_times: ax.plot([t,t],ylim,'k',ls='--')

        if time.shape == meanval.shape:
            ix = np.where(time[:,plotidx[n]]>0)[0]
            simtime = time[ix,plotidx[n]]
            #print "len(simtime)=%d for %s"%(len(simtime),simnames[plotidx[n]])
        elif time.size == meanval[:,0].size:
            ix = np.arange(time.size)
            simtime = time

        ax.errorbar(simtime, meanval[ix,plotidx[n]], 
                    yerr=(meanval[ix,plotidx[n]]-minval[ix,plotidx[n]],
                          maxval[ix,plotidx[n]]-meanval[ix,plotidx[n]]), 
                    color='c')
        ax.plot(simtime, meanval[ix,plotidx[n]], 'k', linewidth=1.5)
        
        if len(oplot_vals)>0:
            if len(oplot_vals_style)!=len(oplot_vals):
                oplot_vals_style=('r',)*len(oplot_vals)
            for i,v in enumerate(oplot_vals): 
                ax.plot(xlim,[v,v],oplot_vals_style[i])
        plt.title(simnames[plotidx[n]])

    fig.subplots_adjust(hspace=0.5,wspace=0.3)
    fig.savefig('%s/%s.eps'%(path,plotname))

def w1w2_scatter_plot(path,plotname,val,w1w2,simnames,
                      plotidx,oplot_val=(),oplot_w1w2=(),
                      xlim=(),ylim=(-0.1,2.5),xlabel='',ylabel='W1-W2'):
    
    nsims = len(plotidx)
    if nsims<1 or nsims>6:
        print "Error: w1w2_scatter_plot() can plot data for 1-6 simulations."
        print simnames[plotidx]
        raise AssertionError

    if len(xlim)!=2: xlim = (0.95*val.min(),1.05*val.max())
    if len(ylim)!=2: ylim = (0.95*w1w2.min(),1.05*w1w2.max())
    
    fig = plt.figure(figsize=(8,8))
    for n in range(nsims):
        if nsims == 5:
            ax = fig.add_subplot(3,2,1) if n==0 else fig.add_subplot(3,2,n+2)
        else:
            ax = fig.add_subplot(3,2,n+1)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xlim(xlim)
        plt.ylim(ylim)

        if len(oplot_val)>0:
            for v in oplot_val: ax.plot([v,v],ylim,'k',ls='--')

        if len(val.shape)==3:
            for icam in range(7):
                ax.plot(val[:,plotidx[n],icam],
                        w1w2[icam,:,plotidx[n]],'ko',ms=0.5)
        elif len(val.shape)==2:
            for icam in range(7):
                ax.plot(val[:,plotidx[n]],w1w2[icam,:,plotidx[n]],'ko',ms=0.5)
        else: raise AssertionError
        plt.title(simnames[plotidx[n]])
    fig.subplots_adjust(hspace=0.6,wspace=0.3)
    fig.savefig('%s/%s.eps'%(path,plotname))


def oplotvals_vs_t(path,plotname,simnames,plotidx,time,val1,val2,val3=(),
                   oplot_times=(),xlim=(),ylim1=(),ylim2=(),ylim3=(),
                   ylabel1='',ylabel2='',ylabel3='',
                   ylog1=True,ylog2=True,ylog3=True,
                   carr=['k','b','m','g','c','r']):
    
    nsims = len(plotidx)
    if nsims<1 or nsims>6:
        print "Error: oplotvals_vs_t() can plot data for 1-6 simulations."
        print [simnames[ix] for ix in plotidx]
        raise AssertionError

    if len(xlim)!=2: xlim = (time.min(),time.max())
    if len(ylim1)!=2: ylim1 = (0.95*val1.min(),1.05*val1.max())
    if len(ylim2)!=2: ylim2 = (0.95*val2.min(),1.05*val2.max())
    if len(val3)>0 and len(ylim3)!=2: ylim3 = (0.95*val3.min(),1.05*val3.max())
    nplots=3 if len(val3)>0 else 2
    
    fig = plt.figure(figsize=(8,8))

    ax1 = fig.add_subplot(nplots,1,1)
    plt.xlabel('time [Gyr]')
    plt.ylabel(ylabel1)
    plt.xlim(xlim)
    plt.ylim(ylim1)
    if ylog1: plt.yscale('log')

    if len(oplot_times)>0:
        for t in oplot_times: ax1.plot([t,t],ylim1,'k',ls='--')
    for n in range(nsims):
        ix = np.where(time[:,plotidx[n]]>0)[0]
        ax1.plot(time[ix,plotidx[n]], val1[ix,plotidx[n]], 'o-',ms=0.5,
                 color=carr[n],label=simnames[plotidx[n]])
    ax1.legend(fontsize=10)

    ax2 = fig.add_subplot(nplots,1,2)
    plt.xlabel('time [Gyr]')
    plt.ylabel(ylabel2)
    plt.xlim(xlim)
    plt.ylim(ylim2)
    if ylog2: plt.yscale('log')

    if len(oplot_times)>0:
        for t in oplot_times: ax2.plot([t,t],ylim2,'k',ls='--')
    for n in range(nsims):
        ix = np.where(time[:,plotidx[n]]>0)[0]
        ax2.plot(time[ix,plotidx[n]], val2[ix,plotidx[n]], 'o-',ms=0.6,
                 color=carr[n],label=simnames[plotidx[n]])
        ax2.legend(fontsize=10)

    if nplots == 3:
        ax3 = fig.add_subplot(nplots,1,3)
        plt.xlabel('time [Gyr]')
        plt.ylabel(ylabel3)
        plt.xlim(xlim)
        plt.ylim(ylim3)
        if ylog3: plt.yscale('log')
        
        if len(oplot_times)>0:
            for t in oplot_times: ax3.plot([t,t],ylim3,'k',ls='--')
        for n in range(nsims):
            #print val3[ix,plotidx[n]].min(),val3[ix,plotidx[n]].max()
            ix = np.where(time[:,plotidx[n]]>0)[0]
            ax3.plot(time[ix,plotidx[n]], val3[ix,plotidx[n]], 'o-',ms=0.6,
                     color=carr[n],label=simnames[plotidx[n]])
            ax3.legend(fontsize=10)
        
    fig.subplots_adjust(hspace=0.5,wspace=0.3)
    fig.savefig('%s/%s.eps'%(path,plotname))
