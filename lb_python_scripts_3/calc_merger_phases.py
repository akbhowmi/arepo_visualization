import matplotlib
matplotlib.use('Agg')
import glob, pylab, numpy, h5py, re, sys, os
from copy import copy
import astro_constants as ac
from pylab import *
import readsnapbhs


def rotate(v,theta,phi):
    
    if len(v) != 3:
        print "Error: function 'rotate' requires a 3-d vector."
        print "v = ",v
        sys.exit()
        
    ## perform a counterclockwise rotation by phi about the z axis
    ## followed by a clockwise rotation by theta about the y axis
    xnew = v[0]*np.cos(theta)*np.cos(phi) - v[2]*np.sin(theta) + v[1]*np.cos(theta)*np.sin(phi)
    ynew = v[1]*np.cos(phi) - v[0]*np.sin(phi)
    znew = v[2]*np.cos(theta) + v[0]*np.cos(phi)*np.sin(theta) + v[1]*np.sin(theta)*np.sin(phi)

    return np.array([xnew,ynew,znew])
                                                    

def rotate_arr(v,theta,phi):
    
    if v.shape[0] != 3:
        print "Error: function 'rotate' requires a 3xN array."
        print "v has shape ", v.shape
        sys.exit()

    if theta.size != phi.size:
        print "Error: array size mismatch b/t theta & phi: ",theta.size,phi.size
        sys.exit()
    
    ## perform a counterclockwise rotation by phi about the z axis
    ## followed by a clockwise rotation by theta about the y axis
    xnew = np.array([ v[0,:]*np.cos(theta[i])*np.cos(phi[i]) - v[2,:]*np.sin(theta[i]) +
                      v[1,:]*np.cos(theta[i])*np.sin(phi[i]) for i in range(len(theta)) ])
    ynew = np.array([ v[1,:]*np.cos(phi[i]) - v[0,:]*np.sin(phi[i])
                      for  i in range(len(theta)) ])
    znew = np.array([ v[2,:]*np.cos(theta[i]) + v[0,:]*np.cos(phi[i])*np.sin(theta[i]) +
                      v[1,:]*np.sin(theta[i])*np.sin(phi[i]) for i in range(len(theta)) ])

    ## returned array has shape (3, ncam, narr):
    return np.array([xnew, ynew, znew])
                                                    

def calc_lbol_edd(mbh):
    return mbh*ac.MSUN * 4*np.pi*ac.G*ac.MP*ac.C / (ac.THOMSON)


def plot(path,tmax_postmrg=0.25,ncam=10,plot_proj_sep=True):

    print "Warning: function 'plot' is deprecated. Please use 'multiplot' instead."
    
    ### load bh data ###
    allbhs = readsnapbhs.all_snapshot_bhs(path,tmax_postmrg=tmax_postmrg)
    
    has2bh = (allbhs.nbh==2)

    lagn = copy(allbhs.lbol[0,:])
    lagn[has2bh] = allbhs.lbol[:,has2bh].sum(axis=0)

    bh_mass = copy(allbhs.mass[0,:])
    bh_mass[has2bh] = allbhs.mass[:,has2bh].sum(axis=0)
    ledd = calc_lbol_edd(bh_mass)
    fedd = np.zeros(lagn.size)
    fedd[ledd>0] = lagn[ledd>0]/ledd[ledd>0]

    lagn_ratio = np.zeros(allbhs.snaps.size)
    lagn_ratio[has2bh] = allbhs.lbol[0,has2bh]/allbhs.lbol[1,has2bh]
    #lagn_ratio[lagn_ratio>1.0] = 1.0/lagn_ratio[lagn_ratio>1.0]
    lg_lagn_ratio = np.log10(lagn_ratio)
    
    fedd1 = allbhs.lbol[0,:]/allbhs.lbol_edd[0,:]
    fedd2 = np.zeros(lagn.size)
    fedd2[has2bh] = allbhs.lbol[1,has2bh]/allbhs.lbol_edd[1,has2bh]

    dt = np.zeros(allbhs.time.size)
    dt[1:] = allbhs.time[1:] - allbhs.time[:-1]
    dt[0] = np.minimum(dt[1],allbhs.time[0])    
    ttot = dt.sum()
    dt_tiled = np.tile(dt, (ncam,1))
    
    bh_3d_sep = np.repeat(np.nan,allbhs.snaps.size)
    bh_3d_sep[has2bh] = np.sqrt( np.sum(((allbhs.pos[1,j,has2bh]-allbhs.pos[0,j,has2bh])**2 
                                        for j in range(3)), axis=0) )
    bh_3d_sep[bh_3d_sep!=bh_3d_sep] = -1.0

    ### get random sight lines ###
    phi = np.random.rand(ncam) * 2*np.pi
    theta = np.random.rand(ncam) * np.pi
    rotated_pos1 = rotate_arr(allbhs.pos[0,:,:], theta, phi)
    rotated_pos2 = rotate_arr(allbhs.pos[1,:,:], theta, phi)
    bh_proj_sep = np.sqrt( np.sum((rotated_pos2[:2,:,:]-rotated_pos1[:2,:,:])**2,axis=0) )
    bh_proj_sep[bh_proj_sep!=bh_proj_sep] = -1.0
    assert np.array([np.allclose(bh_3d_sep[has2bh], np.sqrt(np.sum((rotated_pos2[:,j,has2bh]-rotated_pos1[:,j,has2bh])**2,axis=0)) ) for j in range(ncam)]).all(), \
        "Something went wrong in calculating rotated BH positions! Mismatch in 3D sep calculation."
    #assert np.array_equal(bh_proj_sep[:,has2bh],bh_proj_sep[:,has2bh]), \
    assert np.array_equal(bh_proj_sep,bh_proj_sep), \
        "Something went wrong in calculating rotated BH positions! nan's in bh_proj_sep."
    for i in range(allbhs.snaps.size):
        print "%g %g %g"%(bh_3d_sep[i],bh_proj_sep[:,i].min(),bh_proj_sep[:,i].max())
    
    print "min/max/med lagn:",lagn.min(),lagn.max(),np.median(lagn)
    print "min/max/med bh1_mass:",allbhs.mass[0,:].min(),allbhs.mass[0,:].max(),np.median(allbhs.mass[0,:])
    print "min/max/med bh_mass:",bh_mass.min(),bh_mass.max(),np.median(bh_mass)
    print "min/max/med ledd:",ledd.min(),ledd.max(),np.median(ledd)
    print "min/max/med fedd:",fedd.min(),fedd.max(),np.median(fedd)
    print "min/max/med abs(log(lagn_ratio)):",np.abs(lg_lagn_ratio[lagn_ratio>0]).min(),np.abs(lg_lagn_ratio).max(),np.median(np.abs(lg_lagn_ratio[lagn_ratio>0]))

    print "ttot = ",dt.sum()

    print "\nEdd-ratio AGN definitions:"
    for fedd_agn_lim in [0.01, 0.05, 0.1]:
        ttot_agn = dt[fedd>fedd_agn_lim].sum()
        ttot_dualagn = dt[(fedd1>fedd_agn_lim)&(fedd2>fedd_agn_lim)].sum()
        ttot_dualagn_ratiocut = dt[(fedd1>fedd_agn_lim)&(fedd2>fedd_agn_lim)&
                                   (np.abs(lg_lagn_ratio)<1)].sum()

        print "\n  fedd_agn_lim = %g"%fedd_agn_lim
        print "  ttot_agn = ",ttot_agn
        print "  ttot_dualagn = ",ttot_dualagn
        print "  ttot_dualagn_ratiocut = ",ttot_dualagn_ratiocut

    print "\nLbol AGN definitions:"
    for lbol_agn_lim in [1.0e43, 1.0e44, 1.0e45]:
        ttot_agn = dt[lagn>lbol_agn_lim].sum()
        ttot_dualagn = dt[(allbhs.lbol[0,:]>lbol_agn_lim)&
                          (allbhs.lbol[1,:]>lbol_agn_lim)].sum()
        ttot_dualagn_ratiocut = dt[(allbhs.lbol[0,:]>lbol_agn_lim)&
                                   (allbhs.lbol[1,:]>lbol_agn_lim)&
                                   (np.abs(lg_lagn_ratio)<1)].sum()

        print "\n  lbol_agn_lim = %g"%lbol_agn_lim
        print "  ttot_agn = ",ttot_agn
        print "  ttot_dualagn = ",ttot_dualagn
        print "  ttot_dualagn_ratiocut = ",ttot_dualagn_ratiocut


    dsep=10
    sep_bins = np.arange(-2*dsep,80,dsep)
    #sep_bins = np.arange(0,10,dsep)

    fedd_agn_min = 0.01                                                              
    #fedd_agn_min = 0.1
    lbol_agn_min = 1.0e43

    #dt_3dsep_hist,tmp = np.histogram(bh_3d_sep,sep_bins,weights=dt)
    dt_3dsep_hist,tmp = np.histogram(bh_3d_sep,sep_bins,weights=dt/ttot)
    dt_feddagn_3dsep_hist,tmp = np.histogram(bh_3d_sep[fedd>fedd_agn_min],
                                             sep_bins,weights=dt[fedd>fedd_agn_min]/ttot)
    dt_lbolagn_3dsep_hist,tmp = np.histogram(bh_3d_sep[lagn>lbol_agn_min],
                                             sep_bins,weights=dt[lagn>lbol_agn_min]/ttot)
    dt_dual_feddagn_3dsep_hist,tmp = np.histogram(bh_3d_sep[(fedd1>fedd_agn_min)&(fedd2>fedd_agn_min)],
                                                  sep_bins,
                                                  weights=dt[(fedd1>fedd_agn_min)&(fedd2>fedd_agn_min)]/ttot)
    dt_dual_lbolagn_3dsep_hist,tmp = np.histogram(bh_3d_sep[(allbhs.lbol[0,:]>lbol_agn_min)&
                                                            (allbhs.lbol[1,:]>lbol_agn_min)],
                                                  sep_bins,
                                                  weights=dt[(allbhs.lbol[0,:]>lbol_agn_min)&
                                                             (allbhs.lbol[1,:]>lbol_agn_min)]/ttot)

    print dt_3dsep_hist.shape, sep_bins[:-1].shape,tmp.shape
    print dt[bh_3d_sep>0].sum()/ttot
    print dt[bh_3d_sep<=0].sum()/ttot
    #for i in range(len(dt_3dsep_hist)): print sep_bins[:-1][i], dt_3dsep_hist[i]

    fedd_tiled = np.tile(fedd,(ncam,1)) if plot_proj_sep else fedd
    lagn_tiled = np.tile(lagn,(ncam,1)) if plot_proj_sep else lagn
    fedd_maskarr = [(fedd_tiled>=lim) for lim in [0.01,0.05,0.1]]
    fedd_titles = ['fEdd>%g'%(lim) for lim in [0.01,0.05,0.1]]
    lagn_maskarr = [(lagn_tiled>=lim) for lim in [1.0e43, 1.0e44, 1.0e45]]
    lagn_titles = ['Lagn>%g'%(lim)  for lim in [1.0e43, 1.0e44, 1.0e45]]

    maskarr_list = fedd_maskarr + lagn_maskarr
    titles = fedd_titles + lagn_titles

    print fedd_tiled.shape, lagn_tiled.shape
    print len(fedd_maskarr), len(lagn_maskarr), len(lagn_maskarr[0])
    print np.array(maskarr_list).shape
    print bh_proj_sep.shape, dt_tiled.shape
    
    plt.ioff()
    matplotlib.rcParams.update({'font.size':10})
    
    tothist_kwargs = dict(range=(sep_bins.min(),sep_bins.max()),bins=len(sep_bins)-1,
                          histtype='step')
    
    fig = plt.figure(figsize=(8,6))
    for i in range(6):
        ax = fig.add_subplot(231+i)
        ax.set_ylim(0.0,0.2)
        ax.set_xlabel('3D sep [kpc]') if i in [6,7,8] else ax.set_xlabel('')
        ax.set_ylabel(r'$\Delta$t/t$_{\rm tot}$') if  i in [0,3,6] else ax.set_ylabel('')
        if plot_proj_sep:
            ax.hist(bh_proj_sep.flatten(), weights=dt_tiled.flatten()/(1.0*ncam)/ttot,
                    color='k',**tothist_kwargs)
            ax.hist(bh_proj_sep[maskarr_list[i]], weights=dt_tiled[maskarr_list[i]]/(1.0*ncam)/ttot, 
                    color='m', **tothist_kwargs)
        else:
            ax.hist(bh_3d_sep, weights=dt/ttot, color='k',**tothist_kwargs)
            ax.hist(bh_3d_sep[maskarr_list[i]], weights=dt[maskarr_list[i]]/ttot, 
                    color='m', **tothist_kwargs)
        #data = bh_3d_sep[maskarr_list[i]]
        #wts = dt[maskarr_list[i]]/ttot 
        #n,bins,patches=ax.hist(data, weights=wts, **stackhist_kwargs)
        ax.legend(fontsize=9)
        ax.set_title(titles[i],fontsize=10)
    suptit= "%s (Ncam=%d)"%(path,ncam) if plot_proj_sep else path
    fig.suptitle(suptit)
    fig.subplots_adjust(wspace=0.3,hspace=0.3,left=0.08,right=0.94,bottom=0.06,top=0.92)
    extra = '_projsep_n%d'%ncam if plot_proj_sep else '_3dsep'
    fig.savefig(path+'/merger_agn_phases%s.pdf'%extra,format='pdf')


def makebar(ax,hist,nbins=6,color='c',ecolor='k',width=0.3,xoffset=0.0,
            val='median', errtype='full', verbose=False, nsim=1,
            binpoints_only=False, reverse=False, xmin=0.0, xgap_index=-1,
            elinewidth=1.5,capthick=1.5):

    if hist.shape[-1] != nbins:
        print "Error: last axis of hist array must have length nbins=%d."%nbins
        print "hist.shape: ",hist.shape
        sys.exit()
    
    if len(hist.shape)>2:
        hist = hist.reshape(-1,hist.shape[-1])
    elif len(hist.shape)<2:
        print "Error: expected an n-D array with n>=2 in makebar."
        print hist.shape
        sys.exit()
    
    if val=='median':
        height = np.nanmedian(hist,axis=0)
    elif val=='mean':
        height = np.nanmean(hist,axis=0)
    else:
        print "Error: please enter 'median' or 'mean' for keyword 'val'."
        sys.exit()
        
    if errtype=='full':
        ## error bars give full range of data (max,min)
        ymin = np.nanmin(hist,axis=0)
        ymax = np.nanmax(hist,axis=0)
    elif errtype=='iqr':
        ## error bars give interquartile range of data (0.25-0.75)
        ymin = np.nanpercentile(hist,25,axis=0)
        ymax = np.nanpercentile(hist,75,axis=0)
    elif errtype=='mad':
        ## error bars give median absolute deviation
        mad = np.nanmedian( np.abs( hist - np.nanmedian(hist,axis=0) ), axis=0 ) / np.sqrt(nsim)
        ymin = height - mad
        ymax = height + mad
        if verbose:
            print "in makebar: "
            #print hist.shape
            #print np.abs(hist - np.nanmedian(hist,axis=0)).shape
            #print np.nanmedian(np.abs(hist - np.nanmedian(hist,axis=0))).shape
            #print ymin.shape, ymax.shape
            print "mad: ",mad
            print "ymin: ",ymin
            print "ymax: ",ymax
    else:
        print "Error: please enter 'full' or 'iqr' for keyword 'errtype'."
        sys.exit()
    yerr = (height-ymin, ymax-height)

    if verbose:
        print "bar heights (%s):"%val, height
        print "err min (%s):"%errtype, ymin
        print "err max (%s):"%errtype, ymax

    xvals = xmin + np.arange(nbins) + xoffset
    if xgap_index>=0:
        xvals[xgap_index+1:] = xvals[xgap_index+1:] + 1
    if reverse:
        xvals = -1*xvals[::-1]
        height = height[::-1]
        yerr = (yerr[0][::-1],yerr[1][::-1])

    if binpoints_only:
        return ax.errorbar(xvals, height, yerr=yerr, xerr=width, color=color, 
                           ecolor=ecolor, ls='None',elinewidth=elinewidth,capthick=capthick)
    else:        
        return ax.bar(xvals, height, width=width, color=color, ecolor=ecolor, yerr=yerr,
                      elinewidth=elinewidth,capthick=capthick)


def print_arr_info(arr,nanvals=False):

    if nanvals:
        print "min/max/mean/med: %g %g %g %g"%(np.nanmin(arr),np.nanmax(arr),np.nanmean(arr),np.nanmedian(arr))
    else:
        print "min/max/mean/med: %g %g %g %g"%(arr.min(),arr.max(),arr.mean(),np.median(arr))


def multiplot(maindir='/oasis/projects/nsf/hvd115/lblecha/p3new_merger_sims/',
              subdir_arr=None,tmax_postmrg=0.25,ncam=10,plot_proj_sep=True,
              use_logbins=False, errtype='full'):

    #maindir='/oasis/projects/nsf/hvd115/lblecha/p3new_merger_sims/'

    if use_logbins:
        print "haven't coded histogram arrays to work with use_logbins yet. sorry."
        sys.exit()
        
    if subdir_arr:
        if not isinstance(subdir_arr,list): subdir_arr = [subdir_arr]
        extra = '_'+'_'.join([s for s in subdir_arr])
        print "Using subdir list (%d sims):"%(len(subdir_arr))
        print subdir_arr
    else:
        subdir_arr = ['q1_fg0.3_nomrg/','q1_fg0.1_nomrg/',
                      'q0.5_fg0.3_nomrg/','q0.5_fg0.1_nomrg/',
        #              'q0.5_fg0.1_BT0.2_HRbh_mdyn5e-4_nomrg/',
                      'q0.5_fg0.3_BT0.2_HRbh_mdyn5e-4_nomrg/',
                      'q1_fg0.3_BT0.2_fb0.025_M0.2_largeorb/',
                      'q0.333_fg0.3_allrx10_nomrg','q0.333_fg0.1_allrx10_nomrg',
                      'q0.5_fg0.1_BT0.2_HRbh_mdyn5e-4_allrx10_nomrg',
                      'q0.5_fg0.3_BT0.3_HRbh_mdyn5e-4_allrx10_nomrg']
        newsim_subdir_arr = ['q0.333_fg0.3fg0.1_largeorb',
                             'q0.5_fg0.3fg0.1_BT0.2_largeorb',
                             'q0.5_fg0.3fg0.1_BT0BT0.2_largeorb',
                             'q0.5_fg0.3fg0.1_BT0BT0.2',
                             'q0.5_fg0.3fg0.1_largeorb']
        with_bulges_subdir_arr = ['q0.5_fg0.3_BT0.2_HRbh_mdyn5e-4_nomrg/',
                                  'q1_fg0.3_BT0.2_fb0.025_M0.2_largeorb/',
                                  'q0.5_fg0.1_BT0.2_HRbh_mdyn5e-4_allrx10_nomrg',
                                  'q0.5_fg0.3_BT0.3_HRbh_mdyn5e-4_allrx10_nomrg',
                                  'q0.5_fg0.3fg0.1_BT0.2_largeorb',
                                  'q0.5_fg0.3fg0.1_BT0BT0.2_largeorb']
        fid_orb_subdir_arr = ['q1_fg0.3_nomrg/','q1_fg0.1_nomrg/',
                              'q0.5_fg0.3_nomrg/','q0.5_fg0.1_nomrg/',
                              'q0.5_fg0.3_BT0.2_HRbh_mdyn5e-4_nomrg/',
                              'q0.333_fg0.3_allrx10_nomrg','q0.333_fg0.1_allrx10_nomrg',
                              'q0.5_fg0.1_BT0.2_HRbh_mdyn5e-4_allrx10_nomrg',
                              'q0.5_fg0.3_BT0.3_HRbh_mdyn5e-4_allrx10_nomrg',
                              'q0.5_fg0.3fg0.1_BT0BT0.2',
                              'q0.5_fg0.3fg0.1_largeorb']

        subdir_arr = subdir_arr + newsim_subdir_arr
        extra = ''
        #subdir_arr = fid_orb_subdir_arr
        #extra = '_fid_orb'
        #subdir_arr = with_bulges_subdir_arr
        #extra = '_no_bulgeless'
        minor_mrg_subdir_arr = ['q0.2_fg0.3_BT0.2','q0.2_fg0.3_largeorb',
                                'q0.2_fg0.3_fb0.025','q0.1_fg0.3_BT0.2','q0.1_fg0.3']
        #subdir_arr = minor_mrg_subdir_arr
        #extra = extra+'_minor_only'
        #subdir_arr = subdir_arr+minor_mrg_subdir_arr
        #extra = extra+'_with_minor'

        #test_subdir_arr = ['q1_fg0.3_nomrg/','q1_fg0.1_nomrg/']
        ##test_subdir_arr = ['q1_fg0.3_M0.2/','q0.5_fg0.3fg0.1_BT0.2/']
        #subdir_arr = test_subdir_arr
        #extra = '_test'

        print "Using default subdir list (%d sims):"%len(subdir_arr)
        print subdir_arr

    simpath_arr = ['%s/%s'%(maindir,sub) for sub in subdir_arr]
    nsim = len(simpath_arr)
    
    ### get random sight lines ###
    phi = np.random.rand(ncam) * 2*np.pi
    theta = np.random.rand(ncam) * np.pi

    ### define bins for histograms ###
    #dsep=10
    #sep_bins = np.arange(-1*dsep,80,dsep)
    #dsep=20
    #sep_bins = np.arange(-1*dsep,180,dsep)
    sep_bins = np.array([-3.0, 0, 3, 10, 30, 100, 300])
    #x_sep_bins = sep_bins[:-1]+0.5*(sep_bins[1:]-sep_bins[:-1])
    #print "x_sep_bins: ",x_sep_bins
    dlgsep = 0.5
    #lgsep_bins = np.arange(-1.0,2.5,dlgsep)
    lgsep_bins = np.arange(-4.0,3.0,dlgsep)
    bins = lgsep_bins if use_logbins else sep_bins
    bw=0.3

    ### set limits for agn masks ###
    fedd_lims = [0.01,0.05,0.1]
    lagn_lims = [1.0e43, 1.0e44, 1.0e45]
    hilo_lagn_lims = [10**43.2, 10**44.2]
    #hilo_lagn_lims = [1.0e43, 10**44.2]
    nagnmask = len(fedd_lims)+len(lagn_lims)+(len(hilo_lagn_lims)+1)
    
    ### initialize arrays for combined sims ###
    all_nbh = np.empty((0))
    all_lagn = np.empty((2,0))
    all_mass = np.empty((2,0))
    all_pos = np.empty((2,3,0))
    all_time = np.empty((0))
    all_dt = np.empty((0))
    all_dtfrac = np.empty((0))
    all_ttot = np.empty((0))
    all_tmrg = np.empty((0))
    all_3dsep = np.empty((0))
    all_projsep = np.empty((ncam,0))
    all_totlagn = np.empty((0))
    all_totmass = np.empty((0))
    all_totledd = np.empty((0))
    all_totfedd = np.empty((0))

    ### initialize histogram arrays ###
    hist3d_dt = np.zeros((nsim,len(bins)-1)) + 1.0e-10
    hist3d_dt_frac = np.zeros((nsim,len(bins)-1)) + 1.0e-10
    histproj_dt = np.zeros((nsim,ncam,(len(bins)-1))) + 1.0e-10
    histproj_dt_frac = np.zeros((nsim,ncam,(len(bins)-1))) + 1.0e-10
    agn_hist3d_dt = np.zeros((nsim,nagnmask,len(bins)-1)) + 1.0e-10
    agn_hist3d_dt_frac = np.zeros((nsim,nagnmask,len(bins)-1)) + 1.0e-10
    agn_histproj_dt = np.zeros((nsim,nagnmask,ncam,(len(bins)-1))) + 1.0e-10
    agn_histproj_dt_frac = np.zeros((nsim,nagnmask,ncam,(len(bins)-1))) + 1.0e-10
    
    ### load & process snapshot data ###
    for i_sim,path in enumerate(simpath_arr):

        print "\nloading snapshot data for sim %d (%s)..."%(i_sim,path)
        bhdata = readsnapbhs.all_snapshot_bhs(path,tmax_postmrg=tmax_postmrg)
    
        has2bh = (bhdata.nbh==2)

        dt = np.zeros(bhdata.time.size)
        #dt[1:] = bhdata.time[1:] - bhdata.time[:-1]
        #dt[0] = np.minimum(dt[1],bhdata.time[0])
        dt[:-1] = bhdata.time[1:]-bhdata.time[:-1]
        dt[-1] = 0.0
        ttot = dt.sum()
        dt_tiled = np.tile(dt, (ncam,1))

        bh3dsep = np.repeat(np.nan, bhdata.nbh.size)
        bh3dsep[has2bh] = np.sqrt( np.sum(((bhdata.pos[1,j,has2bh]-bhdata.pos[0,j,has2bh])**2 
                                         for j in range(3)), axis=0) )
        bh3dsep[bh3dsep!=bh3dsep] = -1.0

        rotated_pos1 = rotate_arr(bhdata.pos[0,:,:], theta, phi)
        rotated_pos2 = rotate_arr(bhdata.pos[1,:,:], theta, phi)
        bhprojsep = np.sqrt( np.sum((rotated_pos2[:2,:,:]-
                                     rotated_pos1[:2,:,:])**2,axis=0) )
        bhprojsep[bhprojsep!=bhprojsep] = -1.0
        assert np.array([np.allclose(bh3dsep[has2bh], np.sqrt(np.sum((rotated_pos2[:,j,has2bh]-rotated_pos1[:,j,has2bh])**2,axis=0)) ) for j in range(ncam)]).all(), \
            "Something went wrong in calculating rotated BH positions! Mismatch in 3D sep calculation."
        assert np.array_equal(bhprojsep,bhprojsep), \
            "Something went wrong in calculating rotated BH positions! nan's in bhprojsep."

        totlagn = copy(bhdata.lbol[0,:])
        totlagn[has2bh] = bhdata.lbol[:,has2bh].sum(axis=0)
        
        totmass = copy(bhdata.mass[0,:])
        totmass[has2bh] = bhdata.mass[:,has2bh].sum(axis=0)
        totledd = calc_lbol_edd(totmass)
        totfedd = 0.0*copy(totledd)
        totfedd[totledd>0] = totlagn[totledd>0]/totledd[totledd>0]

        
        ### calc histograms for each simulation ###
        print "calculating histograms for sim %d..."%i_sim
        if len(sep_bins[sep_bins<0])>1 or sep_bins.min()>bhprojsep.min():
            print "Error: invalid definition of sep_bins:",sep_bins
            sys.exit()
        idx_sep_bins = np.digitize(bh3dsep,bins=sep_bins,right=True)
        print "idx_sep_bins:",idx_sep_bins.shape
        print sep_bins.shape, bh3dsep.shape
        idx_projsep_bins = np.array([np.digitize(bhprojsep[j,:],bins=sep_bins,right=True)
                                     for j in range(ncam)])
        print "idx_projsep_bins:",idx_projsep_bins.shape
        print sep_bins.shape, bhprojsep.shape
        ttot_sep_bins = np.histogram(bh3dsep,weights=dt,bins=sep_bins)[0]
        ttot_projsep_bins = np.array([ np.histogram(bhprojsep[j,:],weights=dt,
                                                    bins=sep_bins)[0] for j in range(ncam) ])
        #dt_frac_sep_bins = dt/ttot_sep_bins[idx_sep_bins-1]
        #dt_frac_projsep_bins = dt_tiled/np.array([ttot_projsep_bins[j,idx_projsep_bins[j,:]-1]
        #                                          for j in range(ncam)])    
        print "ttot_sep_bins:",ttot_sep_bins
        print "ttot_sep_bins[idx_sep_bins-1].shape:",ttot_sep_bins[idx_sep_bins-1].shape
        #print "ttot_sep_bins:",ttot_sep_bins[idx_sep_bins-1]
    
        sepvar = np.log10(bh3dsep) if use_logbins else bh3dsep
        projsepvar = np.log10(bhprojsep) if use_logbins else bhprojsep

        hist3d_dt[i_sim,:] = np.histogram(sepvar, bins=bins, weights=dt)[0]
        hist3d_dt_frac[i_sim,:] = hist3d_dt[i_sim,:]/ttot_sep_bins
        #hist3d_dt_frac[i_sim,:] = np.histogram(sepvar, bins=bins,
        #                                       weights=dt_frac_sep_bins)[0]    
        #hist3d_dt[hist3d_dt==0.0] = 1.0e-10
        #hist3d_dt_frac[hist3d_dt_frac==0.0] = 1.0e-10
        print histproj_dt.shape, projsepvar.shape, bins.shape, dt.shape
        histproj_dt[i_sim,:,:] = np.array([np.histogram(projsepvar[j,:],
                                                        bins=bins, weights=dt)[0]
                                           for j in range(ncam)])
        histproj_dt_frac[i_sim,:,:] = histproj_dt[i_sim,:,:]/ttot_projsep_bins
        #histproj_dt_frac[i_sim,:,:] = np.array([np.histogram(projsepvar[j,:],bins=bins,
        #                                                     weights=dt_frac_projsep_bins[j,:])[0]
        #                                        for j in range(ncam)])
        #histproj_dt[histproj_dt==0.0] = 1.0e-10
        #histproj_dt_frac[histproj_dt_frac==0.0] = 1.0e-10

        fedd_maskarr = [(totfedd>=lim) for lim in fedd_lims]
        lagn_maskarr = [(totlagn>=lim) for lim in lagn_lims]
        hilo_lagn_maskarr = [(totlagn<hilo_lagn_lims[0]),
                             ((totlagn>=hilo_lagn_lims[0])&(totlagn<hilo_lagn_lims[1])),
                             (totlagn>=hilo_lagn_lims[1])]
        agn_maskarr_list = fedd_maskarr + lagn_maskarr + hilo_lagn_maskarr

        for i_mask,mask in enumerate(agn_maskarr_list):
            agn_hist3d_dt[i_sim,i_mask,:] = np.histogram(sepvar[mask],bins=bins,weights=dt[mask])[0]
            agn_hist3d_dt_frac[i_sim,i_mask,:] = agn_hist3d_dt[i_sim,i_mask,:]/ttot_sep_bins
            #agn_hist3d_dt_frac[i_sim,i_mask,:] = np.histogram(sepvar[mask],bins=bins,
            #                                                  weights=dt_frac_sep_bins[mask])[0]
            agn_histproj_dt[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,mask],bins=bins,
                                                                       weights=dt[mask])[0] for j in range(ncam)])
            agn_histproj_dt_frac[i_sim,i_mask,:,:] = agn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            #agn_histproj_dt_frac[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,mask], bins=bins, weights=dt_frac_projsep_bins[j,mask])[0] for j in range(ncam)])
            
        
        ### append sim data to combined arrays ###
        print "appending data for sim %d to combined arrays..."%i_sim
        all_nbh = np.append(all_nbh, bhdata.nbh)
        all_lagn = np.append(all_lagn, bhdata.lbol, axis=1)
        all_mass = np.append(all_mass, bhdata.mass, axis=1)
        all_pos = np.append(all_pos, bhdata.pos, axis=2)
        all_time = np.append(all_time, bhdata.time)
        all_tmrg = np.append(all_tmrg, bhdata.tmrg)
        print "max time=%g, ttot=%g"%(bhdata.time.max(), ttot)
        print "tmrg=%g, postmrg time=%g"%(bhdata.tmrg, bhdata.time.max()-bhdata.tmrg)
        
        all_dt = np.append(all_dt, dt)
        all_dtfrac = np.append(all_dtfrac, dt/ttot)
        all_ttot = np.append(all_ttot, ttot)

        all_3dsep = np.append(all_3dsep, bh3dsep)
        all_projsep = np.append(all_projsep, bhprojsep, axis=1)  
        all_totlagn = np.append(all_totlagn, totlagn)
        all_totmass = np.append(all_totmass, totmass)
        all_totledd = np.append(all_totledd, totledd)
        all_totfedd = np.append(all_totfedd, totfedd)


    print "bins:",bins
    print "hist3d_dt:",hist3d_dt.shape
    print "hist3d_dt_frac:",hist3d_dt_frac.shape
    print "histproj_dt:",histproj_dt.shape
    print "histproj_dt_frac:",histproj_dt_frac.shape
    #print "histproj_dt_frac:",histproj_dt_frac

    #print "\nbins:",bins
    print "agn_hist3d_dt:",agn_hist3d_dt.shape
    print "agn_hist3d_dt_frac:",agn_hist3d_dt_frac.shape
    print "agn_histproj_dt:",agn_histproj_dt.shape
    print "agn_histproj_dt_frac:",agn_histproj_dt_frac.shape

    #print "agn_hist3d_dt_frac and agn_histproj_dt_frac hilo lum bin sums:"
    #for k in range(len(bins)-1):
    #    print bins[k]
    #    print agn_hist3d_dt_frac[:,6:,k].sum(axis=1)
    #    print np.array([agn_histproj_dt_frac[:,6:,j,k].sum(axis=1)
    #                    for j in range(ncam)]).transpose()
        
        
    if plot_proj_sep:

        tmp_dttot_ratio_proj = histproj_dt[:,:,2]/histproj_dt[:,:,1]

        print "\n\ntotal dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_dttot_ratio_proj,nanvals=True)
        print np.nanmedian(tmp_dttot_ratio_proj,axis=1)
        print np.nanmin(tmp_dttot_ratio_proj,axis=1)
        print np.nanmax(tmp_dttot_ratio_proj,axis=1)
        print "\ntotal dt[3-10kpc]:"
        print_arr_info(histproj_dt[:,:,2],nanvals=True)
        print np.nanmedian(histproj_dt[:,:,2],axis=1)
        print "\ntotal dt[0-3kpc]:"
        print_arr_info(histproj_dt[:,:,1],nanvals=True)
        print np.nanmedian(histproj_dt[:,:,1],axis=1)

    else:

        tmp_dttot_ratio_3d = hist3d_dt[:,2]/hist3d_dt[:,1]

        print "\n\ntotal dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_dttot_ratio_3d,nanvals=True)
        print tmp_dttot_ratio_3d
        print "total dt[0-3kpc]:"
        print_arr_info(hist3d_dt[:,1],nanvals=True)
        print hist3d_dt[:,1]


    return

    has2bh = (all_nbh==2)
    #all_totlagn = copy(all_lagn[0,:])
    #all_totlagn[has2bh] = all_lagn[:,has2bh].sum(axis=0)

    #all_totmass = copy(all_mass[0,:])
    #all_totmass[has2bh] = all_mass[:,has2bh].sum(axis=0)
    #all_totledd = calc_lbol_edd(all_totmass)
    #all_totfedd = 0.0*copy(all_totledd)
    #all_totfedd[all_totledd>0] = all_totlagn[all_totledd>0]/all_totledd[all_totledd>0]

    all_lagn_ratio = np.zeros(all_nbh.size)
    all_lagn_ratio[has2bh] = all_lagn[0,has2bh]/all_lagn[1,has2bh]
    lg_lagn_ratio = np.log10(all_lagn_ratio)

    all_dt_tiled = np.tile(all_dt, (ncam,1))  
  
    print all_3dsep.shape, all_projsep.shape, all_dtfrac.shape, all_dt_tiled.shape,
    print all_nbh.shape, all_lagn.shape, all_totfedd.shape

    fedd_tiled = np.tile(all_totfedd,(ncam,1)) if plot_proj_sep else all_totfedd
    lagn_tiled = np.tile(all_totlagn,(ncam,1)) if plot_proj_sep else all_totlagn
    #fedd_maskarr = [(fedd_tiled>=lim) for lim in [0.01,0.05,0.1]]
    fedd_titles = ['fEdd>%g'%(lim) for lim in [0.01,0.05,0.1]]
    #lagn_maskarr = [(lagn_tiled>=lim) for lim in [1.0e43, 1.0e44, 1.0e45]]
    lagn_titles = ['Lbol>%g'%(lim)  for lim in [1.0e43, 1.0e44, 1.0e45]]

    #agn_maskarr_list = fedd_maskarr + lagn_maskarr
    titles = fedd_titles + lagn_titles
    
    sepstr = '_bhprojsep' if plot_proj_sep else '_bhsep'
    lgsepstr = '_lgbhprojsep' if plot_proj_sep else '_lgbhsep'
    sepstring = lgsepstr if use_logbins else sepstr


    #plt.close()
    #plt.clf()
    #plt.cla()
    plt.ioff()
    matplotlib.rcParams.update({'font.size':10})                                     

    #tothist_kwargs = dict(range=(sep_bins.min(),sep_bins.max()),bins=len(sep_bins)-1,
    #                      histtype='bar')
    tothist_kwargs = dict(bins=sep_bins, histtype='bar')
    #histtype='step')                                                            
    totlghist_kwargs = dict(range=(lgsep_bins.min(),lgsep_bins.max()),
                            bins=len(lgsep_bins)-1,  histtype='bar')

    #if plot_proj_sep:
    #    hist,tmp = np.histogram(sepvar,bins=sep_bins,weights=all_dt_tiled/(1.0*len(subdir_arr)*ncam))
    #    agn_hist,tmp = np.histogram(sepvar[agn_maskarr_list[i]],bins=sep_bins,
    #                                weights=all_dt_tiled[i]/(1.0*len(subdir_arr)*ncam))

    if plot_proj_sep:
        sepvar = np.log10(all_projsep) if use_logbins else all_projsep
        wvar = all_dt_tiled/(1.0*len(subdir_arr)*ncam)
        xlabel = 'proj. sep [%skpc]'%('log ' if use_logbins else '')
    else:
        sepvar = np.log10(all_3dsep) if use_logbins else all_3dsep
        wvar = all_dt/(1.0*len(subdir_arr))
        xlabel = '3D sep [%skpc]'%('log ' if use_logbins else '')
    idx_all_sep_bins = np.digitize(all_3dsep,bins=sep_bins,right=True)
    print idx_all_sep_bins.shape, all_projsep.shape, sep_bins.shape
    idx_all_projsep_bins = np.array([np.digitize(all_projsep[i,:],bins=sep_bins,right=True)
                                 for i in range(ncam)])
    bins = lgsep_bins if use_logbins else sep_bins
    bw=0.3

    #hist,tmp = np.histogram(sepvar,bins=bins,weights=wvar)
    #print "tmp: ",tmp
    #print "hist: ",hist
    #print "bins: ",bins
    #if use_logbins:
    #    xticklbl = ["%g-%g"%(bins[j],bins[j+1]) for j in np.arange(len(hist))]
    #else:
    #    xticklbl = ['0']+["%g-%g"%(bins[j],bins[j+1]) for j in np.arange(1,len(hist))]
    xticklbl = ["%g-%g"%(bins[j],bins[j+1]) for j in np.arange(len(bins)-1)]
    if not use_logbins: xticklbl[0] = '0'

    #if nagnmask !=6:
    #    print "WARNING!! plotting 6 values for agn mask; actual number is %d."%nagnmask
        
    fig = plt.figure(figsize=(9,6))
    for i in range(6):
        #agn_hist,tmp = np.histogram(sepvar[agn_maskarr_list[i]],bins=bins,
        #                            weights=wvar[agn_maskarr_list[i]])
        ax = fig.add_subplot(231+i)
        #ax.set_yscale('log')
        #ax.set_ylim(1.0e-3,2)
        ax.set_ylim(0.0,1.41)
        ax.set_ylabel(r'$\Delta$t [Gyr]') if  i in [0,3] else ax.set_ylabel('')
        ax.set_xticks(np.arange(len(bins)-1))
        ax.set_xticklabels(xticklbl)

        if plot_proj_sep:
            tot=makebar(ax,histproj_dt,nbins=len(bins)-1,color='c',width=bw,
                        xoffset=-bw,val='median',errtype=errtype)
            agn=makebar(ax,agn_histproj_dt[:,i,:,:],nbins=len(bins)-1,color='m',width=bw,
                        xoffset=0,val='median',errtype=errtype)
        else:
            tot=makebar(ax,hist3d_dt,nbins=len(bins)-1,color='c',width=bw,xoffset=-bw,
                        val='median',errtype=errtype)
            agn=makebar(ax,agn_hist3d_dt[:,i,:],nbins=len(bins)-1,color='m',width=bw,
                        xoffset=0,val='median',errtype=errtype)

        #print hist3d_dt.shape, hist3d_dt.mean(axis=0).shape
        #print histproj_dt.shape, histproj_dt.mean(axis=0).shape
        #print agn_histproj_dt.shape, agn_histproj_dt[:,i,:,:].mean(axis=0).shape
        #ax.hist( [sepvar[agn_maskarr_list[i]].flatten(), sepvar.flatten()],
        #         weights=[all_dt_tiled[agn_maskarr_list[i]].flatten()/(1.0*len(subdir_arr)*ncam),
        #                  all_dt_tiled.flatten()/(1.0*len(subdir_arr)*ncam)],
        #         color=['m','k'],lw=2,stacked=False,
        #         label=['AGN','Total'],**tothist_kwargs)
        ax.set_xlabel(xlabel) if i in [3,4,5] else ax.set_xlabel('')
        ax.legend((tot,agn),('Total','AGN'),fontsize=9,loc='upper left')
        ax.set_title(titles[i],fontsize=10)
    #fig.suptitle(" ".join(subdir_arr),fontsize=8)
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,right=0.94,bottom=0.1,top=0.92)
    if nsim==1:
        plotname='%s/%s/merger_agn_phases_multiplot%s.pdf'%(maindir,subdir_arr[0],sepstring)
    else:
        plotname='%s/merger_agn_phases_multiplot%s%s.pdf'%(maindir,sepstring,extra)
    fig.savefig(plotname,format='pdf')
    plt.cla()
    plt.clf()

    fig = plt.figure(figsize=(9,6))
    for i in range(6):
        ax = fig.add_subplot(231+i)
        #ax.set_yscale('log')
        #ax.set_ylim(1.0e-3,1.5)
        ax.set_ylim(0.0,1.05)
        ax.set_ylabel(r'$\Delta$t/t$_{\rm bin}$',fontsize=11) if  i in [0,3] else ax.set_ylabel('')
        ax.set_xticks(np.arange(len(bins)-1))
        ax.set_xticklabels(xticklbl)

        if plot_proj_sep:
            #tot=makebar(ax,histproj_dt_frac,nbins=len(bins)-1,color='c',width=bw,
            #            xoffset=-1.5*bw,val='median',errtype=errtype)
            agn=makebar(ax,agn_histproj_dt_frac[:,i,:,:],nbins=len(bins)-1,color='m',width=2*bw,
                        xoffset=-bw,val='median',errtype=errtype)
        else:
            #tot=makebar(ax,hist3d_dt_frac,nbins=len(bins)-1,color='c',width=bw,
            #            xoffset=-bw,val='median',errtype=errtype)
            agn=makebar(ax,agn_hist3d_dt_frac[:,i,:],nbins=len(bins)-1,color='m',width=2*bw,
                        xoffset=-bw,val='median',errtype=errtype)

        ax.set_xlabel(xlabel) if i in [3,4,5] else ax.set_xlabel('')
        #ax.legend((tot,agn),('Total','AGN'),fontsize=9,loc='upper right')
        ax.set_title(titles[i],fontsize=10)
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,right=0.94,bottom=0.1,top=0.92)
    if nsim==1:
        plotname='%s/%s/merger_agn_phases_dtfrac_multiplot%s.pdf'%(maindir,subdir_arr[0],sepstring)
    else:
        plotname='%s/merger_agn_phases_dtfrac_multiplot%s%s.pdf'%(maindir,sepstring,extra)
    fig.savefig(plotname,format='pdf')
    plt.cla()
    plt.clf()

    ### separate high- and low-lum AGN phases:

    bw = 0.2
    fig = plt.figure(figsize=(6,6))

    ax1 = fig.add_subplot(211)
    ax1.set_ylim(0.0,1.35)
    ax1.set_ylabel(r'$\Delta$t [Gyr]')
    ax1.set_xticks(np.arange(len(bins)-1))
    ax1.set_xticklabels(xticklbl)

    if plot_proj_sep:
        tot=makebar(ax1,histproj_dt,nbins=len(bins)-1,color='c',width=bw,
                    xoffset=-2*bw,val='median',errtype=errtype)
        hiagn=makebar(ax1,agn_histproj_dt[:,8,:,:],nbins=len(bins)-1,color='b',width=bw,
                      xoffset=-bw,val='median',errtype=errtype)
        loagn=makebar(ax1,agn_histproj_dt[:,7,:,:],nbins=len(bins)-1,color='r',width=bw,
                      xoffset=0,val='median',errtype=errtype)
        noagn=makebar(ax1,agn_histproj_dt[:,6,:,:],nbins=len(bins)-1,color='g',width=bw,
                      xoffset=bw,val='median',errtype=errtype)
    else:
        tot=makebar(ax1,hist3d_dt,nbins=len(bins)-1,color='c',width=bw,xoffset=-2*bw,
                    val='median',errtype=errtype)
        hiagn=makebar(ax1,agn_hist3d_dt[:,8,:],nbins=len(bins)-1,color='b',width=bw,
                      xoffset=-bw,val='median',errtype=errtype)
        loagn=makebar(ax1,agn_hist3d_dt[:,7,:],nbins=len(bins)-1,color='r',width=bw,
                      xoffset=0,val='median',errtype=errtype)
        noagn=makebar(ax1,agn_hist3d_dt[:,6,:],nbins=len(bins)-1,color='g',width=bw,
                      xoffset=bw,val='median',errtype=errtype)
    
    ax1.set_xlabel(xlabel)
    ax1.legend((tot,hiagn,loagn,noagn),('Total','Luminous AGN','Low-lum AGN','Inactive'),
               fontsize=9,loc='upper left')

    ax2 = fig.add_subplot(212)
    ax2.set_ylim(0.0,1.1)
    ax2.set_ylabel(r'$\Delta$t/t$_{\rm bin}$',fontsize=11)
    ax2.set_xticks(np.arange(len(bins)-1))
    ax2.set_xticklabels(xticklbl)

    if plot_proj_sep:
        print "agn_histproj_dt_frac(log(lagn)>=%g):"%hilo_lagn_lims[1]
        hiagn=makebar(ax2,agn_histproj_dt_frac[:,8,:,:],nbins=len(bins)-1,color='b',width=bw,
                      xoffset=-bw,val='median',errtype=errtype,verbose=True)
        print "agn_histproj_dt_frac(%g<=log(lagn)<%g):"%(hilo_lagn_lims[0],hilo_lagn_lims[1])
        loagn=makebar(ax2,agn_histproj_dt_frac[:,7,:,:],nbins=len(bins)-1,color='r',width=bw,
                      xoffset=0,val='median',errtype=errtype,verbose=True)
        print "agn_histproj_dt_frac(log(lagn)<%g):"%hilo_lagn_lims[0]
        noagn=makebar(ax2,agn_histproj_dt_frac[:,6,:,:],nbins=len(bins)-1,color='g',width=bw,
                      xoffset=bw,val='median',errtype=errtype,verbose=True)
    else:
        print "agn_hist3d_dt_frac(log(lagn)>=%g):"%hilo_lagn_lims[1]
        hiagn=makebar(ax2,agn_hist3d_dt_frac[:,8,:],nbins=len(bins)-1,color='b',width=bw,
                      xoffset=-bw,val='median',errtype=errtype,verbose=True)
        print "agn_hist3d_dt_frac(%g<=log(lagn)<%g):"%(hilo_lagn_lims[0],hilo_lagn_lims[1])
        loagn=makebar(ax2,agn_hist3d_dt_frac[:,7,:],nbins=len(bins)-1,color='r',width=bw,
                      xoffset=0,val='median',errtype=errtype,verbose=True)
        print "agn_hist3d_dt_frac(log(lagn)<%g):"%hilo_lagn_lims[0]
        noagn=makebar(ax2,agn_hist3d_dt_frac[:,6,:],nbins=len(bins)-1,color='g',width=bw,
                      xoffset=bw,val='median',errtype=errtype,verbose=True)
    
    ax2.set_xlabel(xlabel)
    ax2.legend((hiagn,loagn,noagn),('Luminous AGN','Low-lum AGN','Inactive'),
               fontsize=9,loc='upper right')

    #fig.suptitle(" ".join(subdir_arr),fontsize=8)
    fig.subplots_adjust(wspace=0.2,hspace=0.25,left=0.1,right=0.96,bottom=0.1,top=0.94)
    if nsim==1:
        plotname='%s/%s/merger_lum_agn_phases%s.pdf'%(maindir,subdir_arr[0],sepstring)
    else:
        plotname='%s/merger_lum_agn_phases%s%s.pdf'%(maindir,sepstring,extra)
    fig.savefig(plotname,format='pdf')
    plt.cla()
    plt.clf()
