import sed, pylab
from pylab import *
import pyfits
import astro_constants as ac
from copy import copy

#snap = 470
#maindir = '/oasis/projects/nsf/hvd115/lblecha/q0.5_fg0.3_allrx10_sunruns/whole_sed_lowres_newcode/'
#path1=maindir
#path2 = maindir+'/test_irequilib0.2'
#compare_fbase='irequlib0.2'
##path2 = maindir+'/test_irequilib0.4'
##compare_fbase='irequlib0.4'
#path3 = maindir+'/test_irequilib0.4'
#compare_fbase='irequlib'

#cam=3


def compare(maindir='/oasis/projects/nsf/hvd115/lblecha/',rundir='q0.5_fg0.3_sunruns/',
            main_subdir='all',compare_subdirs=None,compare_fbase='',snap_arr=[5]):

    if not compare_subdirs: 
        "You must enter at least one subdir for keyword 'compare_subdirs'" 
        return

    maindir = maindir+'/'+rundir 
    if not isinstance(compare_subdirs,list): 
        compare_subdirs = [compare_subdirs]


    #path0=maindir+'/test_irequlib0.05'
    #path1=maindir+'/test_fiducial_hires'
    #path2=maindir+'/test_irequlib0.2'
    #path3=maindir+'/test_irequlib0.4'
    #compare_fbase='irequilib_no0.05'
        
    #for snap in (100,150,250):
    #for snap in (100,150,200,212,250):
    for snap in snap_arr:
        print "making sed plots for snap %d."%snap
        for cam in range(7):
            print "comparison seds for cam %d"%cam
            #sed.plot_seds([(p+'/mcrx_%.3d.fits'%snap,'L_lambda_nonscatter%d'%cam) 
            #               for p in (path0,path1,path2,path3)],legend=False,relative=True)
            sed.plot_seds([(maindir+'/'+p+'/mcrx_%.3d.fits'%snap,'L_lambda_nonscatter%d'%cam) 
                           for p in [main_subdir]+compare_subdirs],legend=False,relative=True)
            pylab.savefig(maindir+'/sed_compare_%s_nonscatter_cam%d_%.3d.pdf'%(compare_fbase,cam,snap))
            sed.plot_seds([(maindir+'/'+p+'/mcrx_%.3d.fits'%snap,'L_lambda_scatter%d'%cam) 
                           for p in [main_subdir]+compare_subdirs],legend=False,relative=True)
            pylab.savefig(maindir+'/sed_compare_%s_scatter_cam%d_%.3d.pdf'%(compare_fbase,cam,snap))
            

def plot_allcam(maindir='/oasis/projects/nsf/hvd115/lblecha/',rundir='q0.5_fg0.3_sunruns/',
                subdir='all',snap_arr=[5],scatter=True):

    path = maindir+'/'+rundir+'/'+subdir
    scatstr='scatter' if scatter else 'nonscatter'

    for snap in snap_arr:
        #print "plotting %s relative sed for snap %d."%(scatstr,snap)
        #sed.plot_seds([(path+'/mcrx_%.3d.fits'%snap,'L_lambda_%s%d'%(scatstr,i)) 
        #               for i in range(7)],legend=False,relative=True)
        #pylab.savefig(path+'/sed_relative_%s_%.3d.pdf'%(scatstr,snap))

        #print "plotting %s sed for snap %d."%(scatstr,snap)
        #sed.plot_seds([(path+'/mcrx_%.3d.fits'%snap,'L_lambda_%s%d'%(scatstr,i)) 
        #               for i in range(7)],legend=False,relative=False)
        #pylab.savefig(path+'/sed_%s_%.3d.pdf'%(scatstr,snap))

        print "plotting test %s sed for snap %d."%(scatstr,snap)
        g=sed.pyxplot_seds([(path+'/mcrx_%.3d.fits'%snap,'L_lambda_%s%d'%(scatstr,i)) 
                            for i in range(7)],range=[3.95e-7,4.0e-7,1e36,5e37],legend=True,
                           relative=False,ytitle='',legendtxt=['c%d'%i for i in range(7)])
        g.writePDFfile(path+'/test_sed_%s_%.3d.pdf'%(scatstr,snap))
        #pylab.savefig(path+'/test_sed_%s_%.3d.pdf'%(scatstr,snap))


 
                                                                                                                         

def plot_allcam_manual(maindir='/oasis/projects/nsf/hvd115/lblecha/',rundir='q0.5_fg0.3_sunruns/',
                       subdir='all',snap_arr=[5],scatter=True,axrange=[3e-7,7e-7,1e36,5e38],extra=''):

    path = maindir+'/'+rundir+'/'+subdir
    scatstr='scatter' if scatter else 'nonscatter'
    if extra!='': extra=extra+'_'

    for snap in snap_arr:

        seds=sed.extract_seddata([(path+'/mcrx_%.3d.fits'%snap,'L_lambda_%s%d'%(scatstr,i)) 
                                  for i in range(7)],False)
 
        f=plt.figure(figsize=(5,3))

        plt.yscale('log')
        plt.xlim([axrange[0],axrange[1]])
        plt.ylim([axrange[2],axrange[3]])

        for l,s,t in seds:
            print l.shape, s.shape, t
            plt.plot(l,s)
            
        f.subplots_adjust(bottom=0.18)
        f.savefig(path+'/test_sed_%s_%s%.3d.pdf'%(scatstr,extra,snap))
        plt.clf()
        plt.close()


def plot_stellarmodel_seds(path='/home/lblecha/sundata/chris_stellarmodel/',
                           fname='logspace-Patrik-imfKroupa-geneva-Zmulti-hires.fits',ylog=True,extra='',
                           nvals=4,xlim=(1e-8,3e-6),ylim=None,relative=False,lambda_ref=1e-6,lambda_units='m'):
   
    ### test_Smodel_hires_testlrange2: 
    ###     stellar model: ~/sundata/chris_stellarmodel/Patrik-imfKroupa-Zmulti-ml.fits
    ###     mappings model ~/sundata/chris_mappings/Smodel_full_hires.fits
    ### hirestest_mappings_kinematics_testlrange: 
    ###     stellar model: ~/sundata/chris_stellarmodel/logspace-Patrik-imfKroupa-geneva-Zmulti-hires.fits
    ###     mappings model: ~/sundata/chris_mappings/hirestest.fits

    assert lambda_units in ('m','A')
    relstr='_relative' if relative else ''
    if extra!='': extra='_'+extra

    f=pyfits.open('%s/%s'%(path,fname))
    print f['SED'].data.shape
    print f['current_mass'].data.shape
    print f['SED'].data.shape
    print f['AXES'].data.shape
    print f['AXES'].data.field('lambda').shape
    time=f['AXES'].data.field('time')
    Z=f['AXES'].data.field('metallicity')
    l=f['AXES'].data.field('lambda')
    print time.shape, Z.shape, l.shape
    time = time[time>0]
    Z = Z[Z>0]
    l = l[l>0]
    if lambda_units=='A': 
        l=l*1e10
        lambda_ref=lambda_ref*1e10
        xlim=[x*1e10 for x in xlim]
    if relative: 
        lref_index=np.where(np.abs(l-lambda_ref)==np.min(np.abs(l-lambda_ref)))[0][0]
        print lref_index
        assert isinstance(lref_index,int)
    print time.shape, Z.shape, l.shape
    seds = f['SED'].data[:,::(Z.size/nvals),::(time.size/nvals)]
    f.close()
    print seds.shape, seds.shape[1],seds.shape[2]



    if not ylim: ylim=(1e30,1e38) if not relative else (1e-10,1e4)
    fig = plt.figure(figsize=(5,3))
    ax = fig.add_subplot(111)
    if xlim[1]/(1.0*xlim[0])>100: ax.set_xscale('log')
    ax.set_xlim(xlim)
    if ylog: ax.set_yscale('log')
    ax.set_ylim(ylim)
    ax.set_xlabel(r'$\lambda$ [%s]'%lambda_units)
    ax.set_ylabel('Flux (relative)') if relative else ax.set_ylabel('Flux [W/m]')
    #plt.plot((3728.30e-10,3728.30e-10),ylim)
    #plt.xlim(3.7e-7,4e-7)
    #plt.ylim(np.log10(ylim))
    #plt.ylim(np.log10(ylim))

    for i in range(seds.shape[1]):
        for j in range(seds.shape[2]):
            s=10**seds[:,i,j]
            if relative: s=s/10**seds[lref_index,i,j]
            ax.plot(l,s)

    #print ax.get_xticks()
    #print ax.get_xlim()
    #print lambda_ref
    #if len(ax.get_xticks())<3: ax.set_xticks(np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/8.))
    #print ax.get_xticks()
    fig.subplots_adjust(bottom=0.18,left=0.16)
    fig.savefig('%s/select_seds%s%s_%s.pdf'%(path,relstr,extra,fname.replace('.fits','')))
    plt.clf()
    plt.close()


def plot_mappings_seds(path='/home/lblecha/sundata/chris_mappings/',
                       fname='hirestest.fits',ylog=True,extra='',nvals=2,xlim=(1e-8,3e-6),
                       ylim=None,relative=False,lambda_ref=1e-6,lambda_units='m'):
   
    ### test_Smodel_hires_testlrange2: 
    ###     mappings model ~/sundata/chris_mappings/Smodel_full_hires.fits
    ### hirestest_mappings_kinematics_testlrange: 
    ###     mappings model: ~/sundata/chris_mappings/hirestest.fits

    assert lambda_units in ('m','A')
    relstr='_relative' if relative else ''
    if extra!='': extra='_'+extra


    f=pyfits.open('%s/%s'%(path,fname))
    print f[1].data.shape
    print f[2].data.shape
    print f[1].data.field('WAVE').shape
    l=f[1].data.field('WAVE').astype('float64')
    ## convert from microns to desired units:
    l = l[:-1]*1.0e4 if lambda_units=='A' else l[:-1]*1.0e-6
    if lambda_units=='A': 
        lambda_ref=lambda_ref*1e10
        xlim=[x*1e10 for x in xlim]
    if relative: 
        lref_index=np.where(np.abs(l-lambda_ref)==np.min(np.abs(l-lambda_ref)))[0][0]
        print lref_index
        assert isinstance(lref_index,int)

    # indices are: WAVE, COVER, PRESSURE, SPARAM, METAL
    dshape = f[2].data.shape[1:]
    print dshape
    #seds = f[2].data[:-1,::np.maximum(dshape[0]/nvals,1),
    #                   ::np.maximum(dshape[1]/nvals,1),
    #                   ::np.maximum(dshape[2]/nvals,1),
    #                   ::np.maximum(dshape[3]/nvals,1)]
    seds = f[2].data[:-1,:]
    #seds = f[2].data[:-1,:,::2,1::2,::2]
    print l.min(),l.max()
    print seds.min(),seds.max()
    print seds.shape, seds.shape[1],seds.shape[2]
    #tmp=f[2].data[:-1,-1,-1,-1,-1]
    #print tmp.min(),tmp.max(),tmp.shape
    f.close()
    seds = seds.reshape(seds.shape[0],seds.shape[1]*seds.shape[2]*seds.shape[3]*seds.shape[4])
    print seds.shape

    if not ylim: ylim=(1e16,1e33) if not relative else (1e-10,1e4)
    fig = plt.figure(figsize=(5,3))
    ax = fig.add_subplot(111)
    if xlim[1]/(1.0*xlim[0])>100: ax.set_xscale('log')
    ax.set_xlim(xlim)
    if ylog: ax.set_yscale('log')
    ax.set_ylim(ylim)
    ax.set_xlabel(r'$\lambda$ [%s]'%lambda_units)
    ax.set_ylabel('Flux (relative)') if relative else ax.set_ylabel('Flux erg/s/Hz')
    #plt.plot((3728.30e-10,3728.30e-10),ylim)
    #plt.xlim(3.7e-7,4e-7)
    #plt.ylim(np.log10(ylim))
    #plt.ylim(np.log10(ylim))

    for i in range(seds.shape[1]):
            s=seds[:,i]
            if relative: s=s/seds[lref_index,i]
            ax.plot(l,s)
    #if relative: tmp=tmp/tmp[lref_index]
    #ax.plot(l,tmp,color='darkorange')

    #print ax.get_xticks()
    #print ax.get_xlim()
    #print lambda_ref
    #if len(ax.get_xticks())<3: ax.set_xticks(np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/8.))
    #print ax.get_xticks()
    fig.subplots_adjust(bottom=0.18,left=0.16)
    fig.savefig('%s/select_seds%s%s_%s.pdf'%(path,relstr,extra,fname.replace('.fits','')))
    plt.clf()
    plt.close()


def compare_manual(maindir='/oasis/projects/nsf/hvd115/lblecha/',rundir='q0.5_fg0.3_sunruns/',
                   main_subdir='all',compare_subdirs=None,compare_fbase='',snap_arr=[5],
                   ylog=True,extra='',nvals=2,xlim=(8.0e-8,1.1e-3),ylim=None,sed_type='scatter',
                   lambda_ref=1e-6,lambda_units='m',flux_units='nuFnu',ncam=7,oplot_wise=False,
                   single_cam=-1,snap_titles=True,snap_labels=[],single_snap_plots=True):

    wcens = [3.368e-6, 4.618e-6, 12.082e-6, 22.194e-6]
    wweff = [0.6626e-6, 1.0423e-6, 5.5056e-6, 4.1017e-6]

    if not compare_subdirs: 
        compare_subdirs=[main_subdir]

    if len(compare_fbase)>0: compare_fbase='_'+compare_fbase

    if not isinstance(compare_subdirs,list): 
        compare_subdirs = [compare_subdirs]

    if snap_labels: 
        assert len(snap_labels)==len(snap_arr), 'Error: snap_labels must have same length as snap_arr.'

    if extra: extra = '_'+extra

    assert lambda_units in ('m','A','microns')
    assert flux_units in ('relative','Flambda','nuFnu')
    assert sed_type in ('scatter','nonscatter','ir','all','tot')
    relstr='_relative' if flux_units=='relative' else ''

    if lambda_units=='A': 
        xlim=[x*1e10 for x in xlim]
        wcens=[w*1e10 for w in wcens]
        wweff=[w*1e10 for w in wweff]
    elif lambda_units=='microns':
        xlim=[x*1e6 for x in xlim]
        wcens=[w*1e6 for w in wcens]
        wweff=[w*1e6 for w in wweff]

    #path0=maindir+'/test_irequlib0.05'
    #path1=maindir+'/test_fiducial_hires'
    #path2=maindir+'/test_irequlib0.2'
    #path3=maindir+'/test_irequlib0.4'
    #compare_fbase='irequilib_no0.05'
        
    #for snap in (100,150,250):
    #for snap in (100,150,200,212,250):


    #if not ylim: ylim=(1e30,1e38) if not relative else (1e-10,1e4)
    if flux_units=='relative': 
        if not ylim: ylim=(0.9,1.1)
        ylbl=('L$_{\lambda}$ [relative]')
    elif flux_units=='Flambda':
        if not ylim: ylim=(1e36,1e45)
        ylbl=(r'L$_{\lambda}$ [W m$^{-1}$]')
    elif flux_units=='nuFnu':
        if not ylim: ylim=(1e31,1e40)
        ylbl=(r'$\nu$ L$_{\nu}$ [W]')
        
    l=np.empty(0)
    #seds=np.empty(0)
    ns_seds=np.empty(0)
    scat_seds=np.empty(0)
    ir_seds=np.empty(0)
    for i_snap,snap in enumerate(snap_arr):
        print "loading sed data for snap %d."%snap

        #l=np.empty(0)
        #seds=np.empty(0)
        for i_dir,subdir in enumerate(compare_subdirs):
            path = maindir+'/'+rundir+'/'+subdir

            f=pyfits.open('%s/mcrx_%.3d.fits'%(path,snap))
            iq=f['INTEGRATED_QUANTITIES']
            if l.size==0: 
                l=iq.data.field('lambda')
                if lambda_units=='A': 
                    l_m = copy(l)
                    l=l*1e10
                    lu_lbl = r'$\AA$'
                elif lambda_units=='microns':
                    l_m = copy(l)
                    l=l*1.0e6
                    lu_lbl = r'$\mu$m'
                else: lu_lbl = lambda_units
            if ns_seds.size==0: ns_seds=np.zeros((l.size,ncam,len(compare_subdirs),len(snap_arr)))
            if scat_seds.size==0: scat_seds=np.zeros((l.size,ncam,len(compare_subdirs),len(snap_arr)))
            if ir_seds.size==0: ir_seds=np.zeros((l.size,ncam,len(compare_subdirs),len(snap_arr)))
            #seds[:,:,i_dir,i_snap]=np.array([iq.data.field('L_lambda_%s%d'%(sed_type,c)) 
            #                                 for c in range(ncam)]).transpose()
            ns_seds[:,:,i_dir,i_snap]=np.array([iq.data.field('L_lambda_%s%d'%('nonscatter',c)) 
                                                for c in range(ncam)]).transpose()
            scat_seds[:,:,i_dir,i_snap]=np.array([iq.data.field('L_lambda_%s%d'%('scatter',c)) 
                                                  for c in range(ncam)]).transpose()
            ir_seds[:,:,i_dir,i_snap]=np.array([iq.data.field('L_lambda_%s%d'%('ir',c)) 
                                                for c in range(ncam)]).transpose()
            #print 'seds.shape:'
            #print seds.shape
            f.close()
            for s in (ns_seds,scat_seds,ir_seds):
                assert s.shape[1]==ncam
                assert s.shape[2]==len(compare_subdirs)
        
    
        #print l.min(),l.max()
        #print seds.min(),seds.max()

    ## make separate plots for each snap, one subplot for each cam:
    for i_snap,snap in enumerate(snap_arr):
        if single_snap_plots:
            print "making subdir compare plots for snap %d"%snap

            plt.cla()
            plt.clf()
            plt.close()
            fig = plt.figure(figsize=(8,8))

            max_diff=0.0
            max_diff_clipped=0.0
            max_diff_oir=0.0
            max_diff_ir=0.0
            lclip=1.0e-7 if lambda_units=='m' else 1.0e3
            loir=3.0e-7 if lambda_units=='m' else 3.0e3
            lir=1.0e-6 if lambda_units=='m' else 1.0e4
            if lambda_units=='A':
                lclip,loir,lir = lclip*1.0e10,loir*1.0e10,lir*1.0e10
            if lambda_units=='microns':
                lclip,loir,lir = lclip*1.0e6,loir*1.0e6,lir*1.0e6

            for i in range(ncam):
            #print 'making sed plot for cam %d.'%i

                ax = fig.add_subplot(331+i)
                if xlim[1]/(1.0*xlim[0])>100: ax.set_xscale('log')
                ax.set_xlim(xlim)
                if ylog: ax.set_yscale('log')
                ax.set_ylim(ylim)
                ax.set_xlabel(r'$\lambda$ [%s]'%lu_lbl)
                ax.set_ylabel(ylbl)
                ax.set_title('cam %d'%i,fontsize=10)
            #plt.ylim(np.log10(ylim))

                if i==0: lh=None

                if oplot_wise: 
                    for k in range(4):  
                        ax.fill_between([wcens[k]-0.5*wweff[k],wcens[k]+0.5*wweff[k]],
                                        [ylim[0]]*2,[ylim[1]]*2,color='lightgray')
                        print wcens[k],0.5*ylim[1],0.9*ylim[1],10**(np.log10(ylim[1]-0.3))
                        print 'W%d'%(k+1)
                        ax.text(wcens[k],10**(np.log10(ylim[1])-0.3),'W%d'%(k+1),fontsize=10)
                    #ax.plot([wcens[k]]*2,ylim,'k-')
                for j in range(len(compare_subdirs)):
                    ns_s=ns_seds[:,i,j,i_snap]
                    scat_s=scat_seds[:,i,j,i_snap]
                    ir_s=ir_seds[:,i,j,i_snap]
                #print np.nanmin(s),np.nanmax(s)
                #print l.min(),l.max()
                    if flux_units=='relative': 
                        ns_s=ns_s/ns_seds[:,i,0,i_snap]
                        scat_s=scat_s/scat_seds[:,i,0,i_snap]
                        ir_s=ir_s/ir_seds[:,i,0,i_snap]
                    #max_diff = np.maximum(max_diff,np.maximum(1.0-np.nanmin(s),
                    #                                          np.nanmax(s)-1.0))
                    #max_diff_clipped = np.maximum(max_diff_clipped,
                    #                              np.maximum(1.0-np.nanmin(s[l>lclip]),
                    #                                         np.nanmax(s[l>lclip])-1.0))
                    #max_diff_oir = np.maximum(max_diff_oir,
                    #                          np.maximum(1.0-np.nanmin(s[l>loir]),
                    #                                     np.nanmax(s[l>loir])-1.0))
                    #max_diff_ir = np.maximum(max_diff_ir,
                    #                          np.maximum(1.0-np.nanmin(s[l>lir]),
                    #                                     np.nanmax(s[l>lir])-1.0))
                    #if max_diff_clipped>0.02:
                    #    print s[:10]
                #print np.nanmin(seds[:,i,0,i_snap]),np.nanmax(seds[:,i,0,i_snap])
                    elif flux_units=='nuFnu':
                        ns_s=l_m*ns_s
                        scat_s=l_m*scat_s
                        ir_s=l_m*ir_s
                    print ns_s.min(), ns_s.max()
                    print scat_s.min(), scat_s.max()
                    print ir_s.min(), ir_s.max()
                    if sed_type=='all':
                        p0,=ax.plot(l,ns_s)
                        p1,=ax.plot(l,scat_s)
                        p2,=ax.plot(l,ir_s)
                        p3,=ax.plot(l,scat_s+ir_s)
                    elif sed_type == 'tot':
                        p0,=ax.plot(l,scat_s+ir_s)
                    elif sed_type == 'nonscatter':
                        p0,=ax.plot(l,ns_s)
                    elif sed_type == 'scatter':
                        p0,=ax.plot(l,scat_s)
                    else:
                        p0,=ax.plot(l,ir_s)
                #p,=ax.plot(l,s)
                    if i==0:
                    #lh=(p,) if lh==None else lh+(p,)
                        if sed_type=='all':
                            lh=(p0,p1,p2,p3)
                            ll=('nonscatter','scatter','ir','tot')
                        else:
                            lh=(p0,) if lh==None else lh+(p0,)
                            ll=compare_subdirs
                if i==0: ax.legend(lh,ll,fontsize=9,loc='upper right',handletextpad=0.08)
            
        #if flux_units=='relative':
        #    print '\nmax_diff: %g'%max_diff
        #    print 'max_diff_clipped: %g'%max_diff_clipped
        #    print 'max_diff_oir: %g'%max_diff_oir
        #    print 'max_diff_ir: %g\n'%max_diff_ir
            fig.subplots_adjust(bottom=0.08,left=0.1,wspace=0.45,hspace=0.4,top=0.95,right=0.97)
            fig.savefig('%s/%s/%s/sed_compare%s_%s_%.3d%s.pdf'%(maindir,rundir,compare_subdirs[0],
                                                                compare_fbase,sed_type,snap,extra))



    ## make separate plots for each sim, one subplot for each snap:
    if len(snap_arr)>10: 
        print "WARNING: truncating number of snaps for snap-comparison plot."
        snap_arr = snap_arr[:10]

    cams=[single_cam] if single_cam>=0 else range(ncam)

    #for i_snap,snap in enumerate(snap_arr):

    if len(compare_subdirs)>1 and len(cams)>1:
        print "***WARNING: plotting multiple subdirs AND multiple cams in same snap panel!***"
    
    plt.cla()
    plt.clf()
    plt.close()
    fig = plt.figure(figsize=(5,8))

    max_diff=0.0
    max_diff_clipped=0.0
    max_diff_oir=0.0
    max_diff_ir=0.0
    lclip=1.0e-7 if lambda_units=='m' else 1.0e3
    loir=3.0e-7 if lambda_units=='m' else 3.0e3
    lir=1.0e-6 if lambda_units=='m' else 1.0e4
    if lambda_units=='A':
        lclip,loir,lir = lclip*1.0e10,loir*1.0e10,lir*1.0e10
    if lambda_units=='microns':
        lclip,loir,lir = lclip*1.0e6,loir*1.0e6,lir*1.0e6
    lsty=('-',':',':-','--',np.repeat('-',6))
    color=('g','k','r','b',np.repeat('k',6))

    matplotlib.rcParams.update({'font.size': 12})
    for i_snap,snap in enumerate(snap_arr):
        print "plotting snap %d"%snap
            #print 'making sed plot for cam %d.'%i

        ax = fig.add_subplot(len(snap_arr)*100+11+i_snap)
        if xlim[1]/(1.0*xlim[0])>100: ax.set_xscale('log')
        ax.set_xlim(xlim)
        if ylog: ax.set_yscale('log')
        ax.set_ylim(ylim)
        ax.set_ylabel(ylbl)
        if snap==snap_arr[-1]:
            ax.set_xlabel(r'$\lambda$ [%s]'%lu_lbl)
        else:
            ax.set_xticklabels(np.repeat('',10))            
        if snap_titles:
            ax.set_title('snap %d'%snap,fontsize=10)
            #plt.ylim(np.log10(ylim))
        if snap_labels: 
            ax.text(xlim[0]+0.65*(xlim[1]-xlim[0]),
                    ylim[0]+0.3*(ylim[1]-ylim[0]),snap_labels[i_snap],fontsize=13)

        if i_snap==0: lh=None
            
        if oplot_wise: 
            for k in range(4):  
                ax.fill_between([wcens[k]-0.5*wweff[k],wcens[k]+0.5*wweff[k]],
                                [ylim[0]]*2,[ylim[1]]*2,color='lightgray')
                    #ax.plot([wcens[k]]*2,ylim,'k-')
            if i_snap==0:
                ax.text(wcens[0]-0.73*wweff[0],10**(np.log10(ylim[1])-0.5),'W1',fontsize=10)
                ax.text(wcens[1]-0.51*wweff[1],10**(np.log10(ylim[1])-0.5),'W2',fontsize=10)
                ax.text(wcens[2]-0.435*wweff[2],10**(np.log10(ylim[1])-0.5),'W3',fontsize=10)
                ax.text(wcens[3]-0.73*wweff[3],10**(np.log10(ylim[1])-0.5),'W4',fontsize=10)
        for j in range(len(compare_subdirs)):
            kwargs = dict(ls=lsty[j],lw=1.9,c=color[j])
            print "plotting subdir %s"%compare_subdirs[j]
            for i in cams:
            #for j in range(len(compare_subdirs)):
                ns_s=ns_seds[:,i,j,i_snap]
                scat_s=scat_seds[:,i,j,i_snap]
                ir_s=ir_seds[:,i,j,i_snap]
                #print np.nanmin(s),np.nanmax(s)
                #print l.min(),l.max()
                if flux_units=='relative': 
                    ns_s=ns_s/ns_seds[:,i,0,i_snap]
                    scat_s=scat_s/scat_seds[:,i,0,i_snap]
                    ir_s=ir_s/ir_seds[:,i,0,i_snap]
                    #max_diff = np.maximum(max_diff,np.maximum(1.0-np.nanmin(s),
                    #                                          np.nanmax(s)-1.0))
                    #max_diff_clipped = np.maximum(max_diff_clipped,
                    #                              np.maximum(1.0-np.nanmin(s[l>lclip]),
                    #                                         np.nanmax(s[l>lclip])-1.0))
                    #max_diff_oir = np.maximum(max_diff_oir,
                    #                          np.maximum(1.0-np.nanmin(s[l>loir]),
                    #                                     np.nanmax(s[l>loir])-1.0))
                    #max_diff_ir = np.maximum(max_diff_ir,
                    #                         np.maximum(1.0-np.nanmin(s[l>lir]),
                    #                                    np.nanmax(s[l>lir])-1.0))
                    #if max_diff_clipped>0.02:
                    #    print s[:10]
                #print np.nanmin(seds[:,i,0,i_snap]),np.nanmax(seds[:,i,0,i_snap])
                elif flux_units=='nuFnu':
                    ns_s=l_m*ns_s
                    scat_s=l_m*scat_s
                    ir_s=l_m*ir_s
                print "ns, scat, ir min/max for subdir %s"%compare_subdirs[j]
                print ns_s.min(), ns_s.max()
                print scat_s.min(), scat_s.max()
                print ir_s.min(), ir_s.max()
                if sed_type=='all':
                    p0,=ax.plot(l,ns_s,**kwargs)
                    p1,=ax.plot(l,scat_s,**kwargs)
                    p2,=ax.plot(l,ir_s,**kwargs)
                    p3,=ax.plot(l,scat_s+ir_s,**kwargs)
                elif sed_type == 'tot':
                    p0,=ax.plot(l,scat_s+ir_s,**kwargs)
                elif sed_type == 'nonscatter':
                    p0,=ax.plot(l,ns_s,**kwargs)
                elif sed_type == 'scatter':
                    p0,=ax.plot(l,scat_s,**kwargs)
                else:
                    p0,=ax.plot(l,ir_s,**kwargs)
                #p,=ax.plot(l,s)
            #if i_snap==0:
            #       lh=(p,) if lh==None else lh+(p,)
            #       ll='%d'%snap
            if i_snap==0:
                    #lh=(p,) if lh==None else lh+(p,)
                if sed_type=='all':
                    lh=(p0,p1,p2,p3)
                    ll=('nonscatter','scatter','ir','tot')
                else:
                    lh=(p0,) if lh==None else lh+(p0,)
                    if compare_subdirs == ['all','all_agnx0']:
                        ll=('fiducial','AGNx0')
                    else:
                        ll=compare_subdirs
        ax.legend(lh,ll,fontsize=11,loc='upper left',handletextpad=0.08)
        #ax.legend(lh,snap_arr,fontsize=9,loc='upper right',handletextpad=0.08)
            
        #if flux_units=='relative':
        #    print '\nmax_diff: %g'%max_diff
        #    print 'max_diff_clipped: %g'%max_diff_clipped
        #    print 'max_diff_oir: %g'%max_diff_oir
        #    print 'max_diff_ir: %g\n'%max_diff_ir
    #fig.subplots_adjust(bottom=0.06,left=0.18,hspace=0.44,top=0.96,right=0.94)
    fig.subplots_adjust(bottom=0.06,left=0.18,hspace=0.05,top=0.96,right=0.94)

    if single_cam>0:
        fig.savefig('%s/%s/%s/sed_compare_snaps%s_%s_cam%d%s.pdf'%(maindir,rundir,compare_subdirs[0],
                                                                   compare_fbase,sed_type,cams[0],extra))
    else:
        fig.savefig('%s/%s/%s/sed_compare_snaps%s_%s_allcam%s.pdf'%(maindir,rundir,compare_subdirs[0],
                                                                    compare_fbase,sed_type,extra))
