import pyfits, glob, re, sys
import numpy as np
import broadband as bb
import astro_constants as ac
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from copy import copy

## Jarrett et al. 2011
W1_Vega_to_AB = 2.699
W2_Vega_to_AB = 3.339
W3_Vega_to_AB = 5.174
W4_Vega_to_AB = 6.620

def wise_colors(path='/oasis/projects/nsf/hvd115/lblecha',
                subdir='q0.5_fg0.3_sunruns/new_fiducial_all',
                snapnum=[],ncam=7,skip_snap0=True,input_type='bb',
                write_to_txt=False,write_restframe_lums=False,write_llambda=False):

    global W1_Vega_to_AB, W2_Vega_to_AB, W3_Vega_to_AB, W4_Vega_to_AB

    assert input_type in ['bb','mcrx']

    if write_restframe_lums and not input_type=='mcrx':
        print "WARNING! Keyword write_restframe_lums requires input_type=='mcrx'. Resetting input type..."
        input_type='mcrx'

    if write_llambda and input_type=='mcrx':
        print "WARNING! Keyword write_llambda requires input_type=='bb'. Resetting input type..."
        input_type='bb'
        if write_restframe_lums: print "WARNING! cannot set both write_restframe_lums and write_llambda. setting write_restframe_lums=False."

    fpath=path+'/'+subdir

    if np.isscalar(snapnum): snapnum=[snapnum]
    snapnum = np.array(snapnum)
    bbfiles = glob.glob(fpath+'/broadband*fits')   
    mcfiles = glob.glob(fpath+'/mcrx*fits')   
    tmp = str.split(str.split(str.split(bbfiles[0],'/')[-1],'.fits')[0],'_')
    bbfbase = '_'.join(tmp[:-1])
    if len(snapnum)==0:
        snapnum = np.sort(np.array([np.int(f.split('/')[-1].split(bbfbase+'_')[-1].split('.fits')[0])
                                    for f in bbfiles]))
        #snapnum = np.sort(np.array([np.int(re.search('(?<=broadband_)\d+',f).group(0))
        #                          for f in files]))
        assert len(snapnum)>0

    if skip_snap0: snapnum = snapnum[snapnum>0]

    print "\n%s:"%fpath
    print "calculating W1W2 for the following snapshot(s):",snapnum

    W1Leff_arr = np.nan*np.ones((ncam,len(snapnum)))
    W2Leff_arr = np.nan*np.ones((ncam,len(snapnum)))
    W3Leff_arr = np.nan*np.ones((ncam,len(snapnum)))
    W4Leff_arr = np.nan*np.ones((ncam,len(snapnum)))
    W1Llambda_arr = np.nan*np.ones((ncam,len(snapnum)))
    W2Llambda_arr = np.nan*np.ones((ncam,len(snapnum)))
    W3Llambda_arr = np.nan*np.ones((ncam,len(snapnum)))
    W4Llambda_arr = np.nan*np.ones((ncam,len(snapnum)))
    W1W2_arr = np.nan*np.ones((ncam,len(snapnum)))
    W2W3_arr = np.nan*np.ones((ncam,len(snapnum)))
    time = np.nan*np.ones(len(snapnum))

    for n in range(len(snapnum)):

        fmc = pyfits.open('%s/mcrx_%.3d.fits'%(fpath,snapnum[n]))
        time[n] = fmc['SFRHIST'].header['SNAPTIME']/1.0e9
        print "snap %d: t=%.3g Gyr."%(snapnum[n],time[n])

        #fbb = pyfits.open('%s/broadband_%.3d.fits'%(fpath,snapnum[n]))
        fbb = pyfits.open('%s/%s_%.3d.fits'%(fpath,bbfbase,snapnum[n]))
        
        filterdata = fbb['FILTERS'].data
        filternames = np.array([ x.strip() for x in filterdata['filter'] ])
        wise_filternames = np.array([ 'WISE-W%d.res'%i if 'WISE-W%d.res'%i in filternames 
                                          else 'wise/WISE-W%d.res'%i for i in np.arange(1,5) ])
        if input_type=='bb':
            #ixw1 = filternames.tolist().index('wise/WISE-W1.res')
            ixw1 = filternames.tolist().index(wise_filternames[0])
            lambda_eff_w1 = filterdata['lambda_eff'][ixw1]
            ixw2 = filternames.tolist().index(wise_filternames[1])
            lambda_eff_w2 = filterdata['lambda_eff'][ixw2]
            ixw3 = filternames.tolist().index(wise_filternames[2])
            lambda_eff_w3 = filterdata['lambda_eff'][ixw3]
            ixw4 = filternames.tolist().index(wise_filternames[3])
            lambda_eff_w4 = filterdata['lambda_eff'][ixw4]
            
            #using 'scatter' values only
            w1Llam = np.array([filterdata['L_lambda_eff%d'%i][ixw1] for i in range(7)])
            w1magAB = np.array([filterdata['AB_mag%d'%i][ixw1] for i in range(7)])

            w2Llam = np.array([filterdata['L_lambda_eff%d'%i][ixw2] for i in range(7)])
            w2magAB = np.array([filterdata['AB_mag%d'%i][ixw2] for i in range(7)])

            w3Llam = np.array([filterdata['L_lambda_eff%d'%i][ixw3] for i in range(7)])
            w3magAB = np.array([filterdata['AB_mag%d'%i][ixw3] for i in range(7)])
                        
            w4Llam = np.array([filterdata['L_lambda_eff%d'%i][ixw4] for i in range(7)])
            w4magAB = np.array([filterdata['AB_mag%d'%i][ixw4] for i in range(7)])
                        

        elif input_type=='mcrx':
            #lbol_tot = fmc[48].header['L_bol_grid']
            #lbol_cam = np.array([fmc[48].header['L_SCAT%d'%i] for i in range(7)])
            fmc = pyfits.open('%s/%s_%.3d.fits'%(fpath,'mcrx',snapnum[n]))

            lam = fmc[5].data.field('lambda')
            llam_cam = np.array([fmc[48].data.field('L_lambda_out%d'%i) 
                                 for i in range(7)])
            w1magAB = np.array([ bb.ABMag(wise_filternames[0].strip('.res'),lam,llam_cam[i],sed_lam_units='m')
                                 for i in range(7) ])
            w2magAB = np.array([ bb.ABMag(wise_filternames[1].strip('.res'),lam,llam_cam[i],sed_lam_units='m')
                                 for i in range(7) ])
            w3magAB = np.array([ bb.ABMag(wise_filternames[2].strip('.res'),lam,llam_cam[i],sed_lam_units='m')
                                 for i in range(7) ])

            ### band luminosity in W ###
            w1Leff=np.array([bb.Leff(wise_filternames[0].strip('.res'),lam,llam_cam[i],sed_lam_units='m')
                             for i in range(7) ])
            w2Leff=np.array([bb.Leff(wise_filternames[1].strip('.res'),lam,llam_cam[i],sed_lam_units='m')
                             for i in range(7) ])
            w3Leff=np.array([bb.Leff(wise_filternames[2].strip('.res'),lam,llam_cam[i],sed_lam_units='m')
                             for i in range(7) ])
            w4Leff=np.array([bb.Leff(wise_filternames[3].strip('.res'),lam,llam_cam[i],sed_lam_units='m')
                             for i in range(7) ])
            fmc.close()

        fbb.close()

        w1magVega = w1magAB - W1_Vega_to_AB
        w2magVega = w2magAB - W2_Vega_to_AB
        w3magVega = w3magAB - W3_Vega_to_AB
        
        w1minusw2_Vega = w1magVega - w2magVega
        w2minusw3_Vega = w2magVega - w3magVega
            
        W1W2_arr[:,n] = w1minusw2_Vega
        W2W3_arr[:,n] = w2minusw3_Vega
        if write_restframe_lums:
            W1Leff_arr[:,n] = w1Leff
            W2Leff_arr[:,n] = w2Leff
            W3Leff_arr[:,n] = w3Leff
            W4Leff_arr[:,n] = w4Leff
        if write_llambda:
            W1Llambda_arr[:,n] = w1Llam
            W2Llambda_arr[:,n] = w2Llam
            W3Llambda_arr[:,n] = w3Llam
            W4Llambda_arr[:,n] = w4Llam

    if write_restframe_lums:
        for band in np.arange(1,5):
            fpl = open(fpath+'/w%d_lums_restframe.txt'%band,'w')
            if band==1: lum_arr = copy(W1Leff_arr)
            if band==2: lum_arr = copy(W2Leff_arr)
            if band==3: lum_arr = copy(W3Leff_arr)
            if band==4: lum_arr = copy(W4Leff_arr)
            for n in np.arange(len(time)):
                pstr=' '.join([str(l) for l in lum_arr[:,n]])
                fpl.write("%g %g %s\n"%(snapnum[n],time[n],pstr))
            fpl.close()

    if write_llambda:
        for band in np.arange(1,5):
            fpll = open(fpath+'/W%d_Llambda.txt'%band,'w')
            if band==1: ll_arr = copy(W1Llambda_arr)
            if band==2: ll_arr = copy(W2Llambda_arr)
            if band==3: ll_arr = copy(W3Llambda_arr)
            if band==4: ll_arr = copy(W4Llambda_arr)
            for n in np.arange(len(time)):
                pstr=' '.join([str(l) for l in ll_arr[:,n]])
                fpll.write("%g %g %s\n"%(snapnum[n],time[n],pstr))
            fpll.close()

    if write_to_txt:
        fp12 = open(fpath+'/w1w2_color.txt','w')
        fp23 = open(fpath+'/w2w3_color.txt','w')
        for n in np.arange(len(time)):
            pstr12 = "%g %g"%(snapnum[n],time[n])
            pstr23 = "%g %g"%(snapnum[n],time[n])
            for i in range(ncam): 
                pstr12 = pstr12+" %g"%W1W2_arr[i,n]
                pstr23 = pstr23+" %g"%W2W3_arr[i,n]
            fp12.write(pstr12+"\n")
            fp23.write(pstr23+"\n")
        fp12.close()
        fp23.close()

    return snapnum, time, W1W2_arr,W2W3_arr


def wise_colors_bhfile(path='/home/lblecha/sundata/', fname='chris_newbhmodel.fits',
                       write_to_txt=False):

    global W1_Vega_to_AB, W2_Vega_to_AB, W3_Vega_to_AB, W4_Vega_to_AB

    fbase=fname.strip('.fits')
    print "calculating W1W2 for the following bhmodel file: %s/%s"%(path,fname)

    f = pyfits.open('%s/%s'%(path,fname))
    sed_arr = f['SED'].data
    lglbol = f['AXES'].data['log_L_bol']
    lglbol = lglbol[:sed_arr.shape[1]]
    lam = f['AXES'].data['lambda']
    
    w1magAB = np.array([ bb.ABMag('wise/WISE-W1',lam,sed_arr[:,n],sed_lam_units='m')
                         for n in range(len(lglbol)) ])
    w2magAB = np.array([ bb.ABMag('wise/WISE-W2',lam,sed_arr[:,n],sed_lam_units='m') 
                         for n in range(len(lglbol)) ])
    w3magAB = np.array([ bb.ABMag('wise/WISE-W3',lam,sed_arr[:,n],sed_lam_units='m')
                         for n in range(len(lglbol)) ])

    w1magVega = w1magAB - W1_Vega_to_AB
    w2magVega = w2magAB - W2_Vega_to_AB
    w3magVega = w3magAB - W3_Vega_to_AB
        
    W1W2 = w1magVega - w2magVega
    W2W3 = w2magVega - w3magVega
    
    if write_to_txt:
        fp12 = open(path+'/%s_w1w2_color.txt'%fbase,'w')
        fp23 = open(path+'/%s_w2w3_color.txt'%fbase,'w')
        for n in np.arange(len(lglbol)):
            fp12.write("%g %g\n"%(lglbol[n],W1W2[n]))
            fp23.write("%g %g\n"%(lglbol[n],W2W3[n]))
        fp12.close()
        fp23.close()

    return lglbol, W1W2, W2W3



def calc_Av(path='/oasis/projects/nsf/hvd115/lblecha',
            subdir='q0.5_fg0.3_sunruns/test_fiducial',
            snapnum=[],skip_snap0=True,ncam=7,
            write_to_txt=False):

    if np.isscalar(snapnum): snapnum=[snapnum]
    if len(snapnum)==0:
        bbfiles = glob.glob(path+'/'+subdir+'/broadband*fits')   
        tmp = str.split(str.split(str.split(bbfiles[0],'/')[-1],'.fits')[0],'_')
        fbase = '_'.join(tmp[:-1])
        snapnum = np.sort(np.array([np.int(f.split('/')[-1].split(fbase+'_')[-1].split('.fits')[0])
                                    for f in bbfiles]))
        #snapnum = np.sort(np.array([np.int(re.search('(?<=broadband_)\d+',f).group(0))
        #                          for f in bbfiles]))
        assert len(snapnum)>0
    if skip_snap0: snapnum = snapnum[snapnum>0]

    print "\n%s/%s:"%(path,subdir)
    print "calculating Av for the following snapshot(s):",snapnum

    Av_arr = np.nan*np.ones((ncam,len(snapnum)))
    time = np.nan*np.ones(len(snapnum))
    for n in range(len(snapnum)):
        
        fmc = pyfits.open('%s/%s/mcrx_%.3d.fits'%(path,subdir,snapnum[n]))
        time[n] = fmc['SFRHIST'].header['SNAPTIME']/1.0e9
        fmc.close()
        print "snap %d: t=%.3g Gyr."%(snapnum[n],time[n])
        fbb = pyfits.open('%s/%s/%s_%.3d.fits'%(path,subdir,fbase,snapnum[n]))
        filterdata = fbb['FILTERS'].data
        filternames = np.array([ x.strip() for x in filterdata['filter'] ])
        ixV = filternames.tolist().index('johnson/V_Johnson.res')
        lambda_eff_V = filterdata['lambda_eff'][ixV]
        
        #'scatter' values 
        Vlum = np.array([filterdata['L_lambda_eff%d'%i][ixV] for i in range(7)])
        VmagAB = np.array([filterdata['AB_mag%d'%i][ixV] for i in range(7)])
        #'nonscatter' values only
        Vlum_nonscat = np.array([filterdata['L_lambda_eff_nonscatter%d'%i][ixV] for i in range(7)])
        VmagAB_nonscat = np.array([filterdata['AB_mag_nonscatter%d'%i][ixV] for i in range(7)])

        fbb.close()

        Av_arr[:,n] = VmagAB - VmagAB_nonscat
        #Av_arr[:,n] = 2.5 * np.log10(Vlum_nonscat/Vlum)
        #print "VmagAB:",VmagAB
        #print "VmagAB_nonscat:",VmagAB_nonscat

    #print Av_arr
    if write_to_txt:
        fp = open(path+'/'+subdir+'/Av.txt','w')
        for n in np.arange(len(time)):
            pstr = "%g %g"%(snapnum[n],time[n])
            for i in range(ncam): 
                pstr = pstr+" %g"%Av_arr[i,n]
            fp.write(pstr+"\n")

    return time, Av_arr
        
def make_plots(path='/oasis/projects/nsf/hvd115/lblecha',
               subdir='q0.5_fg0.3_sunruns/test_fiducial',
               snapnum=[],skip_snap0=True,ncam=7,mycalc=False,
               write_to_txt=False):
    
    if np.isscalar(snapnum): snapnum=[snapnum]
    if len(snapnum)==0:
        bbfiles = glob.glob(path+'/'+subdir+'/broadband*fits')   
        snapnum = np.sort(np.array([np.int(re.search('(?<=broadband_)\d+',f).group(0))
                                  for f in bbfiles]))
        assert len(snapnum)>0
    if skip_snap0: snapnum = snapnum[snapnum>0]

    w1w2_file = glob.glob(path+'/'+subdir+'/w1w2_color.txt')
    w2w3_file = glob.glob(path+'/'+subdir+'/w2w3_color.txt')
    Av_file = glob.glob(path+'/'+subdir+'/Av.txt')

    if not write_to_txt and len(w1w2_file)==1 and len(w2w3_file)==1:
        tmpdata = np.loadtxt(w1w2_file[0],unpack=True)
        print tmpdata.shape
        snapnum = tmpdata[0,:]
        time = tmpdata[1,:]
        w1w2 = tmpdata[2:,:]
        tmpdata = np.loadtxt(w2w3_file[0],unpack=True)
        print tmpdata.shape
        snapnum_23 = tmpdata[0,:]
        time_23 = tmpdata[1,:]
        w2w3 = tmpdata[2:,:]
        if len(snapnum) != len(snapnum_23) or len(time) != len(time_23):
            print "Error: snapnum or time arrays for w1w2 and w2w3 do not match."
            print len(snapnum),len(time)
            print len(snapnum_23),len(time_23)
            return
    else:
        print "Could not read WISE colors from file. Calculating..."
        snapnum,time,w1w2,w2w3 = wise_colors(path=path,subdir=subdir,snapnum=snapnum,
                                             ncam=ncam,skip_snap0=skip_snap0,
                                             write_to_txt=write_to_txt)
    print time.shape
    print w1w2.shape

    if not write_to_txt and len(Av_file)==1:
        tmpdata = np.loadtxt(Av_file[0],unpack=True)
        print tmpdata.shape
        snapnum_Av = tmpdata[0,:]
        time_Av = tmpdata[1,:]
        Av = tmpdata[2:,:]
        if len(snapnum) != len(snapnum_Av) or len(time) != len(time_Av):
            print "Error: snapnum or time arrays for w1w2 and Av do not match."
            print len(snapnum),len(time)
            print len(snapnum_Av),len(time_Av)
            return
    else:
        print "Could not read Av from file. Calculating..."
        time_Av,Av = calc_Av(path=path,subdir=subdir,snapnum=snapnum,ncam=ncam,
                             skip_snap0=skip_snap0,write_to_txt=write_to_txt)

    print "min/max w1w2:",w1w2.min(),w1w2.max()
    print "min/max w2w3:",w2w3.min(),w2w3.max()
    print "min/max Av:",Av.min(),Av.max()

    plot_LOSerr_vs_t(path,subdir,'Av',time,np.mean(Av,axis=0),
                     np.min(Av,axis=0),np.max(Av,axis=0),ylabel='Av')
    plot_LOSerr_vs_t(path,subdir,'W1W2',time,np.mean(w1w2,axis=0),
                     np.min(w1w2,axis=0),np.max(w1w2,axis=0),ylabel='W1-W2')
    plot_LOSerr_vs_t(path,subdir,'W2W3',time,np.mean(w2w3,axis=0),
                     np.min(w2w3,axis=0),np.max(w2w3,axis=0),ylabel='W2-W3')


def plot_LOSerr_vs_t(path,subdir,plotname,time,meanval,minval,maxval,
                     meanval2=(),minval2=(),maxval2=(),
                     lcolor='k',ecolor='gray',lcolor2='b',ecolor2='c',
                     oplot_times=(),oplot_vals=(),oplot_vals_style=(),
                     xlim=(),ylim=(),ylabel='',plt_fmt='pdf'):

    if any(len(lst)!=len(time) for lst in (meanval,minval,maxval)):
        print "Error: array size mismatch."
        print time.size,meanval.size,minval.size,maxval.size
        return

    if len(meanval2)>0 and any(len(lst)!=len(time) for lst in (meanval2,minval2,maxval2)):
        print "Error: array size mismatch for secondary values."
        print time.size,meanval2.size,minval2.size,maxval2.size
        return        
        
    if plt_fmt not in ('pdf','eps'):
        print "WARNING: unrecognized plt_fmt: %s. Defaulting to pdf."%plt_fmt
        plt_fmt='pdf'

    if len(xlim)!=2: xlim = (0.95*time.min(),1.05*time.max())
    if len(ylim)!=2: ylim = (0.95*minval.min(),1.05*maxval.max())

    plt.clf()
    plt.close()
    fig = plt.figure(figsize=(5,3))
    
    plt.xlabel('time [Gyr]')
    plt.ylabel(ylabel)
    plt.xlim(xlim)
    plt.ylim(ylim)

    if len(oplot_times)>0:
        for t in oplot_times: plt.plot([t,t],ylim,'k',ls='--')

    plt.errorbar(time, meanval, yerr=(meanval-minval,maxval-meanval),color=ecolor)
    plt.plot(time, meanval, lcolor, linewidth=1.5) 

    if len(meanval2)>0:
        plt.errorbar(time, meanval2, yerr=(meanval2-minval2,maxval2-meanval2),color=ecolor2)
        plt.plot(time, meanval2, lcolor2, linewidth=1.5) 


    if len(oplot_vals)>0:                                        
        if len(oplot_vals_style)!=len(oplot_vals):               
            oplot_vals_style=('r',)*len(oplot_vals)
        for i,v in enumerate(oplot_vals):
            plt.plot(xlim,[v,v],oplot_vals_style[i])

    #plt.title(plotname)
    #fig.subplots_adjust(hspace=0.5,wspace=0.3)
    plt.tight_layout()
    plt.savefig('%s/%s/%s.%s'%(path,subdir,plotname,plt_fmt),format=plt_fmt)

    
