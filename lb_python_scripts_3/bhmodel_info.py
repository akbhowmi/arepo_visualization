import numpy as np
import pyfits, glob, re, sys
import numpy as np
import broadband as bb
import astro_constants as ac
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt


## Jarrett et al. 2011
W1_Vega_to_AB = 2.699
W2_Vega_to_AB = 3.339
W3_Vega_to_AB = 5.174
W4_Vega_to_AB = 6.620


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
        
    f.close()
    return lglbol, W1W2, W2W3


def plot_seds(path='/home/lblecha/sundata/',fname='chris_newbhmodel.fits',
              relative=False):

    fbase=fname.strip('.fits')
    print "plotting seds for the following bhmodel file: %s/%s"%(path,fname)

    f = pyfits.open('%s/%s'%(path,fname))
    sed_arr = f['SED'].data
    lglbol = f['AXES'].data['log_L_bol']
    lglbol = lglbol[:sed_arr.shape[1]]
    lam = f['AXES'].data['lambda']
    fac=10**lglbol[0]/10**lglbol
    #fac=(10**lglbol[0]/10**lglbol[i]) if relative else 1.0

    f.close()

    plt.close()
    plt.clf()
    fig = plt.figure()
    #plt.xlim(1e-12,1e-3)
    #plt.plot([3.4e-6,3.4e-6],[30,55],'k:')
    #plt.plot([4.5e-6,4.5e-6],[30,55],'k:')

    ax1=fig.add_subplot(211)
    for i in range(0,len(lglbol),10):
        ax1.set_xlim(5e-13,1e-3)
        ax1.set_ylim(3e30,3e42)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r'$\lambda$ [m]')
        ax1.set_ylabel(r'$\lambda$ L$_{\lambda}$ [W]')
        ax1.plot(lam,lam*10**sed_arr[:,i],lw=1+i*0.025)
        plt.title('min/max Lbol = %.3g,%.3g [log W]'%(lglbol.min(),lglbol.max()))

    ax2=fig.add_subplot(212)
    for i in range(0,len(lglbol),10):
        print "lglbol[%d]=%g (fac=%g)"%(i,lglbol[i],fac[i])
        ax2.set_xlim(5e-13,1e-3)
        ax2.set_ylim(1e31,1e35)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel(r'$\lambda$ [m]')
        ax2.set_ylabel(r'$\lambda$ L$_{\lambda}$ [W]')
        ax2.plot(lam,lam*10**sed_arr[:,i]*fac[i],lw=1+i*0.025)
        plt.title('scaled to min Lbol')
        #plt.plot(lam,sed_arr[:,i] - np.log10(fac[i])

    fig.subplots_adjust(hspace=0.4)
    fig.show()
    fig.savefig('%s/%s_seds.pdf'%(path,fbase))

    return
